"""
Experimental single-pass giant-cube combiner using reproject's zarr return mode.

This is an alternative to the channel-by-channel `make_giant_mosaic_cube` path
in `make_mosaic.py`.  Instead of mosaicking 350 channels independently and then
stitching, it hands the full list of field cubes to a single
`reproject.mosaicking.reproject_and_coadd` call with ``return_type='zarr'``,
which computes the 3D mosaic block-by-block into a zarr store (bounded memory,
truly parallel I/O) and returns dask arrays reading from that store.  The store
is then streamed to FITS.

Requires:
  * reproject with the ``return_type='zarr'`` feature (astropy/reproject#621,
    branch ``coadd-zarr-return``; merged into main 2026-07-06).
  * These run in the dedicated ``python312-zarr`` conda env with editable
    worktrees of reproject / spectral-cube / reduction_ACES, so the production
    ``python312`` env (used by running jobs) is untouched.

Settings recommended by astrofrog (gist 81729066204b0332735b9461d015ee7c),
tuned for the ACES 12m giant cubes:

    parallel=16, block_size=(1, 2048, -1),
    return_type='zarr', zarr_path=..., zarr_batch_size=64

The cube/weight/beam loading mirrors `make_mosaic.make_giant_mosaic_cube`
(including the channel-independent flat/`min_weight_fraction` handling) so any
molecule that has a `make_giant_mosaic_cube_*` builder can be routed here by
passing the same filelist / reference_frequency / cdelt_kms / nchan.
"""
import os
import warnings

import numpy as np
from astropy import units as u
from astropy import wcs
from astropy.io import fits

from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from spectral_cube import SpectralCube

from aces.imaging.make_mosaic import (make_giant_mosaic_cube_header,
                                      get_common_beam,
                                      _weight_is_flat,
                                      _flat_weightcube)
from aces import conf

basepath = conf.basepath


def _load_cubes_and_weights(filelist, weightfilelist, reference_frequency,
                            beam_threshold, min_weight_fraction=0.05,
                            use_beams=True, verbose=True):
    """Load data + weight cubes exactly as make_giant_mosaic_cube does: spectral
    unit -> km/s at reference_frequency, drop cubes with beams above threshold,
    broadcast flat weights, and apply a CHANNEL-INDEPENDENT min_weight_fraction
    footprint mask.  Returns (cubes, weightcubes)."""
    reference_frequency = u.Quantity(reference_frequency, u.Hz)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cubes = [SpectralCube.read(fn,
                                   format='fits' if fn.endswith('fits') else 'casa_image',
                                   use_dask=True).with_spectral_unit(
                                       u.km / u.s, velocity_convention='radio',
                                       rest_value=reference_frequency)
                 for fn in filelist]

        weightcubes = []
        for fn, refcube in zip(weightfilelist, cubes):
            if _weight_is_flat(fn):
                weightcubes.append(_flat_weightcube(fn, refcube))
            else:
                weightcubes.append(
                    SpectralCube.read(fn, format='fits' if fn.endswith('fits') else 'casa_image',
                                      use_dask=True).with_spectral_unit(
                                          u.km / u.s, velocity_convention='radio',
                                          rest_value=reference_frequency))

        if min_weight_fraction is not None:
            # channel-independent footprint from the mid plane (see the matching
            # fix in make_mosaic.make_giant_mosaic_cube)
            masked = []
            for weightcube in weightcubes:
                midplane = weightcube[weightcube.shape[0] // 2, :, :]
                footprint = np.asarray(midplane > min_weight_fraction * midplane.max())
                masked.append(weightcube.with_mask(footprint[None, :, :]))
            weightcubes = masked

    # normalise timesys casing (some FITS write UTC in caps)
    for cube in cubes + weightcubes:
        cube._wcs.wcs.timesys = cube.wcs.wcs.timesys.lower()
        if hasattr(cube.mask, '_wcs'):
            cube.mask._wcs.wcs.timesys = cube.wcs.wcs.timesys.lower()

    # beam filter
    if use_beams:
        beams = [get_common_beam(cube.beams) if hasattr(cube, 'beams') else cube.beam
                 for cube in cubes]
        ok = [beam.major < beam_threshold for beam in beams]
        if verbose and not all(ok):
            dropped = [fn for fn, k in zip(filelist, ok) if not k]
            print(f"Beam filter dropped {len(dropped)} cubes > {beam_threshold}: {dropped}", flush=True)
        cubes = [c for k, c in zip(ok, cubes) if k]
        weightcubes = [c for k, c in zip(ok, weightcubes) if k]

    return cubes, weightcubes


def giant_mosaic_cube_zarr(filelist,
                           weightfilelist,
                           reference_frequency,
                           cdelt_kms,
                           cubename,
                           nchan,
                           output_file=None,
                           zarr_path=None,
                           target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr',
                           beam_threshold=2.75 * u.arcsec,
                           min_weight_fraction=0.05,
                           parallel=16,
                           block_size=(1, 2048, -1),
                           zarr_batch_size=64,
                           combine_function='mean',
                           use_beams=True,
                           overwrite=False,
                           verbose=True):
    """
    Build a giant mosaic cube in one reproject_and_coadd(return_type='zarr') call.

    Parameters mirror make_giant_mosaic_cube where they overlap.  ``zarr_path``
    must NOT already exist (reproject enforces this).  Writes the coadded array
    to ``output_file`` (FITS) by streaming from the zarr store.
    """
    if output_file is None:
        output_file = f'{basepath}/mosaics/cubes/{cubename}_CubeMosaic.fits'
    if zarr_path is None:
        zarr_path = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/{cubename}_CubeMosaic.zarr'

    if os.path.exists(zarr_path):
        if overwrite:
            import shutil
            print(f"Removing existing zarr store {zarr_path}", flush=True)
            shutil.rmtree(zarr_path)
        else:
            raise ValueError(f"zarr_path {zarr_path} already exists (pass overwrite=True to replace)")

    header = make_giant_mosaic_cube_header(target_header=target_header,
                                           reference_frequency=u.Quantity(reference_frequency, u.Hz),
                                           cdelt_kms=cdelt_kms,
                                           nchan=nchan)
    target_wcs = wcs.WCS(header)
    shape_out = (header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])

    if verbose:
        print(f"Loading {len(filelist)} cubes + weights for {cubename}", flush=True)
    cubes, weightcubes = _load_cubes_and_weights(
        filelist, weightfilelist, reference_frequency,
        beam_threshold=beam_threshold, min_weight_fraction=min_weight_fraction,
        use_beams=use_beams, verbose=verbose)

    if verbose:
        print(f"Coadding {len(cubes)} cubes into {shape_out} via return_type='zarr' "
              f"(parallel={parallel}, block_size={block_size}, batch={zarr_batch_size})", flush=True)

    # input_data: (array, wcs) tuples on each field's own 3D grid; reproject_interp
    # interpolates all three axes, so per-field spectral-grid differences are
    # handled here (no pre-regridding needed).
    input_data = [(cube.unitless_filled_data[:], cube.wcs) for cube in cubes]
    # input_weights: finite arrays, 0 outside the footprint mask
    input_weights = [np.nan_to_num(wc.unitless_filled_data[:]) for wc in weightcubes]

    array, footprint = reproject_and_coadd(
        input_data,
        target_wcs,
        shape_out=shape_out,
        input_weights=input_weights,
        reproject_function=reproject_interp,
        combine_function=combine_function,
        roundtrip_coords=False,
        parallel=parallel,
        block_size=block_size,
        return_type='zarr',
        zarr_path=zarr_path,
        zarr_batch_size=zarr_batch_size,
    )

    if verbose:
        print(f"Streaming coadded array {array.shape} -> {output_file}", flush=True)
    header['BUNIT'] = 'Jy/beam'
    fits.PrimaryHDU(data=array, header=header).writeto(output_file, overwrite=overwrite)
    if verbose:
        print(f"Wrote {output_file}", flush=True)
    return output_file

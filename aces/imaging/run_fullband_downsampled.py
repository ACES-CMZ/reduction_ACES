"""Driver for the full-band spatially-downsampled per-spw cubes.

Pick a spw with the SPW env var (25/27/33/35); optional overrides:
PIXSCALE_AS, BEAM_AS, PARALLEL, ZARR_BATCH_SIZE, BLOCK_SPATIAL, OVERWRITE.

    SPW=27 python -m aces.imaging.run_fullband_downsampled
"""
import os

from astropy import units as u

from aces.imaging.fullband_downsampled_cube import make_fullband_downsampled_cube, SPWS


def main():
    spw = int(os.environ['SPW'])
    if spw not in SPWS:
        raise ValueError(f"SPW={spw} not in {SPWS}")
    kwargs = dict(
        pixscale=float(os.environ.get('PIXSCALE_AS', 5)) * u.arcsec,
        target_beam=float(os.environ.get('BEAM_AS', 10)) * u.arcsec,
        parallel=int(os.environ.get('PARALLEL', 16)),
        block_size=(1, int(os.environ.get('BLOCK_SPATIAL', 2048)), -1),
        zarr_batch_size=int(os.environ.get('ZARR_BATCH_SIZE', 64)),
        overwrite=os.environ.get('OVERWRITE', 'False') == 'True',
    )
    print(f"[fullband] SPW{spw} pixscale={kwargs['pixscale']} beam={kwargs['target_beam']}", flush=True)
    make_fullband_downsampled_cube(spw, **kwargs)


if __name__ == '__main__':
    main()

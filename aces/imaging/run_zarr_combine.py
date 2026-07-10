"""
Driver for the experimental zarr giant-cube combiner (mosaic_cube_zarr).

Pick a cube with the MOLNAME env var; parameters (input glob, rest frequency,
channel width, nchan, spw, beam threshold) come from MOLCONFIG below, which
mirrors the corresponding make_giant_mosaic_cube_* builders in mosaic_12m.py.

    MOLNAME=CH3CHO_3m13 python -m aces.imaging.run_zarr_combine

Env overrides (all optional):
    ZARR_PATH, OUTPUT_FILE, PARALLEL, ZARR_BATCH_SIZE, BLOCK_SPATIAL (the middle
    block dim; block is (1, BLOCK_SPATIAL, -1)), OVERWRITE=True
"""
import glob
import os

from astropy import units as u

from aces.imaging.mosaic_12m import get_weightfile
from aces.imaging.mosaic_cube_zarr import giant_mosaic_cube_zarr

FEATHER = '/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP'


def _spw_glob(spw):
    return sorted(glob.glob(f'{FEATHER}/SPW{spw}/cubes/'
                            f'Sgr_A_st_*.TP_7M_12M_feather_all.SPW_{spw}.image.statcont.contsub.fits'))


# molname -> dict(spw, restfrq_Hz, cdelt_kms, nchan, beam_threshold_arcsec)
# values copied from the make_giant_mosaic_cube_* builders in mosaic_12m.py
MOLCONFIG = {
    'CH3CHO_3m13': dict(spw=35, restfrq=101.343448e9, cdelt_kms=1.47015502, nchan=350, beam_threshold=2.75),
    'HC3N': dict(spw=35, restfrq=100.0763e9, cdelt_kms=1.47015502, nchan=350, beam_threshold=2.75),
    'NSplus': dict(spw=35, restfrq=100.198e9, cdelt_kms=1.47015502, nchan=350, beam_threshold=2.75),
    'CH3OH_72-63': dict(spw=27, restfrq=86.902947e9, cdelt_kms=0.84455895, nchan=600, beam_threshold=3.3),
    'CS21': dict(spw=33, restfrq=97.98095330e9, cdelt_kms=1.4844932, nchan=350, beam_threshold=2.9),
}


def main():
    molname = os.environ['MOLNAME']
    if molname not in MOLCONFIG:
        raise KeyError(f"{molname} not in MOLCONFIG; add it (keys: {list(MOLCONFIG)})")
    cfg = MOLCONFIG[molname]

    filelist = _spw_glob(cfg['spw'])
    if len(filelist) == 0:
        raise FileNotFoundError(f"No SPW{cfg['spw']} feather cubes found for {molname}")
    weightfilelist = [get_weightfile(fn, spw=cfg['spw']) for fn in filelist]
    for wf in weightfilelist:
        assert os.path.exists(wf), f"missing weight {wf}"

    block_spatial = int(os.environ.get('BLOCK_SPATIAL', 2048))
    kwargs = dict(
        parallel=int(os.environ.get('PARALLEL', 16)),
        block_size=(1, block_spatial, -1),
        zarr_batch_size=int(os.environ.get('ZARR_BATCH_SIZE', 64)),
        beam_threshold=cfg['beam_threshold'] * u.arcsec,
        overwrite=os.environ.get('OVERWRITE', 'False') == 'True',
    )
    if os.environ.get('ZARR_PATH'):
        kwargs['zarr_path'] = os.environ['ZARR_PATH']
    if os.environ.get('OUTPUT_FILE'):
        kwargs['output_file'] = os.environ['OUTPUT_FILE']

    print(f"[zarr-combine] {molname}: {len(filelist)} cubes, spw{cfg['spw']}, "
          f"restfrq={cfg['restfrq']/1e9:.4f} GHz, nchan={cfg['nchan']}", flush=True)
    giant_mosaic_cube_zarr(filelist, weightfilelist,
                           reference_frequency=cfg['restfrq'],
                           cdelt_kms=cfg['cdelt_kms'],
                           cubename=molname,
                           nchan=cfg['nchan'],
                           **kwargs)


if __name__ == '__main__':
    main()

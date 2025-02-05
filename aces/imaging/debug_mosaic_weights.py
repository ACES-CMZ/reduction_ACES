import shutil
import glob
import radio_beam
import os
from astropy import units as u
from astropy.io import fits
from astropy import log
from astropy.wcs import WCS
from astropy.table import Table
from aces.imaging.make_mosaic import (read_as_2d, get_peak,
                                      get_m0,
                                      make_downsampled_cube,
                                      make_giant_mosaic_cube,
                                      rms
                                      )
from aces.imaging.make_mosaic import make_mosaic as make_mosaic_, all_lines as all_lines_
# import os
# from functools import partial
from multiprocessing import Process, Pool
import multiprocessing
import numpy as np
from spectral_cube.cube_utils import mosaic_cubes
import dask
from functools import partial
import multiprocessing.pool as mpp
from itertools import repeat

from aces import conf
import uvcombine


#if os.getenv('NO_PROGRESSBAR') is None and not (os.getenv('ENVIRON') == 'BATCH'):
#    from dask.diagnostics import ProgressBar
#    pbar = ProgressBar(minimum=20)  # don't show pbar <5s
#    pbar.register()

basepath = conf.basepath

from aces.imaging.mosaic_12m import get_weightfile

def main():
    kwargs = {}

    filelist = sorted(glob.glob(f'{basepath}/upload/Feather_12m_7m_TP/HNCO/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.hnco43.image.statcont.contsub.fits'))
    filelist = [x for x in filelist if any(y in x for y in ('_al.', '_ar.', '_z.', '_as.'))]
    # what if we reverse order?
    filelist = filelist[::-1]
    print(f"Found {len(filelist)} HNCO 7m+12m+TP FITS files")
    log.debug(f"Full names of files: {filelist}")

    weightfilelist = [get_weightfile(fn, spw=31) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)
    print(f"Found {len(weightfilelist)} HNCO 7m+12m+TP FITS weight files")
    log.debug(f"Full names of weights: {weightfilelist}")
    assert len(weightfilelist) == len(filelist)
    for xx, yy in zip(filelist, weightfilelist):
        #print(f'Beginning of filenames: {os.path.basename(xx.split(".")[0])}, {os.path.basename(yy.split("/")[-1].split(".")[0])}')
        print(f"Basenames: {os.path.basename(xx)}, {os.path.basename(yy)}")
        assert os.path.exists(xx)
        assert os.path.exists(yy)

    restfrq = 87.925238e9
    #cdelt_kms = 0.10409296373
    cdelt_kms = 0.20818593  # smooth by 2 chans
    for mode in ('mean', 'sum'):
        for ii in range(4):
            if False:
                # experiment for comparison: unweighted combination
                print("\nUnweighted mosaic")
                make_giant_mosaic_cube(filelist,
                                       reference_frequency=restfrq,
                                       cdelt_kms=cdelt_kms,
                                       cubename=f'HNCO_7m12mTP_alarzas_unweighted_{mode}',
                                       nchan=4,
                                       beam_threshold=3.3 * u.arcsec,
                                       target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                                       channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                                       weightfilelist=None,
                                       #use_reproject_cube=True,
                                       #parallel=4,
                                       channels=[ii],
                                       verbose=True,
                                       fail_if_cube_dropped=True,
                                       reproject_kwargs={'combine_function': mode})

            flatvalue_filelist = [fn.replace(".fits", ".flat-value.fits") for fn in filelist]
            for jj, (oldfn, fn) in enumerate(zip(filelist, flatvalue_filelist)):
                assert fn.endswith(".flat-value.fits")
                if not os.path.exists(fn):
                    shutil.copy(oldfn, fn)
                    with fits.open(fn, mode='update') as fh:
                        fh[0].data[:] = jj
                        fh.flush()

            if True:
                print("\nFlat (ones) mosaic")
                make_giant_mosaic_cube(flatvalue_filelist,
                                       reference_frequency=restfrq,
                                       cdelt_kms=cdelt_kms,
                                       cubename=f'HNCO_7m12mTP_alarzas_ones_{mode}',
                                       nchan=4,
                                       beam_threshold=3.3 * u.arcsec,
                                       target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                                       channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                                       weightfilelist=weightfilelist,
                                       #use_reproject_cube=True,
                                       #parallel=4,
                                       channels=[ii],
                                       fail_if_cube_dropped=True,
                                       verbose=True,
                                       reproject_kwargs={'combine_function': mode},
                                       )

            if True:
                print("\nOrder 1")
                make_giant_mosaic_cube(filelist,
                                    reference_frequency=restfrq,
                                    cdelt_kms=cdelt_kms,
                                    cubename=f'HNCO_7m12mTP_alarzas_order1_{mode}',
                                    nchan=4,
                                    beam_threshold=3.3 * u.arcsec,
                                    target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                                    channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                                    weightfilelist=weightfilelist,
                                    #use_reproject_cube=True,
                                    #parallel=4,
                                    channels=[ii],
                                    fail_if_cube_dropped=True,
                                    verbose=True,
                                    reproject_kwargs={'combine_function': mode},
                                    )


                # experiment: reverse order.  Result: null
                # make_giant_mosaic_cube(filelist[::-1],
                #                         reference_frequency=restfrq,
                #                         cdelt_kms=cdelt_kms,
                #                         cubename='HNCO_7m12mTP_alarzas_order2',
                #                         nchan=4,
                #                         beam_threshold=3.3 * u.arcsec,
                #                         target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                #                         channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                #                         weightfilelist=weightfilelist[::-1],
                #                         #use_reproject_cube=True,
                #                         #parallel=4,
                #                         channels=[ii],
                #                         fail_if_cube_dropped=True,
                #                         **kwargs,)

                # use_reproject_cube does not complete in finite time
                if False:
                    print("\nUse reproject cube")
                    make_giant_mosaic_cube(filelist,
                                        reference_frequency=restfrq,
                                        cdelt_kms=cdelt_kms,
                                        cubename=f'HNCO_7m12mTP_alarzas_usereprojectcube_{mode}',
                                        nchan=4,
                                        beam_threshold=3.3 * u.arcsec,
                                        target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                                        channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                                        weightfilelist=weightfilelist,
                                        use_reproject_cube=True,
                                        parallel=4,
                                        channels=[ii],
                                        fail_if_cube_dropped=True,
                                        verbose=True,
                                        reproject_kwargs={'combine_function': mode},
                                        )

                # no need to do two weights
                # make_giant_mosaic_cube(weightfilelist,
                #                         reference_frequency=restfrq,
                #                         cdelt_kms=cdelt_kms,
                #                         cubename='HNCO_7m12mTP_alarzas_order1_weights',
                #                         nchan=4,
                #                         beam_threshold=3.3 * u.arcsec,
                #                         target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                #                         channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                #                         weightfilelist=weightfilelist,
                #                         #use_reproject_cube=True,
                #                         #parallel=4,
                #                         channels=[ii],
                #                         fail_if_cube_dropped=True,
                #                         use_beams=False,
                #                         **kwargs,)

                print("\nMosaiced weighted-average weights")
                make_giant_mosaic_cube(weightfilelist[::-1],
                                    reference_frequency=restfrq,
                                    cdelt_kms=cdelt_kms,
                                    cubename=f'HNCO_7m12mTP_alarzas_order2_weights_{mode}',
                                    nchan=4,
                                    beam_threshold=3.3 * u.arcsec,
                                    target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_alarzas.hdr',
                                    channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                                    weightfilelist=weightfilelist[::-1],
                                    #use_reproject_cube=True,
                                    #parallel=4,
                                    channels=[ii],
                                    fail_if_cube_dropped=True,
                                    verbose=True,
                                    use_beams=False,
                                    reproject_kwargs={'combine_function': mode},
                                    )

if __name__ == "__main__":
    main()
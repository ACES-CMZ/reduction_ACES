import os
import shutil
from aces.imaging.make_mosaic import downsample_spectrally
from aces import conf

basepath = conf.basepath

dsfactor_dict = {
    "HNCO_7m12mTP": 7,
    "HCOP_noTP": 7,
    #"CH3CHO": 1,
    #"NSplus": 1,
    #"H40a": 1,
    "HC15N": 2,
    "SO21": 2,
    "H13CN": 2,
    "HN13C": 2,
    "H13COp": 2,
    #"CS21": 1,
    #"HC3N": 1,
    "HCOP": 7,
    "HCOP_mopra": 3,
    "SiO21": 2,
    #"SO32": 1,
}


def main():
    if os.getenv('MOLNAME'):
        molnames = [os.getenv('MOLNAME')]
    else:
        molnames = dsfactor_dict.keys()
    if os.getenv('USE_DASK'):
        use_dask = os.getenv('USE_DASK').lower() == 'true'
    else:
        use_dask = False

    if os.getenv('NUM_CORES'):
        numcores = int(os.getenv('NUM_CORES'))
    elif os.getenv('SLURM_NTASKS') and os.getenv('SLURM_JOB_CPUS_PER_NODE'):
        numcores = int(os.getenv('SLURM_NTASKS')) * int(os.getenv('SLURM_JOB_CPUS_PER_NODE'))
    else:
        numcores = 1

    print(f"Using Dask: {use_dask}")
    print(f"Number of cores: {numcores}")
    print(f"molnames: {molnames}")

    for molname in molnames:
        print(molname)
        target_path = '/red/adamginsburg/ACES/workdir/'
        cubename = f"{basepath}/mosaics/cubes/{molname}_CubeMosaic.fits"
        outname = f"{target_path}/mosaics/cubes/{molname}_CubeMosaic_spectrally.fits"
        downsample_spectrally(cubename,
                              outname,
                              factor=dsfactor_dict[molname],
                              use_dask=use_dask,
                              num_cores=numcores,
                              )
        if os.path.exists(outname):
            destination = f"{basepath}/mosaics/cubes/{os.path.basename(outname)}"
            if os.path.exists(destination):
                print(f"Overwriting destionation {destination}")
            print(f"Moving {outname} to {destination}")
            shutil.move(outname, destination)


if __name__ == "__main__":
    main()

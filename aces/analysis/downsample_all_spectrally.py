import os
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
"HCOP": 2,
"SiO21": 2,
#"SO32": 1,
}

    
def main():
    if os.getenv('MOLNAME'):
        molnames = [os.getenv('MOLNAME')]
    else:
        molnames = dsfactor_dict.keys()

    for molname in molnames:
        print(molname)
        downsample_spectrally(f"{basepath}/mosaics/cubes/{molname}_CubeMosaic.fits",
                              f"{basepath}/mosaics/cubes/{molname}_CubeMosaic_spectrally.fits",
                              factor=dsfactor_dict[molname],
                              )


if __name__ == "__main__":
    main()
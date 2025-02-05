from astropy.io import fits
import itertools
import numpy as np

import aces
basepath = aces.conf.basepath

molecules = ['CS21', 'HC3N', 'HCOP', 'HNCO_7m12mTP', 'SiO21']

if __name__ == '__main__':
    for suffix in ('max', 'mom0'):
        for mol1, mol2 in list(itertools.combinations(molecules, 2)):
            fh1 = fits.open(f'{basepath}/mosaics/cubes/moments/{mol1}_CubeMosaic_{suffix}.fits')
            fh2 = fits.open(f'{basepath}/mosaics/cubes/moments/{mol2}_CubeMosaic_{suffix}.fits')
            if np.all(fh2[0].data.shape == fh1[0].data.shape):
                rr = fh1[0].data / fh2[0].data
                fh1[0].data = rr
                fh1.writeto(f'{basepath}/mosaics/cubes/moments/ratios/{mol1}_over_{mol2}_{suffix}.fits', overwrite=True)
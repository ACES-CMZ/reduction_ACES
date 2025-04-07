import os
from glob import glob
from casatasks import importfits, feather, exportfits

os.system('rm -rf *.last')
os.system('rm -rf *.log')
os.system('rm -rf *.img')

inputfile = './../data/processed/cont_tp_cropped.fits'
print(inputfile)
importfits(inputfile, inputfile.replace('fits', 'img'), overwrite=True)

inputfile = './../data/processed/cont_12m_cropped.fits'
print(inputfile)
importfits(inputfile, inputfile.replace('fits', 'img'), overwrite=True)

feather(imagename='./../data/feathered/cont_12mtp.img',
        highres='./../data/processed/cont_12m_cropped.img',
        lowres='./../data/processed/cont_tp_cropped.img')

exportfits(imagename='./../data/feathered/cont_12mtp.img',
           fitsimage='./../data/feathered/cont_12mtp.fits',
           overwrite=True)
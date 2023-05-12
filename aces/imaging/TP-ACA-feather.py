# This script is to feather 7m and TP cubes.


# A template CMZ cube for mosaic
#templateIM='/almalustre/home/pyhsieh/ESO/spw24/template/CS2-1-temlate-whole-band.image'   # whole band
#templateIM='/almalustre/home/pyhsieh/ESO/spw24/template/CS2-1-temlate-narrow-band.image'  # only 0~1500 channel
templateIM='CS2-1-temlate-whole-band.image'

import glob
linename=glob.glob('*.image.pbcor')
linewt=glob.glob('*.weight')


fieldnum=38

# update rest frequency of 12m cubes (line and weight data) with CS(2-1) line
  
for i in range(len(linename)):
  imhead( imagename=linename[i],
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)

for i in range(len(linewt)):
  imhead( imagename=linewt[i],
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)



# reading spectral line and weight cubes of 7m and TP from dict.txt
# in dict-aca.txt, adding the pairs of 7m and TP mous id

daca={}
dtp={}
with open("dict-aca.txt") as f:
  for line in f:
    (key,valaca,valtp)=line.split()
    daca[int(key)]=valaca
    dtp[int(key)]=valtp


imageListACA=[]
imageListTP=[] 
imageListACAwt=[]


for i in range(fieldnum):
  imageListACAwt.append(glob.glob("*"+daca[i]+"*.weight"))

for i in range(fieldnum):
  imageListACA.append(glob.glob("*"+daca[i]+"*.image.pbcor"))

for i in range(fieldnum):
  imageListTP.append(glob.glob("*"+dtp[i]+"*.sd.fits.image"))                                                                                                                                     


# preparing a template CMZ cube
## using the mosaic template from globus under /mosaics
importfits(fitsimage = '7m_CS_2-1_m0_mosaic.fits',
           imagename = '7m_CS_2-1_m0_mosaic.image')

imhead(imagename='7m_CS_2-1_m0_mosaic.image',
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)
    
## regrid the mosaic template to Celestial coordiate (it appears that feather is not working well in Galactic coordinate)
imregrid(imagename = '7m_CS_2-1_m0_mosaic.image',
         template='J2000',
         output='7m_CS_2-1_m0_mosaic.radec.image')

## regrid (any) 7m line cube to the CMZ template (only X and Y)
imregrid(imagename = 'uid___A001_X1590_X30ac.s38_0.Sgr_A_star_sci.spw24.cube.I.iter1.image.pbcor',
         template='7m_CS_2-1_m0_mosaic.radec.image',
         output='CS2-1-temlate-whole-band.image', axes=[0, 1])

# Imtrans TP map, here needs to swap velocity and stokes axes to match the 7m data header

for i in range(fieldnum):
  imtrans(imagename=imageListTP[i][0],
    outfile=imageListTP[i][0]+'.imtrans',
    order="0132")

# regrid TP data to match 7m data
for i in range(fieldnum):
  imregrid(imagename=imageListTP[i][0]+'.imtrans',
    template=imageListACA[i][0],
    output=imageListTP[i][0]+'.imtrans.regrid',
    overwrite=False,
    axes=[-1])

# feather 7m and TP data
for i in range(fieldnum):
  feather(imagename=imageListACA[i][0]+'.feather',
    highres=imageListACA[i][0],
    lowres=imageListTP[i][0]+'.imtrans.regrid',
    sdfactor=1.0,effdishdiam=-1.0,lowpassfiltersd=False)

# smooth feather cubes to the same beam size (Here I check the maximum beam by visual inspection of 7m cube).
# To speed up this process, I only smooth data from channel 0 to 1550
for i in range(fieldnum):
  imsmooth(imagename=imageListACA[i][0]+'.feather',
    kernel="gauss",
    targetres=True,
    major="20arcsec",
    minor= "20arcsec",
    pa= "0deg",
    chans="0~1550", # updated this number depending on your need, or commend out this line if you want to use whole band
    outfile=imageListACA[i][0]+'.feather.imsmooth')

# create a narrower band to speed up mosaic, you don't need this step if you want to use whole band
for i in range(fieldnum):
  imsubimage(imagename=imageListACA[i][0]+'.feather.imsmooth',
    chans="369~935",
    outfile=imageListACA[i][0]+'.feather.imsmooth.imsubimage')

for i in range(fieldnum):
  imsubimage(imagename=imageListACAwt[i][0],
    chans="369~935",
    outfile=imageListACAwt[i][0]+'.imsubimage')
        
## regrid feather cubes to match the template for mosaic 

for i in range(fieldnum):
  imregrid(imagename=imageListACA[i][0]+'.feather.imsmooth.imsubimage',
    template=templateIM,
    output=imageListACA[i][0]+'.feather.imsmooth.imsubimage.regrid',
    overwrite=False,
    axes=[-1])

# deleting masks created during regriding, this step is necessary.
for i in range(fieldnum):
  ia.open(imageListACA[i][0]+'.feather.imsmooth.imsubimage.regrid')
  ia.maskhandler('delete',['mask0'])
  ia.close()


for i in range(fieldnum):
  imregrid(imagename=imageListACAwt[i][0]+'.imsubimage',
    template=templateIM,
    output=imageListACAwt[i][0]+'.regrid',
    overwrite=False,
    axes=[-1])

# deleting masks created during regriding,  this step is necessary.
for i in range(fieldnum):
  ia.open(imageListACAwt[i][0]+'.regrid')
  ia.maskhandler('delete',['mask0'])
  ia.close()

# mosacing maps with weighted averaged
# expor = (IM0*IM1+IM2*IM3+...)/(IM1+IM3+...)
# IM0, 2, 4,... are line cubes, IM1, 3, 5, ... are weight cubes
imageList=[] 
num = '(' 
den = '(' 
for i in range(fieldnum): 
  imageList.append(imageListACA[i][0]+'.feather.imsmooth.imsubimage.regrid')
  imageList.append(imageListACAwt[i][0]+'.regrid')
  tempNum = 'IM'+str(2*i)+'*IM'+str(2*i+1) 
  tempDen = 'IM'+str(2*i+1) 
  if i < len(imageListACA) -1: 
      num = num + tempNum + '+' 
      den = den + tempDen + '+' 
  else: 
      num = num + tempNum + ')' 
      den = den + tempDen + ')' 
expr = num + '/' + den  

immath(imagename=imageList,
  outfile='CS21-feather-Mosaic.image',
  expr=expr)

# convert coordinate from Celestial to Galactic coordinate.
imregrid(imagename='CS21-feather-Mosaic.image',
    template='Galactic',
    output='CS21-feather-Mosaic.gal.image')
    
# trim the blank area to reduce cube size
imsubimage(imagename='CS21-feather-Mosaic.gal.image',
  outfile='CS21-feather-Mosaic.gal.imsub.image',
  box='0,555,1595,1241')
      
# make mom8 map of individual field
for i in range(fieldnum):
  immoments(imagename=imageListACA[i][0]+'.feather',
    includepix=[20e-2,100],
    moments=[8],
    outfile=imageListACA[i][0]+'.feather.mom8',
    axis='spectral')

# convert coordinate from Celestial to Galactic coordinate of individual field
for i in range(fieldnum):
  imregrid(imagename=imageListACA[i][0]+'.feather.mom8',
    template='Galactic',
    output=imageListACA[i][0]+'.gal.feather.mom8')
    
for i in range(fieldnum):
  exportfits(imagename=imageListACA[i][0]+'.gal.feather.mom8',
    fitsimage=imageListACA[i][0]+'.gal.feather.mom8.fits')


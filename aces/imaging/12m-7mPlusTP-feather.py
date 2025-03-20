# This script is to feather 12m data with existing 7m+TP cubes (made in another script). The next step is to combine these two scripts. Basically they are doing the same thing.

# template image for CMZ mosaic (made the template cube below)
#templateIM='CS2-1-template-pm300-band.12mosaic.image'
#templateIM='CS2-1-temlate-pm180-band.image.zbin4'
templateIM='CS2-1-temlate-whole-band.image'


# reading spectral line and weight cubes from dict.txt
# in dict.txt, adding the pairs of 7m and 12m mous id
#--- example to dict.txt ---
#0 X1590_X30ac X1590_X30aa
#1 X15a0_Xae   X15a0_Xac
#2 X15a0_Xc0   X15a0_Xbe
#3 X15a0_Xc6   X15a0_Xc4
#4 ...         ...

fieldnum=25
daca={}
d12m={}
with open("dict.txt") as f:
  for line in f:
    (key,valaca,val12m)=line.split()
    daca[int(key)]=valaca
    d12m[int(key)]=val12m

imageListACA=[]
imageList12M=[] 
imageList12Mwt=[]


# 12m cubes of weight data
for i in range(fieldnum):
  imageList12Mwt.append(glob.glob("*"+d12m[i]+"*spw33*.weight.imsubimage"))

# feathered 7m+TP (here I trimmed the cube from +/- 480 km/s to speed up mosaic process, using task "imsubimage")
for i in range(fieldnum):
  imageListACA.append(glob.glob("*"+daca[i]+"*feather.pm400.imsubimage"))

# 12m cubes of target (here I trimmed the cube from +/- 180 km/s to speed up mosaic process, using task "imsubimage")
for i in range(fieldnum):
  imageList12M.append(glob.glob("*"+d12m[i]+"*spw33.cube.I.iter1*image.imsubimage"))

# update rest frequency of 12m cubes (line and weight data) with CS(2-1) line

for i in range(fieldnum):
  imhead(imagename=imageList12M[i][0],
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)

for i in range(fieldnum):
  imhead(imagename=imageList12Mwt[i][0],
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)

# preparing a template CMZ cube
## using the mosaic template from globus under /mosaics
importfits(fitsimage = '12m_CS_2-1_max_mosaic.fits',
           imagename = '12m_CS_2-1_max_mosaic.image')

imhead(imagename='12m_CS_2-1_max_mosaic.image',
    mode='put',
    hdkey='restfreq',hdvalue='97.980953GHz',verbose=False)
    
## regrid the mosaic template to Celestial coordiate (it appears that feather is not working well in Galactic coordinate)
imregrid(imagename = '12m_CS_2-1_max_mosaic.image',
         template='J2000',
         output='12m_CS_2-1_max_mosaic.radec.image')

## regrid one 12m line cube to the mosaic template, and using this regrided line cube as a template for stitching fields.
## Here the template has 2GHz bandwidth, but you could use "imsubimage" to create a narrow band image to speed up mosaic process

imregrid(imagename="uid___A001_X1590_X30aa.s38_0.Sgr_A_star_sci.spw33.cube.I.iter1.image.pbcor",
         template="12m_CS_2-1_max_mosaic.radec.image",
         output="CS2-1-temlate-whole-band.image", axes=[0, 1])


# regrid 7m+TP cubes with 12m cubes for each field  
for i in range(fieldnum):
  imregrid(imagename=imageListACA[i][0],
    template=imageList12M[i][0],
    output=imageListACA[i][0]+'.regrid',
    overwrite=False,
    axes=[-1])

# feather 7m+TP regrided cubes and 12m cubes
for i in range(fieldnum):
  feather(imagename=imageList12M[i][0]+'.ALL.feather',
    highres=imageList12M[i][0],
    lowres=imageListACA[i][0]+'.regrid',
    sdfactor=1.0,effdishdiam=-1.0,lowpassfiltersd=False)

# making moment maps
for i in range(fieldnum):
  immoments(imagename=imageList12M[i][0]+'.ALL.feather',
    includepix=[0,100],
    moments=[8],
    outfile=imageList12M[i][0]+'.feather.mom8',
    axis='spectral')

for i in range(fieldnum):
  exportfits(imagename=imageList12M[i][0]+'.feather.mom8',
    fitsimage=imageList12M[i][0]+'.feather.mom8.fits')

# smooth the feather cubes to the same resolution
for i in range(fieldnum):
  imsmooth(imagename=imageList12M[i][0]+'.ALL.feather',
    kernel="gauss",
    targetres=True,
    major="2.5arcsec",
    minor= "2.5arcsec",
    pa= "0deg",
    outfile=imageList12M[i][0]+'.ALL.feather.imsmooth')

# regrid the feathered cubes to a template map (CMZ)
for i in range(fieldnum):
  imregrid(imagename=imageList12M[i][0]+'.ALL.feather.imsmooth',
    template=templateIM,
    output=imageList12M[i][0]+'.ALL.feather.imsmooth.regrid',
    overwrite=False,
    axes=[-1])

# deleting masks created during regriding
for i in range(fieldnum):
  ia.open(imageList12M[i][0]+'.ALL.feather.imsmooth.regrid')
  ia.maskhandler('delete',['mask0'])
  ia.close()

# regrid the feathered weight cubes to a template map (CMZ)
for i in range(fieldnum):
  imregrid(imagename=imageList12Mwt[i][0],
    template=templateIM,
    output=imageList12Mwt[i][0]+'.regrid',
    overwrite=False,
    axes=[-1])

# deleting masks created during regriding
for i in range(fieldnum):
  ia.open(imageList12Mwt[i][0]+'.regrid')
  ia.maskhandler('delete',['mask0'])
  ia.close()

# mosacing maps with weighted averaged
# expor = (IM0*IM1+IM2*IM3+...)/(IM1+IM3+...)
# IM0, 2, 4,... are line cubes, IM1, 3, 5, ... are weight cubes

imageList=[] 
num = '(' 
den = '(' 
for i in range(fieldnum): 
  imageList.append(imageList12M[i][0]+'.ALL.feather.imsmooth.regrid')
  imageList.append(imageList12Mwt[i][0]+'.regrid')
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
  outfile='CS21-ALL-feather-Mosaic.SB25.image',
  expr=expr)

# convert the coordinate from celestial to galactic coordinate

imregrid(imagename='CS21-ALL-feather-Mosaic.SB25.image',
  template='Galactic',
  output='CS21-ALL-feather-Mosaic.SB25.gal.image')
  
# trim the blanked area to reduce size
imsubimage(imagename='CS21-ALL-feather-Mosaic.SB25.gal.image',
  box='0,4855,13684,10354',
  output='CS21-ALL-feather-Mosaic.SB25.gal.imsub.image')
  
# export fits of mosaic map

exportfits(imagename='CS21-ALL-feather-Mosaic.SB25.image',
    fitsimage='CS21-ALL-feather-Mosaic.SB25.image.fits')

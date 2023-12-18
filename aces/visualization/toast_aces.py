import os
from astropy import coordinates
from astropy import units as u
from radio_beam import Beam
from astropy.io import fits
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['figure.dpi'] = 300
from astropy import visualization
from astropy import wcs
import PIL
PIL.Image.MAX_IMAGE_PIXELS = 933120000


import matplotlib.colors as mcolors
import numpy as np
from astropy.visualization import simple_norm
import pyavm
import glob

from toasty import study, image as timage, pyramid, builder, merge
from wwt_data_formats import write_xml_doc, folder
from astropy.coordinates import SkyCoord

from astropy.wcs.utils import fit_wcs_from_points

def fits_to_avmpng(fitsfn, outpng):
    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ACES = fits.open(fitsfn)
    ACESheader = ACES[0].header
    ACESdata = ACES[0].data
    ACESwcs = wcs.WCS(ACES[0].header)
    ACESdata *= (1*u.Jy).to(u.K, Beam.from_fits_header(ACESheader).jtok_equiv(95*u.GHz)).value
    ACESheader['BUNIT'] = 'K'
    
    norm = simple_norm([0], min_cut=0.0001, max_cut=1.5, stretch='log')
    
    norm = simple_norm(mustangdata, min_cut=0.0001, max_cut=1.5, stretch='log')
    colordata_tens = mymap(norm(mustangdata))
    ct = (colordata_tens[::-1,:,:3] * 256).astype('uint8')
    ct[(colordata_tens[::-1,:,:3] * 256) > 255] = 255
    img_tens = PIL.Image.fromarray(ct)
    img_tens.save(outpng)
    avm = pyavm.AVM.from_wcs(mustangwcs)
    avm.embed(outpng, outpng)

def toast(imfn, targetdir='/orange/adamginsburg/web/public/ACES/toasts/'):

    img = PIL.Image.open(imfn)
    avm = pyavm.AVM.from_image(imfn)
    wcs = avm.to_wcs()

    tim = timage.Image.from_pil(img)
    data = np.array(img)

    #points = np.array([(0,0), (0, rslt.data.shape[0]), (rslt.data.shape[1], 0), (rslt.data.shape[1], rslt.data.shape[0])])
    #points = np.array([np.linspace(0, rslt.data.shape[1]), np.linspace(0, rslt.data.shape[0])]).T
    # overkill: 
    points = np.mgrid[0:data.shape[1]:1000, 0:data.shape[0]:1000]
    points = points.reshape(2, np.prod(points.shape[1:])).T
    wpoints = wcs.pixel_to_world(points[:,0], points[:,1])

    # +1 to convert from py->FITS
    wcsfk5 = fit_wcs_from_points((points[:,0]+1, points[:,1]+1), wpoints.fk5)
    #print(wcs, wcsfk5)

    # really dumb, definitely wrong hack to make the builder not complain about parity
    wcsfk5.wcs.cd[:,1] *= -1

    bui = builder.Builder(pyramid.PyramidIO(targetdir))
    stud = bui.prepare_study_tiling(tim)
    height, width, _ = np.array(img).shape
    if not os.path.exists(f'{targetdir}/0/0/0_0.png'):
        bui.execute_study_tiling(tim, stud)
        merge.cascade_images(
            bui.pio, start=7, merger=merge.averaging_merger, cli_progress=True
        )
    assert os.path.exists(f'{targetdir}/0/0/0_0.png'), "Failed"
    url = targetdir.replace("/orange/adamginsburg/web/public/", "https://data.rc.ufl.edu/pub/adamginsburg/")
    buisuf = bui.imgset.url
    bui.imgset.url = url + '/' + buisuf
    # not sure I have this right yet... there is a more sensible way to construct this URL
    bui.imgset.credits_url = bui.imgset.url.replace("pub", "secure").replace('/' + buisuf, '.fits')
    bui.apply_wcs_info(wcsfk5, width=width, height=height)
    bui.imgset.thumbnail_url = bui.imgset.url.format(0, 0, 0, 0)
    bui.imgset.name = os.path.basename(targetdir)

    fldr = bui.create_wtml_folder()
    
    # write both the 'rel' and 'full' URL versions
    bui.write_index_rel_wtml()
    with open(os.path.join(bui.pio._base_dir, "index.wtml"), 'w') as fh:
        write_xml_doc(fldr.to_xml(), dest_stream=fh)          
    print("Wrote ", os.path.join(bui.pio._base_dir, "index.wtml"))
    return os.path.join(bui.pio._base_dir, "index.wtml")

def make_all_indexes():
    import glob
    indexes = []
    for imfn in glob.glob("/orange/adamginsburg/ACES/mosaics/12m_flattened/*noaxes.png"):
        tdr = os.path.basename(imfn).replace("_noaxes.png", "")
        print(imfn, tdr)
        try:
            ind = toast(imfn, targetdir=f'/orange/adamginsburg/web/public/ACES/mosaics/12m_flattened/{tdr}')
        except Exception as ex:
            print(f'{imfn} failed')
            print(ex)
        indexes.append(ind)

    return indexes


def make_joint_index(indexes):
    fld = folder.Folder()
    fld.browseable = True
    fld.group = 'Explorer'
    fld.name = 'ACES'
    fld.searchable = True
    fld.thumbnail = 'https://data.rc.ufl.edu/pub/adamginsburg/ACES/mosaics/12m_flattened/12m_continuum_reimaged_spw33_35_mosaic/0/0/0_0.png'
    fld.type = 'SKY'
    fld.children = []

    acestens = folder.Folder.from_file("/orange/adamginsburg/web/public/ACES/MUSTANG_Feather/index_rel.wtml")
    acestens.children[0].name = 'ACES+TENS (toasty)'
    fld.children.extend(acestens.children)

    for ind in indexes:
        newfld = folder.Folder.from_file(ind)
        name = ind.split(os.sep)[-2].replace("_mosaic", "")
        assert '/' not in name
        newfld.children[0].name = name
        for child in newfld.children:
            child.thumbnail = child.foreground_image_set.url.format("", "0", "0", "0")
        fld.children.extend(newfld.children)
        with open('/orange/adamginsburg/web/public/ACES/mosaics/mosaics.wtml', 'w') as fh:
            write_xml_doc(fld.to_xml(), dest_stream=fh)


if __name__ == "__main__":
    indexes = make_all_indexes()

    make_joint_index(indexes)

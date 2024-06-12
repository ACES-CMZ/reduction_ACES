import os
from astropy import coordinates
from astropy import units as u
from radio_beam import Beam
from astropy.io import fits
import pylab as pl
from astropy import visualization
from astropy import wcs
import PIL

import matplotlib.colors as mcolors
import numpy as np
from astropy.visualization import simple_norm
import pyavm
import glob

from toasty import study, image as timage, pyramid, builder, merge
from wwt_data_formats import write_xml_doc, folder
from astropy.coordinates import SkyCoord

from astropy.wcs.utils import fit_wcs_from_points

pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['figure.dpi'] = 300
PIL.Image.MAX_IMAGE_PIXELS = 933120000


def fits_to_avmpng(fitsfn, outpng):
    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    ACES = fits.open(fitsfn)
    ACESheader = ACES[0].header
    ACESdata = ACES[0].data
    ACESwcs = wcs.WCS(ACES[0].header)
    ACESdata *= (1 * u.Jy).to(u.K, Beam.from_fits_header(ACESheader).jtok_equiv(95 * u.GHz)).value
    ACESheader['BUNIT'] = 'K'

    norm = simple_norm([0], min_cut=0.0001, max_cut=1.5, stretch='log')

    norm = simple_norm(ACESdata, min_cut=0.0001, max_cut=1.5, stretch='log')
    colordata_tens = mymap(norm(ACESdata))
    ct = (colordata_tens[::-1, :, :3] * 256).astype('uint8')
    ct[(colordata_tens[::-1, :, :3] * 256) > 255] = 255
    img_tens = PIL.Image.fromarray(ct)
    img_tens.save(outpng)
    avm = pyavm.AVM.from_wcs(ACESwcs)
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
    wpoints = wcs.pixel_to_world(points[:, 0], points[:, 1])

    wcsfk5 = fit_wcs_from_points((points[:, 0], points[:, 1]), wpoints.fk5)

    # sanity check
    gc = SkyCoord(0 * u.deg, 0 * u.deg, frame='galactic')
    print(f"Sanity check: 0,0 gal -> pix in orig: {wcs.world_to_pixel(gc)}, in fk5: {wcsfk5.world_to_pixel(gc)}"
          f" should be the same! diff is {np.array(wcs.world_to_pixel(gc)) - np.array(wcsfk5.world_to_pixel(gc))}")

    height, width, _ = np.array(img).shape

    # redo: using world_to_pixel makes the cd flip not work (probably a wcs.set call populates pc instead of cd)
    wcsfk5 = timage._flip_wcs_parity(fit_wcs_from_points((points[:, 0], points[:, 1]), wpoints.fk5),
                                     height - 1)
    # really dumb, definitely wrong hack to make the builder not complain about parity
    # wcsfk5.wcs.cd[:,1] *= -1
    print(f"Sanity check: 0,0 gal from the flipped fk5: {wcsfk5.world_to_pixel(gc)[1]} should match {wcs.world_to_pixel(gc)[1]} if we do naxis-y: {height - np.array(wcsfk5.world_to_pixel(gc)[1])}")

    bui = builder.Builder(pyramid.PyramidIO(targetdir))
    stud = bui.prepare_study_tiling(tim)
#    if 'HNCO' in targetdir or not os.path.exists(f'{targetdir}/0/0/0_0.png'):
    #if not os.path.exists(f'{targetdir}/0/0/0_0.png'):
    if True:  # always redo
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

    ctr = SkyCoord(0.1189 * u.deg, -0.05505 * u.deg, frame='galactic').fk5
    bui.place.ra_hr = ctr.ra.hourangle
    bui.place.dec_deg = ctr.dec.deg
    bui.place.zoom_level = 4

    fldr = bui.create_wtml_folder()

    # write both the 'rel' and 'full' URL versions
    bui.write_index_rel_wtml()
    with open(os.path.join(bui.pio._base_dir, "index.wtml"), 'w') as fh:
        write_xml_doc(fldr.to_xml(), dest_stream=fh)
    print("Wrote ", os.path.join(bui.pio._base_dir, "index.wtml"))
    return os.path.join(bui.pio._base_dir, "index.wtml")


def make_all_indexes():
    indexes = []
    for imfn in (glob.glob("/orange/adamginsburg/ACES/mosaics/continuum/*continuum*noaxes.png") +
                 glob.glob("/orange/adamginsburg/ACES/mosaics/cubes/moments/*png")):
        if 'residual' in imfn:
            continue
        tdr = os.path.basename(imfn).replace("_noaxes.png", "").replace(".png", "")
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
    acestens.children[0].thumbnail = 'https://data.rc.ufl.edu/pub/adamginsburg/ACES/MUSTANG_Feather/0/0/0_0.png'

    ctr = SkyCoord(0.1189 * u.deg, -0.05505 * u.deg, frame='galactic').fk5
    acestens.children[0].ra_hr = ctr.ra.hourangle
    acestens.children[0].dec_deg = ctr.dec.deg
    acestens.children[0].zoom_level = 4

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


def main():
    indexes = make_all_indexes()

    make_joint_index(indexes)


if __name__ == "__main__":
    main()

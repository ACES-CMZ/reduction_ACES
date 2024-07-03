import numpy as np
from astropy.table import Table
import scipy.stats
import scipy
from radio_beam import Beam
from casatools import msmetadata, ms as mstool
from casatasks import imhead
import glob
import pylab as pl
from astropy import units as u
from astropy.io import fits

import sys,os

def savefig(path, bbox_inches='tight', **kwargs):
    pl.savefig(path, bbox_inches=bbox_inches, **kwargs)
    pl.savefig(path.replace(".pdf", ".png"), bbox_inches=bbox_inches, **kwargs)

def make_figure(data, wavelength, beam, bins=50):
    if len(data) == 0:
        print("FAILURE: data were empty")
        return
    uvcts = np.concatenate([data[spw]['uvdist'][~data[spw]['flag'].any(axis=(0,1))] for spw in data]).ravel()
    uvwts = np.concatenate([data[spw]['weight'].mean(axis=0)[~data[spw]['flag'].any(axis=(0,1))] for spw in data]).ravel()

    beam_to_bl = (wavelength / beam).to(u.m, u.dimensionless_angles())

    pl.figure(figsize=(8,4))
    ax1 = pl.subplot(1,2,1)
    _=pl.hist(uvcts, bins=bins)
    _=pl.xlabel('Baseline Length (m)')
    _=pl.ylabel("Number of Visibilities")
    yl = pl.ylim()

    #try:
    #    pl.fill_betweenx(yl, beam_to_bl[0].value, beam_to_bl[1].value, zorder=-5, color='orange', alpha=0.5)
    #except TypeError:
    #    pl.axvline(beam_to_bl.value, color='orange', zorder=-5, alpha=0.5)
    pl.fill_betweenx(yl, np.percentile(uvcts, 25), np.percentile(uvcts, 75), zorder=-5, color='red', alpha=0.25)

    pl.ylim(yl)
    ax1t = ax1.secondary_xaxis('top', functions=(lambda x: x/1e3/wavelength.to(u.m).value, lambda x:x/1e3/wavelength.to(u.m).value))
    ax1t.set_xlabel("Baseline Length (k$\lambda$)")
    #ax1t.set_ticks(np.linspace(1000,100000,10))
    ax2 = pl.subplot(1,2,2)
    _=pl.hist(uvcts,
            weights=uvwts,
            bins=bins, density=True)
    _=pl.xlabel('Baseline Length (m)')
    _=pl.ylabel("Fractional Weight")
    def forward(x):
        return (wavelength.to(u.m)/(x*u.arcsec)).to(u.m, u.dimensionless_angles()).value
    def inverse(x):
        return (wavelength.to(u.m)/(x*u.m)).to(u.arcsec, u.dimensionless_angles()).value
    ax2t = ax2.secondary_xaxis('top', functions=(forward, inverse))
    ax2t.set_xlabel("Angular size $\lambda/D$ (arcsec)")
    if ax2.get_xlim()[1] > 1000:
        ax2t.set_ticks([10,1,0.5,0.4,0.3,0.2,0.1])
    elif ax2.get_xlim()[1] > 600:
        ax2t.set_ticks([10,2,1,0.6,0.5,0.4,0.3])
    else:
        ax2t.set_ticks([10,2,1,0.8,0.7,0.6,0.2])
    yl = pl.ylim()
    try:
        ax2.fill_betweenx(yl, beam_to_bl[0].value, beam_to_bl[1].value, zorder=-5, color='orange', alpha=0.5)
    except TypeError:
        ax2.axvline(beam_to_bl.value, color='orange', zorder=-5, alpha=0.5)
    ax2.fill_betweenx(yl, np.percentile(uvcts, 25), np.percentile(uvcts, 75), zorder=-5, color='red', alpha=0.25)
    ax2.set_ylim(yl)
    #pl.subplots_adjust(wspace=0.3)
    pl.tight_layout()

    print(f"25th pctile={forward(np.percentile(uvcts, 25))}, 75th pctile={forward(np.percentile(uvcts, 75))}")
    try:
        return (forward(np.percentile(uvcts,
                                      [1,5,10,25,50,75,90,95,99])),
                scipy.stats.percentileofscore(uvcts, beam_to_bl[0].value),
                scipy.stats.percentileofscore(uvcts, beam_to_bl[1].value))
    except TypeError:
        return (forward(np.percentile(uvcts,
                                      [1,5,10,25,50,75,90,95,99])),
                scipy.stats.percentileofscore(uvcts, beam_to_bl.value),
                scipy.stats.percentileofscore(uvcts, beam_to_bl.value))

def main():
    basepath = '/orange/adamginsburg/ACES/'
    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')

    # /orange/adamginsburg/ACES//data//2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X1*/calibrated/working/*ms
    mslist = {row['Obs ID']:
              glob.glob(f'{basepath}/data//2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{row["12m MOUS ID"]}/calibrated/working/*target.ms')
              for row in tbl}

    msmd = msmetadata()
    ms = mstool()

    uvdata = []

    for row in tbl:
        region = row['Obs ID']
        datapath = f'{basepath}/data//2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{row["12m MOUS ID"]}/calibrated/working'
        data = {}
        for msname in mslist[region]:

            msmd.open(msname)
            spws = msmd.spwsforfield('Sgr_A_star')
            freqs = np.concatenate([msmd.chanfreqs(spw) for spw in spws])
            freqweights = np.concatenate([msmd.chanfreqs(spw) for spw in spws])
            msmd.close()
            print(msname)

            avfreq = np.average(freqs, weights=freqweights)
            wavelength = (avfreq*u.Hz).to(u.m, u.spectral())

            for spw in spws:
                ms.open(msname)
                ms.selectinit(spw)
                newdata = ms.getdata(items=['weight', 'uvdist', 'flag'])
                if len(newdata) > 0:
                    if spw in data:
                        for key in data[spw]:
                            data[spw][key] = np.concatenate([data[spw][key], newdata[key]], axis=-1)
                            print(f'{key}: {data[spw][key].shape}')
                    else:
                        data[spw] = newdata
                ms.close()

        if len(data) == 0:
            print(f"FAILURE FOR REGION {region}: len(data)=0")
            continue

        try:
            fname = glob.glob(f'{datapath}/*.spw25_27_29_31_33_35.cont.I.iter1.image.tt0.pbcor.fits')[0]
            beam = Beam.from_fits_header(fits.getheader(fname))
        except IndexError:
            fname = glob.glob(f'{datapath}/*.spw25_27_29_31_33_35.cont.I.iter1.image.tt0.pbcor')[0]
            bmaj = imhead(fname, mode='get', hdkey='beammajor')['value']
            bmin = imhead(fname, mode='get', hdkey='beamminor')['value']
            bpa = imhead(fname, mode='get', hdkey='beampa')['value']
            beam = Beam(major=bmaj*u.arcsec, minor=bmin*u.arcsec, pa=bpa*u.deg)

        print('beam: ',beam)
        with np.errstate(divide='ignore'):
            pctiles,majpct,minpct = make_figure(data, wavelength, beam)
        pl.suptitle(f"{region}")
        savefig(f'{basepath}/diagnostic_plots/uvhistograms/{region}_uvhistogram.pdf', bbox_inches='tight')

        uvdata.append({
                       'region': region,
                       '1%': pctiles[0],
                       '5%': pctiles[1],
                       '10%': pctiles[2],
                       '25%': pctiles[3],
                       '50%': pctiles[4],
                       '75%': pctiles[5],
                       '90%': pctiles[6],
                       '95%': pctiles[7],
                       '99%': pctiles[8],
                       'beam_major': beam.major,
                       'beam_minor': beam.minor,
                       'beam_major_pctile': majpct,
                       'beam_minor_pctile': minpct,
                       'wavelength': wavelength.to(u.um).value,
                       })
        print(uvdata[-1])
        # debug
        #break
    uvtbl = Table(uvdata, units={'beam_major':u.arcsec, 'beam_minor':u.arcsec, 'wavelength':u.um, '1%':u.arcsec, '5%':u.arcsec, '10%':u.arcsec, '25%':u.arcsec, '50%':u.arcsec, '75%':u.arcsec, '90%':u.arcsec, '95%':u.arcsec, '99%':u.arcsec})
    uvtbl.write(f'{basepath}/reduction_ACES/aces/data/tables/uvspacings.ecsv', overwrite=True)

    fontsize=16
    pl.rcParams['font.size'] = fontsize

    rows = {key: uvtbl[(uvtbl['region']==key)] for key in sorted(uvtbl['region'])}
    stats = [{
        "label": key,
        "med": row['50%'][0],
        "q1": row['25%'][0],
        "q3": row['75%'][0],
        "whislo": row['5%'][0],
        "whishi": row['95%'][0],
        "fliers": [],
        } for key,row in rows.items() if len(row)>0][::-1]

    fig, axes = pl.subplots(nrows=1, ncols=1, figsize=(12, 12), sharey=True)
    axes.bxp(stats, vert=False)
    #axes.set_title(f'{band} UV distribution overview', fontsize=fontsize)
    axes.set_xlabel("Angular Scale (\")", fontsize=fontsize)
    axes.set_xlim(0.1,20)
    rad_to_as = u.radian.to(u.arcsec)
    def fcn(x):
        return rad_to_as / x / 1000
    ax1t = axes.secondary_xaxis('top', functions=(fcn, fcn))
    ax1t.set_xlabel("Baseline Length (k$\lambda$)")
    ax1t.set_ticks([15, 20, 30, 50, 100, 400])
    savefig(f'{basepath}/diagnostic_plots/uvhistograms/summary_uvdistribution.pdf', bbox_inches='tight')

    distance = 8.5

    stats = [{
        "label": key,
        "med": row['50%'][0] * distance*1000,
        "q1": row['25%'][0] * distance*1000,
        "q3": row['75%'][0] * distance*1000,
        "whislo": row['5%'][0] * distance*1000,
        "whishi": row['95%'][0] * distance*1000,
        "fliers": [],
        } for key,row in rows.items() if len(row)>0][::-1]

    fig, axes = pl.subplots(nrows=1, ncols=1, figsize=(12, 12), sharey=True)
    axes.bxp(stats, vert=False)
    #axes.set_title(f'{band} UV distribution overview', fontsize=fontsize)
    axes.set_xlabel("Physical Scale (au)", fontsize=fontsize)
    axes.set_xlim(0, 200000)
    axes.set_xticks([0, 2000, 10000, 20000, 30000, 40000, 50000])
    rad_to_as = u.radian.to(u.arcsec)
    def fcn(x):
        return x / rad_to_as
    ax1t = axes.secondary_xaxis('top', functions=(fcn, fcn))
    ax1t.set_xlabel("Physical Scale (pc)")
    ax1t.set_ticks([0.005, 0.1, 0.2])
    savefig(f'{basepath}/diagnostic_plots/uvhistograms/summary_uvdistribution_physicalscale.pdf', bbox_inches='tight')

if  __name__ == "__main__":
    main()
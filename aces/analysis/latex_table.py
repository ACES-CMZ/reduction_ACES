import numpy as np
from astropy.table import Table, Column
from astropy import table
import requests
import keyring
import datetime
from astropy import units as u
import sigfig
from radio_beam import Beam

from aces import conf
basepath = conf.basepath

from aces.analysis.latex_info import (latexdict, format_float, round_to_n, rounded,
                        rounded_arr, strip_trailing_zeros, exp_to_tex)

latexdict = latexdict.copy()

def make_latex_table(savename='continuum_data_summary'):

    tbl = Table.read(f'{basepath}/tables/metadata_image.tt0.ecsv')
    tbl['mad'] = tbl['mad'] * 1000
    tbl['peak'] = (tbl['peak/mad'] * tbl['mad'])
    tbl['casaversion'] = [[x for x in y.split() if 'CASA' not in x][0] for y in tbl['casaversion']]
    tbl['spws'] = ['all' if x == '25,27,29,31,33,35' else x for x in tbl['spws']]

    good = ((np.char.find(tbl['filename'], 'preQA3')==-1) |
            (np.char.find(tbl['filename'], 'X15a0_X178')==-1) |
            (np.char.find(tbl['filename'], 'obsolete')==-1) #| False # TODO: check for more
            ) & (
            (np.char.find(tbl['spws'], 'v1') == -1) &
            ~(tbl['pbcor']) &
            (tbl['array'] == '12M')
            )

    tbl = tbl[good]
    assert len(tbl) > 0

    # filter out v1 first
    nueff = {'all': 97.271*u.GHz, '25,27': 86.539*u.GHz, '33,35': 99.543*u.GHz}
    jtok = u.Quantity([Beam(x['bmaj']*u.arcsec, x['bmin']*u.arcsec, x['bpa']).jtok(nueff[x['spws']]) for x in tbl])
    tbl['peak_K'] = (tbl['peak'] * jtok / 1000)
    tbl['mad_K'] = (tbl['mad'] * jtok)

    cols_to_keep = {'region':'Region',
                    'bmaj':r'$\theta_{maj}$',
                    'bmin':r'$\theta_{min}$',
                    'bpa':'BPA',
                    'spws':'SPWs',
                    #'Req_Res': r"$\theta_{req}$",
                    #'BeamVsReq': r"$\Omega_{syn}^{1/2}/\Omega_{req}^{1/2}$",
                    #'peak/mad': "DR",
                    'casaversion': 'CASA Version',
                    'peak':'$S_{peak}$',
                    'mad':r'$\sigma_{MAD}$',
                    'peak_K':'$T_{B,peak}$',
                    'mad_K':r'$\sigma_{MAD,mK}$',
                    #'Req_Sens': r"$\sigma_{req}$",
                    #'SensVsReq': r"$\sigma_{MAD}/\sigma_{req}$",
                    #'dr_pre': "DR$_{pre}$",
                    #'dr_post': "DR$_{post}$",
                    #'dr_improvement': "DR$_{post}$/DR$_{pre}$"
                    }

    units = {'$S_{peak}$':(u.mJy/u.beam).to_string(u.format.LatexInline),
             '$T_{B,peak}$':(u.K).to_string(u.format.LatexInline),
             r'$\sigma_{MAD}$':(u.mJy/u.beam).to_string(u.format.LatexInline),
             r'$\sigma_{MAD,mK}$':(u.mK).to_string(u.format.LatexInline),
             r'$\sigma_{req}$':(u.mJy/u.beam).to_string(u.format.LatexInline),
             r'$\theta_{req}$':u.arcsec.to_string(u.format.LatexInline),
             r'$\theta_{maj}$':u.arcsec.to_string(u.format.LatexInline),
             r'$\theta_{min}$':u.arcsec.to_string(u.format.LatexInline),
             r'PA':u.deg.to_string(u.format.LatexInline),
             r'BPA':u.deg.to_string(u.format.LatexInline),
            }
    latexdict['units'] = units

    ftbl = tbl[list(cols_to_keep.keys())]


    for old, new in cols_to_keep.items():
        if old in tbl.colnames:
            #tbl[old].meta['description'] = description[old]
            ftbl.rename_column(old, new)
            if new in units:
                ftbl[new].unit = units[new]

    float_cols =  ['$\\theta_{maj}$',
     '$\\theta_{min}$',
     'BPA',
     '$S_{peak}$',
     '$T_{B,peak}$',
     '$\\sigma_{MAD}$',
     '$\\sigma_{MAD,mK}$',
     '$\\theta_{req}$',
     '\\sigma_{req}$',
     '$\\sigma_{MAD}/\\sigma_{req}$',
    # '$\\theta_{req}/\\theta_{maj}$',
     #r"$\Omega_{syn}^{1/2}/\Omega_{req}^{1/2}$",
     #'DR$_{pre}$',
     #'DR$_{post}$',
     #'DR$_{post}$/DR$_{pre}$'
     ]

    # ALREADY IN mJy
    # convert to mJy
    #ftbl['$\sigma_{MAD}$'] *= 1000


    formats = {key: lambda x: strip_trailing_zeros('{0:0.3f}'.format(round_to_n(x,2), nsf=2))
               for key in float_cols}
    formats = {key: lambda x: str(sigfig.round(str(x), sigfigs=2))
               for key in float_cols}

    ftbl.write('continuum_data_summary_table.ecsv', format='ascii.ecsv', overwrite=True)



    # caption needs to be *before* preamble.
    #latexdict['caption'] = 'Continuum Source IDs and photometry'
    latexdict['header_start'] = '\\label{tab:continuum_data_summary}'#\n\\footnotesize'
    latexdict['preamble'] = '\\caption{Continuum Image Summary}\n\\resizebox{\\textwidth}{!}{'
    latexdict['col_align'] = 'l'*len(ftbl.columns)
    latexdict['tabletype'] = 'table*'
    latexdict['tablefoot'] = ("}\\par\n"
                             )

    ftbl.sort(['Region'])

    ftbl.write(f"{basepath}/papers/continuum_data/{savename}.tex", formats=formats,
               overwrite=True, latexdict=latexdict)

    return ftbl

if __name__ == "__main__":

    result = make_latex_table()
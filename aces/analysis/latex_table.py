import numpy as np
from astropy.table import Table, Column
from astropy import table
import requests
import keyring
import datetime
from astropy import units as u
import sigfig
from radio_beam import Beam

from aces.analysis.latex_info import (latexdict, format_float, round_to_n, rounded,
                                      rounded_arr, strip_trailing_zeros, exp_to_tex)

from aces import conf
basepath = conf.basepath

latexdict = latexdict.copy()


def make_latex_table(savename='continuum_data_summary'):

    tbl = Table.read(f'{basepath}/tables/metadata_image.tt0.ecsv')
    tbl['mad'] = tbl['mad'] * 1000
    tbl['peak'] = (tbl['peak/mad'] * tbl['mad'])
    tbl['casaversion'] = [[x for x in y.split() if 'CASA' not in x][0] for y in tbl['casaversion']]
    tbl['spws'] = ['all' if x == '25,27,29,31,33,35' else x for x in tbl['spws']]

    good = (((np.char.find(tbl['filename'], 'preQA3') == -1) |
            (np.char.find(tbl['filename'], 'X15a0_X178') == -1) |
             (np.char.find(tbl['filename'], 'obsolete') == -1) #| False # TODO: check for more
             ) & (
             (np.char.find(tbl['spws'], 'v1') == -1) &
             ~(tbl['pbcor']) &
             (tbl['array'] == '12M'))
    )

    tbl = tbl[good]
    assert len(tbl) > 0

    # filter out v1 first
    nueff = {'all': 97.271 * u.GHz, '25,27': 86.539 * u.GHz, '33,35': 99.543 * u.GHz}
    jtok = u.Quantity([Beam(x['bmaj'] * u.arcsec, x['bmin'] * u.arcsec, x['bpa']).jtok(nueff[x['spws']]) for x in tbl])
    tbl['peak_K'] = (tbl['peak'] * jtok / 1000)
    tbl['mad_K'] = (tbl['mad'] * jtok)

    cols_to_keep = {'region': 'Region',
                    'bmaj': r'$\theta_{maj}$',
                    'bmin': r'$\theta_{min}$',
                    'bpa': 'BPA',
                    'robust': 'Robust',
                    'spws': 'SPWs',
                    #'Req_Res': r"$\theta_{req}$",
                    #'BeamVsReq': r"$\Omega_{syn}^{1/2}/\Omega_{req}^{1/2}$",
                    #'peak/mad': "DR",
                    'casaversion': 'CASA Version',
                    'peak': '$S_{peak}$',
                    'mad': r'$\sigma_{MAD}$',
                    'peak_K': '$T_{B,peak}$',
                    'mad_K': r'$\sigma_{MAD,mK}$',
                    #'Req_Sens': r"$\sigma_{req}$",
                    #'SensVsReq': r"$\sigma_{MAD}/\sigma_{req}$",
                    #'dr_pre': "DR$_{pre}$",
                    #'dr_post': "DR$_{post}$",
                    #'dr_improvement': "DR$_{post}$/DR$_{pre}$"
                    }

    units = {'$S_{peak}$': (u.mJy / u.beam).to_string(u.format.LatexInline),
             '$T_{B,peak}$': (u.K).to_string(u.format.LatexInline),
             r'$\sigma_{MAD}$': (u.mJy / u.beam).to_string(u.format.LatexInline),
             r'$\sigma_{MAD,mK}$': (u.mK).to_string(u.format.LatexInline),
             #r'$\sigma_{req}$': (u.mJy / u.beam).to_string(u.format.LatexInline),
             #r'$\theta_{req}$': u.arcsec.to_string(u.format.LatexInline),
             r'$\theta_{maj}$': u.arcsec.to_string(u.format.LatexInline),
             r'$\theta_{min}$': u.arcsec.to_string(u.format.LatexInline),
             r'PA': u.deg.to_string(u.format.LatexInline),
             r'BPA': u.deg.to_string(u.format.LatexInline),
            }
    latexdict['units'] = units

    ftbl = tbl[list(cols_to_keep.keys())]

    for old, new in cols_to_keep.items():
        if old in tbl.colnames:
            #tbl[old].meta['description'] = description[old]
            ftbl.rename_column(old, new)
            if new in units:
                ftbl[new].unit = units[new]

    float_cols = ['$\\theta_{maj}$',
                  '$\\theta_{min}$',
                  'BPA',
                  '$S_{peak}$',
                  '$T_{B,peak}$',
                  '$\\sigma_{MAD}$',
                  '$\\sigma_{MAD,mK}$',
                  #'$\\theta_{req}$',
                  #'\\sigma_{req}$',
                  #'$\\sigma_{MAD}/\\sigma_{req}$',
                 ]

    formats = {key: lambda x: strip_trailing_zeros('{0:0.3f}'.format(round_to_n(x, 2)), nsf=2)
               for key in float_cols}
    formats = {key: lambda x: str(sigfig.round(str(x), sigfigs=2))
               for key in float_cols}

    ftbl.write(f'{basepath}/tables/continuum_data_summary_table.ecsv', format='ascii.ecsv', overwrite=True)

    # caption needs to be *before* preamble.
    #latexdict['caption'] = 'Continuum Source IDs and photometry'
    latexdict['header_start'] = '\\label{tab:continuum_data_summary}'#\n\\footnotesize'
    latexdict['preamble'] = '\\caption{Continuum Image Summary}\n\\resizebox{\\textwidth}{!}{'
    latexdict['col_align'] = 'l' * len(ftbl.columns)
    latexdict['tabletype'] = 'table*'
    latexdict['tablefoot'] = ("}\\par\n"
                             )

    ftbl.sort(['Region'])

    ftbl.write(f"{basepath}/papers/continuum_data/{savename}.tex", formats=formats,
               overwrite=True, latexdict=latexdict)

    return ftbl


def make_observation_table():
    from astroquery.alma import Alma
    from astropy.time import Time
    import numpy as np
    from astropy.table import Table

    metadata = Alma.query(payload={'project_code': '2021.1.00172.L'})
    metadata['t_min'] = Time(metadata['t_min'], format='mjd')
    metadata['Observation Start Time'] = Time(metadata['t_min'], format='mjd')
    metadata['Observation Start Time'].format = 'fits'
    metadata['Observation End Time'] = Time(metadata['t_max'], format='mjd')
    metadata['Observation End Time'].format = 'fits'

    sub_meta = metadata['schedblock_name', 'Observation Start Time',
                        'Observation End Time', 'pwv', 't_exptime',
                        's_resolution', 'spatial_scale_max']
    usub_meta = Table(np.unique(sub_meta))

    tm = np.char.find(usub_meta['schedblock_name'], '_TM') > 0
    sm = np.char.find(usub_meta['schedblock_name'], '_7M') > 0
    tp = np.char.find(usub_meta['schedblock_name'], '_TP') > 0

    configuration_schedules = requests.get('https://almascience.eso.org/observing/observing-configuration-schedule/prior-cycle-observing-and-configuration-schedule')

    from astropy.time import Time
    from astropy.io import ascii

    def try_read_table(index):
        try:
            return ascii.read(configuration_schedules.text, format='html', htmldict={'table_id': index})
        except Exception as ex:  # noqa
            #print(ex)
            return

    config_sched_tbls = {ii:
                         try_read_table(ii)
                         for ii in range(1, 10)
                         if try_read_table(ii) is not None
                        }
    assert len(config_sched_tbls) > 0

    for ii, config_sched_tbl in config_sched_tbls.items():
        config_sched_tbl['Start date'] = Time(config_sched_tbl['Start date'])
        config_sched_tbl['End date'] = Time(config_sched_tbl['End date'])

    usub_meta['cycle_id'] = np.zeros(len(usub_meta), dtype='int')
    usub_meta['config_id'] = ['x' * 11 for x in range(len(usub_meta))]
    usub_meta['block'] = ['xx' for x in range(len(usub_meta))]
    for row in usub_meta:
        for cycle_id, config_sched_tbl in config_sched_tbls.items():
            match = ((config_sched_tbl['Start date'] < row['Observation Start Time']) &
                     (config_sched_tbl['End date'] >= row['Observation Start Time']))
            if match.sum() > 0:
                row['cycle_id'] = cycle_id
                row['config_id'] = config_sched_tbl[match]['Approx Config.'][0]
                row['block'] = config_sched_tbl[match]['Block'][0]
                #print(f"Matched {row['schedblock_name']} to cycle {cycle_id} config {row['config_id']} block {row['block']}")

    usub_meta['Field'] = [x.split("_")[3] for x in usub_meta['schedblock_name']]
    usub_meta.rename_column('pwv', 'PWV')
    usub_meta.rename_column('t_exptime', 'Exposure Time')
    usub_meta.rename_column('s_resolution', 'Resolution')
    usub_meta.rename_column('spatial_scale_max', 'LAS')
    usub_meta.rename_column('config_id', 'Configuration')

    colnames = ['Field', 'Observation Start Time', 'Configuration', 'PWV', 'Exposure Time', 'Resolution', 'LAS']
    usub_meta['LAS'].unit = u.arcsec
    usub_meta['Resolution'].unit = u.arcsec
    usub_meta['Exposure Time'].unit = u.s
    usub_meta['PWV'].unit = u.mm
    usub_meta[colnames][tm].write(f'{basepath}/tables/observation_metadata_12m.ecsv', format='ascii.ecsv', overwrite=True)
    usub_meta[colnames][sm].write(f'{basepath}/tables/observation_metadata_7m.ecsv', format='ascii.ecsv', overwrite=True)
    usub_meta[colnames][tp].write(f'{basepath}/tables/observation_metadata_TP.ecsv', format='ascii.ecsv', overwrite=True)

    #latexdict['header_start'] = '\\label{tab:observation_metadata_12m}'#\n\\footnotesize'
    #latexdict['preamble'] = '\\caption{12m Observation Metadata}\n\\resizebox{\\textwidth}{!}{'
    latexdict['col_align'] = 'l' * len(usub_meta.columns)
    latexdict['tabletype'] = 'table*'
    latexdict['tablefoot'] = ("}\\par\n"
                             )

    float_cols = ['PWV', 'Exposure Time', 'Resolution', 'LAS']

    formats = {key: lambda x: strip_trailing_zeros('{0:0.3f}'.format(round_to_n(x, 2)), nsf=2)
               for key in float_cols}
    formats = {key: lambda x: str(sigfig.round(str(x), sigfigs=2))
               for key in float_cols}

    for sel, nm in ((tm, '12m'), (sm, '7m'), (tp, 'TP')):
        latexdict['header_start'] = f'\\label{{tab:observation_metadata_{nm}}}'#\n\\footnotesize'
        latexdict['preamble'] = f'\\caption{{{nm} Observation Metadata}}\n\\resizebox{{\\textwidth}}{{!}}{{'
        usub_meta[colnames][sel].write(f"{basepath}/papers/continuum_data/observation_metadata_{nm}.tex",
                                       formats=formats,
                                       overwrite=True,
                                       latexdict=latexdict)

    return usub_meta, (tm, sm, tp)


if __name__ == "__main__":

    result = make_latex_table()

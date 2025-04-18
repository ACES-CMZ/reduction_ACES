import os
import numpy as np
import json
from astropy.table import Table, Column
from astropy import table
import requests
import keyring
import datetime
from astropy import units as u
import sigfig
from radio_beam import Beam
from tqdm.auto import tqdm
import casa_formats_io
import radio_beam
from astropy.io import fits

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

    good = (((np.char.find(tbl['filename'], 'preQA3') == -1) &
             (np.char.find(tbl['filename'], 'X15a0_X178') == -1) &
             (np.char.find(tbl['filename'], 'X15a0_Xd6') == -1) &
             (np.char.find(tbl['filename'], 'X15a0_Xe8') == -1) &
             (((np.char.find(tbl['filename'], 'X15a0_X184') >= 0) &
               (np.char.find(tbl['filename'], 'lower_plus_upper') >= 0)) |
              ((np.char.find(tbl['filename'], 'X15a0_X184') == -1))) &
             (np.char.find(tbl['filename'], 'obsolete') == -1) #| False # TODO: check for more
             ) & (
            (np.char.find(tbl['spws'], 'v1') == -1) &
            (np.char.find(tbl['filename'], 'v1.1') == -1) &
            ~(tbl['pbcor']) &
            (tbl['array'] == '12M'))
    ) # noqa

    tbl = tbl[good]
    assert len(tbl) > 0
    assert not any(['v1' in x for x in tbl['spws']])

    # filter out v1 first
    # 'cont' is a special-case for field am; it's... ?!!?
    nueff = {'all': 97.271 * u.GHz, '25,27': 86.539 * u.GHz, '33,35': 99.543 * u.GHz, 'cont': 97.271 * u.GHz}
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


def find_schedblock_elements(data, results=None):
    """
    Recursively search through JSON data to find elements with {'className': 'SchedBlock'}.

    :param data: The JSON object (dict or list) to search.
    :param results: A list to store the results (optional).
    :return: A list of matching elements.
    """
    if results is None:
        results = []

    if isinstance(data, dict):
        # Check if the current dictionary contains the key-value pair
        if data.get('className') == 'SchedBlock':
            results.append(data)
        # Recursively search within each value
        for value in data.values():
            find_schedblock_elements(value, results)

    elif isinstance(data, list):
        # Recursively search within each item in the list
        for item in data:
            find_schedblock_elements(item, results)

    return results


def retrieve_execution_metadata(project_uid='uid://A001/X1525/X290',
                                username='keflavich', timeout=30, retry=5, access_token=None,
                                cache='/red/adamginsburg/ACES/logs/execution_metadata.json'
                                ):

    if cache and os.path.exists(cache):
        with open(cache, 'r') as fh:
            schedblock_meta = json.load(fh)
    else:
        schedblock_meta = {}

    # hard-coded to avoid re-logging-in
    if len(schedblock_meta) == 152:
        return schedblock_meta

    import requests
    if access_token is None:
        from astroquery.alma import Alma
        import time
        Alma.login(username)
        self = Alma._auth
        login_url = f'https://{self.get_valid_host()}{self._LOGIN_ENDPOINT}'

        password = Alma._get_auth_info('keflavich')[1]
        data = {'username': username,
                'password': password,
                'grant_type': self._GRANT_TYPE,
                'client_id': self._CLIENT_ID}

        login_response = self._request('POST', login_url, data=data, cache=False)
        access_token = login_response.json()["access_token"]
        print(login_response.json()['access_token'])
        time.sleep(1)
    project_uidstr = project_uid.replace("/", "%7C")
    resp2 = requests.get(f'https://asa.alma.cl/snoopi/restapi/project/{project_uidstr}',
                         params={'authormode': 'null'}, headers={'Authorization': f'Bearer {access_token}'})
    resp2.raise_for_status()
    project_meta = resp2.json()
    schedblocks = {x['description']: x['archiveUID'] for x in find_schedblock_elements(project_meta)}

    for key, sbuid in tqdm(schedblocks.items()):
        uidstr = sbuid.replace("/", "%7C")

        if key not in schedblock_meta:

            resp2 = requests.get(f'https://asa.alma.cl/snoopi/restapi/schedblock/{uidstr}',
                                 params={'authormode': 'pi'},
                                 headers={'Authorization': f'Bearer {access_token}'},
                                 timeout=timeout
                                 )
            resp2.raise_for_status()
            rslt = resp2.json()
            schedblock_meta[key] = rslt

    if cache:
        with open(cache, 'w') as fh:
            json.dump(schedblock_meta, fh)

    """totaltime, number_pointings, executions"""
    return schedblock_meta


def make_observation_table(access_token=None):
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
    metadata = metadata[metadata['target_name'] == 'Sgr_A_star']

    sub_meta = metadata['schedblock_name', 'Observation Start Time',
                        'Observation End Time', 'pwv', 't_exptime',
                        's_resolution', 'spatial_scale_max', 'member_ous_uid']
    usub_meta = Table(np.unique(sub_meta))

    schedblock_meta = retrieve_execution_metadata(access_token=access_token)

    usub_meta['executions'] = [', '.join([x['date']
                                         for x in schedblock_meta[schedblock_name]['executions']
                                         if x['status'] == 'Pass'
                                         ])
                               for schedblock_name in usub_meta['schedblock_name']]
    usub_meta['npointings'] = [int(schedblock_meta[schedblock_name]['number_pointings'])
                               for schedblock_name in usub_meta['schedblock_name']]
    timefmt = lambda x: f"{x:0.1f}"  # noqa
    usub_meta['time_per_eb'] = [rf'\makecell{{{",\\\\".join([timefmt(float(x['time']))
                                                             for x in schedblock_meta[schedblock_name]['executions']
                                                             if x['status'] == 'Pass'
                                                             ])}}}'
                                for schedblock_name in usub_meta['schedblock_name']]

    # get PWVs
    execution_uids = {schedblock_name: [x['execblockuid'] for x in schedblock_meta[schedblock_name]['executions']
                                        if x['status'] == 'Pass']
                      for schedblock_name in usub_meta['schedblock_name']}
    # schedblock_to_mous = {schedblock_name: schedblock_meta[schedblock_name]['parentObsUnitSetStatusId']
    #                       for schedblock_name in usub_meta['schedblock_name']}
    mous_to_schedblock = {schedblock_meta[schedblock_name]['parentObsUnitSetStatusId']: schedblock_name
                          for schedblock_name in usub_meta['schedblock_name']}

    pwv_cache_fn = f'{basepath}/reduction_ACES/aces/data/pwv_measurements.json'
    if os.path.exists(pwv_cache_fn):
        with open(pwv_cache_fn, 'r') as fh:
            pwv_measurements = json.load(fh)
    else:
        pwv_measurements = {}
    pwvs = []
    pwvfmt = lambda x: f"{x:0.2f}"  # noqa
    for mous in tqdm(usub_meta['member_ous_uid'], desc='PWV'):
        sbname = mous_to_schedblock[mous]
        mousstr = mous.replace('/', '_').replace(':', '_')
        these_pwvs = []
        if sbname in pwv_measurements:
            these_pwvs = pwv_measurements[sbname]
        else:
            for xid in execution_uids[sbname]:
                xidstr = xid.replace("/", "_").replace(":", "_")
                path = f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.{mousstr}/calibrated/working/{xidstr}.ms'
                try:
                    tbl = Table.read(path + "/ASDM_CALWVR")
                    if len(tbl) == 0:
                        print(f"No water table found for {sbname} {mous} {xid} in filename {path}")
                        raise ValueError("No water table found")
                    med_pwv = np.median(tbl['water'])
                    if np.isnan(med_pwv):
                        raise ValueError("WTF?")
                    these_pwvs.append(med_pwv * 1000)
                except ValueError:
                    lltbl = casa_formats_io.table_reader.CASATable.read(path + "/ASDM_CALWVR")
                    #watercol = [x for x in lltbl.column_set.columns if x.name == 'water']
                    tbl = lltbl.as_astropy_table(include_columns=['water'])
                    med_pwv = np.median(tbl['water'])
                    if np.isnan(med_pwv):
                        raise ValueError("WTF?")
                    these_pwvs.append(med_pwv * 1000)
                except IOError as ex:
                    if 'TP' not in sbname:
                        path = f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.{mousstr}/calibrated/{xidstr}.ms'
                        if os.path.exists(path):
                            lltbl = casa_formats_io.table_reader.CASATable.read(path + "/ASDM_CALWVR")
                            tbl = lltbl.as_astropy_table(include_columns=['water'])
                            med_pwv = np.median(tbl['water'])
                            if np.isnan(med_pwv):
                                raise ValueError("WTF?")
                            these_pwvs.append(med_pwv * 1000)
                        else:
                            these_pwvs.append(np.nan)
                            print("No CALWVR table found for ", sbname, mous, xid, ex)
                    print(f"Skipped TP {path}")
                except Exception as ex:
                    these_pwvs.append(np.nan)
                    print("Unknown exception for ", sbname, mous, xid, ex)
            pwv_measurements[sbname] = these_pwvs
        pwvs.append(fr"\makecell{{{',\\\\ '.join(map(pwvfmt, these_pwvs))}}}".replace('nan', '-'))
    with open(pwv_cache_fn, 'w') as fh:
        json.dump(pwv_measurements, fh)
    del usub_meta['pwv']
    usub_meta['pwv'] = pwvs

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
    usub_meta.rename_column('s_resolution', 'Res.')
    usub_meta.rename_column('spatial_scale_max', 'LAS')
    usub_meta.rename_column('config_id', 'Configuration')
    usub_meta.rename_column('executions', 'Execution Dates')
    usub_meta.rename_column('npointings', 'N(P)')
    usub_meta.rename_column('time_per_eb', 'Time')

    updated = (np.char.find(usub_meta['schedblock_name'], 'updated') >= 0)
    print(f"Total rows for tm={tm.sum()} tp={tp.sum()} 7m={sm.sum()}")
    for fieldname in set(usub_meta['Field']):
        for sub in (tm, sm, tp):
            match = (usub_meta['Field'][sub] == fieldname)
            if match.sum() > 1:
                #print(fieldname, match.sum(), sub[sub].sum())
                #sub[sub][match] = updated[sub][match]
                # hard to parse logic, eh?  Just... talk yourself through it 5-10 times...
                sub[sub] &= (~match) | (match & (updated[sub]))
                #print(fieldname, match.sum(), sub[sub].sum())
    print(f"Total rows for tm={tm.sum()} tp={tp.sum()} 7m={sm.sum()} (after filtering for updated)")

    colnames = ['Field', 'Execution Dates', 'Configuration', 'PWV', 'Time', 'N(P)', 'Res.', 'LAS']
    usub_meta['LAS'].unit = u.arcsec
    usub_meta['Res.'].unit = u.arcsec
    #usub_meta['Exposure Time'].unit = u.s
    usub_meta['Time'].unit = u.min
    usub_meta['PWV'].unit = u.mm

    colnames_noconfig = colnames[:2] + colnames[3:]

    usub_meta[colnames][tm].write(f'{basepath}/tables/observation_metadata_12m.ecsv', format='ascii.ecsv', overwrite=True)
    usub_meta[colnames_noconfig][sm].write(f'{basepath}/tables/observation_metadata_7m.ecsv', format='ascii.ecsv', overwrite=True)
    usub_meta[colnames_noconfig][tp].write(f'{basepath}/tables/observation_metadata_TP.ecsv', format='ascii.ecsv', overwrite=True)

    usub_meta['Execution Dates'] = [fr'\makecell{{{",\\\\ ".join([x for x in ed.split(", ")])}.}}' for ed in usub_meta['Execution Dates']]

    #latexdict['header_start'] = '\\label{tab:observation_metadata_12m}'#\n\\footnotesize'
    #latexdict['preamble'] = '\\caption{12m Observation Metadata}\n\\resizebox{\\textwidth}{!}{'
    latexdict['col_align'] = 'l' * len(usub_meta.columns)
    latexdict['tabletype'] = 'table*'
    latexdict['tablefoot'] = (
        "}\\par\n"
        """
        The \\emph{Time} column gives the execution time of each execution block.
        The \\emph{PWV} is the median of the water column in he \\texttt{ASDM\\_CALWVR} table.
        \\emph{Res.} is the resolution in arcseconds.
        \\emph{N(P)} is the number of pointings.
        """
    )

    float_cols = ['Res.', 'LAS']

    formats = {key: lambda x: strip_trailing_zeros('{0:0.3f}'.format(round_to_n(x, 2)), nsf=2)
               for key in float_cols}
    formats = {key: lambda x: str(sigfig.round(str(x), sigfigs=2))
               for key in float_cols}

    for sel, nm in ((tm, '12m'), (sm, '7m'), (tp, 'TP')):
        latexdict['header_start'] = f'\\label{{tab:observation_metadata_{nm}}}'#\n\\footnotesize'
        latexdict['preamble'] = f'\\caption{{{nm} Observation Metadata}}\n\\resizebox{{\\textwidth}}{{!}}{{'
        usub_meta[colnames][sel].write(f"{basepath}/papers/continuum_data/tables/observation_metadata_{nm}.tex",
                                       formats=formats,
                                       overwrite=True,
                                       latexdict=latexdict)

    return usub_meta, (tm, sm, tp)


def make_spw_table():

    linetbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/linelist.csv')
    for fcol in ['Bandwidth', 'F_Lower', 'F_Upper', 'F_Resolution']:
        linetbl[fcol].unit = u.GHz
    linetbl['F_Resolution'] = linetbl['F_Resolution'].to(u.MHz)

    units = {'F_Lower': u.GHz.to_string(u.format.LatexInline),
             'F_Upper': u.GHz.to_string(u.format.LatexInline),
             'F_Resolution': u.MHz.to_string(u.format.LatexInline),
             'Bandwidth': u.GHz.to_string(u.format.LatexInline),
            }
    latexdict['units'] = units

    cols_to_keep = {'12m SPW': "SPW",
                    'F_Lower': r'$\nu_{L}$',
                    'F_Upper': r'$\nu_{U}$',
                    'F_Resolution': r'$\Delta\nu$',
                    'Bandwidth': 'Bandwidth',
                    }
    ftbl = linetbl[list(cols_to_keep.keys())]
    ftbl = Table(np.unique(ftbl))

    for old, new in cols_to_keep.items():
        if old in linetbl.colnames:
            #tbl[old].meta['description'] = description[old]
            ftbl.rename_column(old, new)
            if old in units:
                ftbl[new].unit = units[old]

    float_cols = ['Bandwidth',
                  r'$\nu_{L}$',
                  r'$\nu_{U}$',
                  r'$\Delta\nu$']

    formats = {key: lambda x: strip_trailing_zeros('{0:0.5f}'.format(round_to_n(x, 6)), nsf=4)
               for key in float_cols}
    #formats = {key: lambda x: str(sigfig.round(str(x), sigfigs=2))
    #           for key in float_cols}

    lines = [", ".join([x for x in linetbl['Line'][(linetbl['12m SPW'] == row['SPW']) & (np.char.find(linetbl['col9'], '**') != -1)]])
             for row in ftbl]
    lines = [x.replace("HC3N", "HC$_3$N").replace("13", "$^{13}$")
             .replace(" Alpha", r"$\alpha$").replace("15", "$^{15}$")
             .replace("CH3CHO", "CH$_3$CHO")
             for x in lines]

    linefreqs = [", ".join(
                           [str(x)
                            for x in linetbl['Rest (GHz)'][((linetbl['12m SPW'] == row['SPW']) &
                                                            (np.char.find(linetbl['col9'], '**') != -1))]])
                 for row in ftbl]

    ftbl['Lines'] = [f'{x}\\\\ \n &&&&& {y}' for x, y in zip(lines, linefreqs)]

    # caption needs to be *before* preamble.
    #latexdict['caption'] = 'Continuum Source IDs and photometry'
    latexdict['header_start'] = '\\label{tab:spectral_setup}'#\n\\footnotesize'
    latexdict['preamble'] = '\\caption{ACES Spectral Configuration}\n\\resizebox{\\textwidth}{!}{'
    latexdict['col_align'] = 'l' * len(ftbl.columns)
    latexdict['tabletype'] = 'table*'
    latexdict['tablefoot'] = (
        "}\\par\n"
        "ACES Spectral Configuration, including a non-exhaustive list of prominent, "
        "potentially continuum-affecting, lines.  The included lines are those that are, "
        "in at least some portion of the survey, masked out (see Section \\ref{sec:continuum_selection}).  "
        "The rest frequencies of the targeted lines are given in GHz in the row below their names."
    )

    ftbl.write(f"{basepath}/papers/continuum_data/tables/spectral_setup.tex", formats=formats,
               overwrite=True, latexdict=latexdict)

    return ftbl

    
def make_table_3():
    mosaic_path = f'{basepath}/mosaics/continuum/'
    filenames = {'spw25_27': '12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits',
                 'spw33_35': '12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits',
                 'agg': '12m_continuum_commonbeam_circular_reimaged_mosaic.fits'}
    beams = {key: radio_beam.Beam.from_fits_header(fits.getheader(f'{mosaic_path}/{fn}')).major for key, fn in filenames.items()}
    spw25beam = beams['spw25_27'].to(u.arcsec).value
    spw33beam = beams['spw33_35'].to(u.arcsec).value
    aggbeam = beams['agg'].to(u.arcsec).value

    table3 = fr"""
    \begin{{table}}[h]
    \centering
    \caption{{Continuum Image Products}}
    \begin{{tabular}}{{c c c}}
    \hline
    \hline
    \label{{tab:imagetypes}}
        Image Name                & SPWs & Beam \\
        \hline
        Low, circular & 25, 27 & {spw25beam:0.2f}\arcsec \\
        High, circular & 33, 35 & {spw33beam:0.2f}\arcsec \\
        All data, circular & 25, 27, 33, 35 & {aggbeam:0.2f}\arcsec \\
        All data, best resolution & 25, 27, 33, 35 & n/a \\
        \hline
    \end{{tabular}}
    \end{{table}}
    """ 

    with open(f"{basepath}/papers/continuum_data/tables/table_3.tex", 'w') as fh:
        fh.write(table3)

    return table3


if __name__ == "__main__":

    result = make_latex_table()

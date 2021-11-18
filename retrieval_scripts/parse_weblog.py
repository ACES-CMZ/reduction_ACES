import os
import numpy as np
from astropy import table
from astropy.table import Table,Column
from astropy import units as u
from astropy.utils.console import ProgressBar
from astroquery.alma import Alma
from bs4 import BeautifulSoup
import re

flux_scales = {'Jy': 1,
               'mJy': 1e-3,
               'µJy': 1e-6,
              }

def get_mous_to_sb_mapping(project_code):

    tbl = Alma.query(payload={'project_code': project_code}, cache=False,
                     public=False)['member_ous_uid','schedblock_name', 'qa2_passed']
    mapping = {row['member_ous_uid']: row['schedblock_name'] for row in tbl if row['qa2_passed'] == 'T'}
    return mapping

def get_human_readable_name(weblog, mapping=None):
    print("Reading weblog {0}".format(weblog))
    for directory, dirnames, filenames in os.walk(weblog):
        if 't2-1_details.html' in filenames:
            #print("Found {0}:{1}".format(directory, "t2-1_details.html"))
            with open(os.path.join(directory, 't2-1_details.html')) as fh:
                txt = fh.read()

            max_baseline = re.compile("<th>Max Baseline</th>\s*<td>([0-9a-z\. ]*)</td>").search(txt).groups()[0]
            max_baseline = u.Quantity(max_baseline)

            array_name = ('7MorTP' if max_baseline < 100*u.m else 'TM2'
                          if max_baseline < 1000*u.m else 'TM1')
            #print("array_name = {0}".format(array_name))
            break


    if mapping is None:
        for directory, dirnames, filenames in os.walk(weblog):
            if 't2-2-3.html' in filenames:
                with open(os.path.join(directory, 't2-2-3.html')) as fh:
                    txt = fh.read()
                array_table = table.Table.read(txt, format='ascii.html')
                antenna_size, = map(int, set(array_table['Diameter']))
                break

        for directory, dirnames, filenames in os.walk(weblog):
            if 't2-2-2.html' in filenames:
                with open(os.path.join(directory, 't2-2-2.html')) as fh:
                    txt = fh.read()

                array_table = table.Table.read(txt, format='ascii.html')
                band_string, = set(array_table['Band'])
                band = int(band_string.split()[-1])
                break

        for directory, dirnames, filenames in os.walk(weblog):
            if 't2-2-1.html' in filenames:
                with open(os.path.join(directory, 't2-2-1.html')) as fh:
                    txt = fh.read()

                array_table = table.Table.read(txt, format='ascii.html')
                mask = np.array(['TARGET' in intent for intent in array_table['Intent']], dtype='bool')
                source_name, = set(array_table[mask]['Source Name'])
                break

        if array_name == '7MorTP':
            if antenna_size == 7:
                array_name = '7M'
            elif antenna_size == 12:
                array_name = 'TP'
            else:
                raise

        sbname = "{0}_a_{1:02d}_{2}".format(source_name, band, array_name, )

        print(sbname, max_baseline)

    else:
        for directory, dirnames, filenames in os.walk(weblog):
            if 't1-1.html' in filenames:
                with open(os.path.join(directory, 't1-1.html')) as fh:
                    txt = fh.read()

                soup = BeautifulSoup(txt, 'html5lib')
                overview_tbls = [xx for xx in soup.findAll('table')
                                 if 'summary' in xx.attrs and
                                 xx.attrs['summary'] == 'Data Details']
                assert len(overview_tbls) == 1
                overview_table = overview_tbls[0]

                for row in overview_table.findAll('tr'):
                    if 'OUS Status Entity id' in row.text:
                        for td in row.findAll('td'):
                            if 'uid' in td.text:
                                uid = td.text

                sbname = mapping[uid]
#                try:
#                    sbname = mapping[uid]
#                except:
#                    sbname = 'fail'
#                    print('fail = {0}'.format(directory))

    return sbname, max_baseline

def get_matching_text(list_of_elts, text):
    if hasattr(text, 'search'):
        match = [xx.text for xx in list_of_elts if text.search(xx.text)]
    else:
        match = [xx.text for xx in list_of_elts if text in xx.text]
    if len(match) >= 1:
        return match[0]

def get_calibrator_fluxes(weblog):

    for directory, dirnames, filenames in os.walk(weblog):
        if 't1-1.html' in filenames:
            with open(os.path.join(directory, 't1-1.html')) as fh:
                txt = fh.read()

            soup = BeautifulSoup(txt, 'html5lib')
            date_tbls = [xx for xx in soup.findAll('table')
                         if 'summary' in xx.attrs
                         and xx.attrs['summary'] == 'Measurement Set Summaries']
            assert len(date_tbls) == 1
            date_tbl = date_tbls[0]

            date_map = {}
            for row in date_tbl.findAll('tr'):
                if 'uid___' in row.text:
                    uid = row.find('td').find('a').text
                    date = row.findAll('td')[3].text.split()[0]
                    date_map[uid] = date

        if 't2-4m_details.html' in filenames and 'stage15' in directory:
            with open(os.path.join(directory, 't2-4m_details.html')) as fh:
                txt = fh.read()

            soup = BeautifulSoup(txt, 'html5lib')

            tbls = [xx for xx in soup.findAll('table')
                    if 'summary' in xx.attrs
                    and xx.attrs['summary'] == 'Flux density results']
            if len(tbls) != 1:
                raise ValueError("No flux density data found in pipeline run "
                                 "{0}.".format(weblog))
            tbl = tbls[0]
            rows = tbl.findAll('tr')

            uid, source, freq, spw = None,None,None,None

            data = {}
            for row_a,row_b in zip(rows[3::2],rows[4::2]):
                uid = get_matching_text(row_a.findAll('td'), 'uid') or uid
                source = get_matching_text(row_a.findAll('td'), 'PHASE') or source
                freqstr = get_matching_text(row_a.findAll('td'), 'GHz') or freq
                spw = get_matching_text(row_a.findAll('td'), re.compile('^[0-9][0-9]$')) or spw
                flux_txt = get_matching_text(row_a.findAll('td'), 'Jy')
                catflux_txt = get_matching_text(row_b.findAll('td'), 'Jy')

                assert spw is not None

                fscale = flux_scales[flux_txt.split()[1]]
                efscale = flux_scales[flux_txt.split()[4]]
                cscale = flux_scales[catflux_txt.split()[1]]

                flux = float(flux_txt.split()[0]) * fscale
                eflux = float(flux_txt.split()[3]) * efscale
                catflux = float(catflux_txt.strip().split()[0]) * cscale

                date = date_map[uid]

                freq = float(freqstr.split()[0])
                #freqres = float(freqstr.split()[2])

                data[(source, uid, spw, freq, date)] = {'measured':flux,
                                                        'error': eflux,
                                                        'catalog': catflux}

            return data
    raise ValueError("{0} is not a valid weblog (it may be missing stage15)".format(weblog))

def get_all_fluxes(weblog_list, mapping=None):

    data_dict = {}
    for weblog in ProgressBar(weblog_list):
        try:
            data = get_calibrator_fluxes(weblog)
            name,_ = get_human_readable_name(weblog, mapping=mapping)
            data_dict[name] = data
        except ValueError:
            continue

    flux_data = {name:{ii:
                       {'date': key[4],
                        'ms': key[1],
                        'calibrator': key[0],
                        'spw': key[2],
                        'freq': key[3],
                        'measurement': value}
                       for ii,(key,value) in enumerate(data_.items())
                      }
                 for name,data_ in data_dict.items()
                }

    return flux_data

def fluxes_to_table(flux_dict):

    sbname = Column(name='schedblock_name', data=[name for name,item in flux_dict.items() for row in item])
    uid = Column(name='UID', data=[data['ms'] for name,item in flux_dict.items() for num,data in item.items()])
    calname = Column(name='Calibrator', data=[data['calibrator'] for name,item in flux_dict.items() for num,data in item.items()])
    spw = Column(name='SPW', data=[data['spw'] for name,item in flux_dict.items() for num,data in item.items()])
    date = Column(name='Date', data=[data['date'] for name,item in flux_dict.items() for num,data in item.items()])
    freq = Column(name='Frequency', data=[data['freq'] for name,item in flux_dict.items() for num,data in item.items()])
    flux = Column(name='Flux', data=[data['measurement']['measured'] for name,item in flux_dict.items() for num,data in item.items()])
    eflux = Column(name='Flux error', data=[data['measurement']['error'] for name,item in flux_dict.items() for num,data in item.items()])
    catflux = Column(name='Catalog flux', data=[data['measurement']['catalog'] for name,item in flux_dict.items() for num,data in item.items()])

    tbl = Table([sbname, uid, calname, spw, date, freq, flux, eflux, catflux])

    return tbl


def weblog_names(list_of_weblogs, mapping):

    data = [(get_human_readable_name(weblog, mapping), weblog)
            for weblog in list_of_weblogs]
    hrns = [x[0][0] for x in data]
    if len(set(hrns)) < len(data):
        for nm in set(hrns):
            if hrns.count(nm) > 1:
                print("Fixing {0}".format(nm))
                dupes = [ii for ii,x in enumerate(hrns) if x==nm]
                assert len(dupes) == 2
                bl1 = data[dupes[0]][0][1]
                bl2 = data[dupes[1]][0][1]
                if bl1 < bl2:
                    # short = TM2
                    # long = TM1
                    data[dupes[0]] = \
                        (((data[dupes[0]][0][0].replace('TM1','TM2'),
                           data[dupes[0]][0][1])), data[dupes[0]][1])
                    data[dupes[1]] = \
                        (((data[dupes[1]][0][0].replace('TM2','TM1'),
                           data[dupes[1]][0][1]), data[dupes[1]][1]))
                else:
                    data[dupes[1]] = \
                        (((data[dupes[1]][0][0].replace('TM1','TM2'),
                           data[dupes[1]][0][1]), data[dupes[1]][1]))
                    data[dupes[0]] = \
                        (((data[dupes[0]][0][0].replace('TM2','TM1'),
                           data[dupes[0]][0][1]), data[dupes[0]][1]))

    rslt = {x[0][0]: x[1] for x in data}
    return rslt

def make_links(weblog_maps):
    reverse_map = {v:k for k,v in weblog_maps.items()}
    assert len(reverse_map) == len(weblog_maps)

    for k,v in weblog_maps.items():
        try:
            os.symlink('../{0}'.format(v), 'humanreadable/{0}'.format(k))
        except FileExistsError:
            pass

import numpy as np
import difflib
import os
from ghapi.all import GhApi
import re
import glob

data_dir = '/orange/adamginsburg/ACES/rawdata'

api = GhApi(repo='reduction_ACES', owner='ACES-CMZ')

issues = api('/repos/ACES-CMZ/reduction_ACES/issues')

# uid://A001/X15a0/X17a
uid_re = re.compile("uid://A[0-9]*/X[a-z0-9]*/X[a-z0-9]*")

# Sgr_A_st_ak_03_7M
sb_re = re.compile('Sgr_A_st_([a-z]*)_03_(7M|12M|TP)')


sb_searches = [sb_re.search(issue.title) for issue in issues]
sb_names = [search.group() for search in sb_searches if search]
sb_arrays = {search.group(): search.groups()[1] for search in sb_searches if search}


uid_searches = [uid_re.search(issue.title) for issue in issues]
uid_names = [search.group() for search in uid_searches if search]

uids_to_sbs = {uid: sb for uid, sb in zip(uid_names, sb_names)}
sbs_to_uids = {sb: uid for uid, sb in zip(uid_names, sb_names)}
sbs_to_issues = {sbsearch.group(): issue for issue, sbsearch
                 in zip(issues, sb_searches)
                 if sbsearch}

assert set(sbs_to_issues.keys()) == set(sb_names)


from astroquery.alma import Alma
alma = Alma()
alma.archive_url = 'https://almascience.eso.org'
alma.dataarchive_url = 'https://almascience.eso.org'
results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=False, cache=False)

unique_oids = np.unique(results['obs_id'])
unique_sbs = np.unique(results['schedblock_name'])

assert unique_oids.size == unique_sbs.size

new_oids = set(unique_oids) - set(uid_names)
new_sbs = set(unique_sbs) - set(sb_names)

results.add_index('schedblock_name')

underscore_uid_re = re.compile("uid___A[0-9]*_X[a-z0-9]*_X[a-z0-9]*")
downloaded_uids = [underscore_uid_re.search(x).group()
                   for x in glob.glob(f'{data_dir}/2021.1.00172.L_*_001_of_001.tar')]
weblog_names = [os.path.basename(x) for x in glob.glob('/orange/adamginsburg/web/secure/ACES/weblogs/humanreadable/*')]


sb_status = {}

for new_sb in unique_sbs:
    matches = results.loc[new_sb]
    new_uid = matches['obs_id'][0]
    delivered = '3000' not in matches['obs_release_date'][0]
    if new_sb in sbs_to_uids:
        uuid = sbs_to_uids[new_sb].replace("/","_").replace(":","_")
        downloaded = uuid in downloaded_uids
    else:
        downloaded = False
    weblog_url = f'https://data.rc.ufl.edu/secure/adamginsburg/ACES/weblogs/humanreadable/{new_sb}/html/'
    print(new_sb, downloaded, delivered, weblog_url, new_sb in weblog_names)
    sb_status[new_sb] = {'downloaded': downloaded,
                         'delivered': delivered,
                         'weblog_url': weblog_url}
    if new_sb in new_sbs:
        issuebody = f"""
{new_sb}
[{new_uid}](https://almascience.org/aq/?result_view=observation&mous={new_uid})

* [x] Observations completed?
* [{'x' if delivered else ' '}] Delivered?
* [{'x' if downloaded else ' '}] Downloaded? (specify where)
  * [{'x' if downloaded else ' '}] hipergator
* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked
* [ ] [Weblog]({weblog_url}) Quality Assessment?
* [ ] Imaging: Continuum
* [ ] Imaging: Lines
  * [ ] HCO+
  * [ ] HNCO
  * [ ] H13CN
  * [ ] H13CO+/SiO
  * [ ] Cont 1
  * [ ] Cont 2
        """.replace("\r","")

        array = sb_re.search(new_sb).groups()[1]

        print(f"Posting new issue for {new_sb}")

        title=f"Execution Block ID {new_uid} {new_sb}"

        new_issue = api.issues.create(title=title,
                                      body=issuebody,
                                      labels=['EB', array])
    else:
        issue = sbs_to_issues[new_sb]
        body = issue.body

        assert 'Quality Assessment?' in body

        lines = body.split("\n")
        for ii, line in enumerate(lines):
            if 'Delivered?' in line:
                lines[ii] = f"* [{'x' if delivered else ' '}] Delivered?"
            elif 'Downloaded?' in line:
                lines[ii] = f"* [{'x' if downloaded else ' '}] Downloaded? (specify where)"
                insert_hipergator_at = ii+1
            elif 'Quality Assessment?' in line:
                if 'unpacked' not in body:
                    insert_weblog_at = ii
                else:
                    insert_weblog_at = False
                lines[ii] = f"* [ ] [Weblog]({weblog_url}) Quality Assessment?"
            elif 'unpacked' in line and '[Weblog]' in line:
                lines[ii] = f"* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked"

        # we never want to insert at 0, so it's OK if 0 evaluates to False
        if insert_weblog_at:
            lines.insert(insert_weblog_at, f"* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked")

        if 'hipergator' not in lines[insert_hipergator_at]:
            lines.insert(insert_hipergator_at, f"  * [{'x' if downloaded else ' '}] hipergator")

        issuebody = "\n".join(lines)

        if re.sub('\s', '', issue.body) != re.sub('\s', '', issuebody):
            print(f"Updating issue for {new_sb}")
            if False:
                print('\n'.join(difflib.ndiff(issuebody.split("\n"),
                                              issue.body.split("\n"))
                             ))
            api.issues.update(issue_number=issue.number,
                              title=issue.title,
                              body=issuebody)

# TODO: finish this
# issues = api('/repos/ACES-CMZ/reduction_ACES/issues')
# projects = api('/repos/ACES-CMZ/reduction_ACES/projects')
# columns = api(projects[0].columns_url)
# coldict = {column.name: column for column in columns}
# cards = [api(col.cards_url) for col in columns]
# issue_urls_in_cards = [card.content_url for ccard in cards for card in ccard]
# issue_urls = [issue.url for issue in issues]
# 
# for issue in issues:
#     if 'EB' in [label.name for label in issue.labels]:
#         if issue.url not in issue_urls_in_cards:
#             # need to add it
#             completed = '[x] Observations completed' in issue.body
#             delivered = '[x] Delivered' in issue.body
#             if completed and not delivered:
#                 col = coldict['Completed but not delivered/downloaded']
#             elif completed and delivered:
#                 col = coldict['Delivered Execution Blocks']
#             else:
#                 continue
#             api.projects.create_card(col.id, note=f'{issue.title} {issue.html_url}', content_id=issue.id)
# 
# 

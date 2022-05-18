import numpy as np
import difflib
import os
from ghapi.all import GhApi, paged
import re
import glob

dryrun = os.getenv("DRYRUN") or (dryrun if "dryrun" in locals() else False)
print(f"Dryrun={dryrun}")

def all_flat(apicall, **kwargs):
    paged_stuff = paged(apicall, **kwargs)
    return [x for page in paged_stuff for x in page]

data_dir = '/orange/adamginsburg/ACES/rawdata'

api = GhApi(repo='reduction_ACES', owner='ACES-CMZ')

#paged_issues = paged(api('/repos/ACES-CMZ/reduction_ACES/issues', query={'state': 'all'}))
#paged_issues = paged(api.issues.list_for_repo, state='all')
#issues = [x for page in paged_issues for x in page]
issues = all_flat(api.issues.list_for_repo, state='all')
assert len(issues) > 30

# example: uid://A001/X15a0/X17a
uid_re = re.compile("uid://A[0-9]*/X[a-z0-9]*/X[a-z0-9]*")

# example: Sgr_A_st_ak_03_7M
sb_re = re.compile('Sgr_A_st_([a-z]*)_03_(7M|12M|TP|TM1|TM2).*')


sb_searches = [sb_re.search(issue.title) for issue in issues]
issue_sb_names = [search.group() for search in sb_searches if search]
sb_arrays = {search.group(): search.groups()[1] for search in sb_searches if search}


uid_searches = [uid_re.search(issue.title) for issue in issues]
uid_names = [search.group() for search in uid_searches if search]

assert len(uid_names) == len(issue_sb_names)

uids_to_sbs = {uid: sb for uid, sb in zip(uid_names, issue_sb_names)}
sbs_to_uids = {sb: uid for uid, sb in zip(uid_names, issue_sb_names)}
sbs_to_issues = {sbsearch.group(): issue for issue, sbsearch
                 in zip(issues, sb_searches)
                 if sbsearch}

assert set(sbs_to_issues.keys()) == set(issue_sb_names)


from astroquery.alma import Alma
alma = Alma()
alma.archive_url = 'https://almascience.eso.org'
alma.dataarchive_url = 'https://almascience.eso.org'
results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=None, cache=False)
results.add_index('obs_id')

# .data required b/c making it an index breaks normal usage?
unique_oids = np.unique(results['obs_id'].data)
unique_sbs = np.unique(results['schedblock_name'])
# remember that .unique sorts the data so you can't zip these together!


assert unique_oids.size == unique_sbs.size

new_oids = set(unique_oids) - set(uid_names)
new_sbs = set(unique_sbs) - set(issue_sb_names)

assert len(unique_oids) == len(unique_sbs)

def get_sb_name_from_oid(oid):
    sbnames = np.unique(results.loc[oid]['schedblock_name'])
    assert len(sbnames) == 1
    return sbnames[0]


# add mapping for new (not already in issues) SBs
# explicitly avoid overwriting any existing entries b/c we want the issue names
# to override the archive names
# (order matters since we're updating uids_to_sbs below)
sbs_to_uids.update({get_sb_name_from_oid(uid): uid for uid in unique_oids if uid not in uids_to_sbs})
uids_to_sbs.update({uid: get_sb_name_from_oid(uid) for uid in unique_oids if uid not in uids_to_sbs})

# we know this one is already done, this is a sanity check to make sure we
# don't get caught out by pagination again
assert 'Sgr_A_st_aq_03_7M' not in new_sbs


underscore_uid_re = re.compile("uid___A[0-9]*_X[a-z0-9]*_X[a-z0-9]*")
downloaded_uids = [underscore_uid_re.search(x).group()
                   for x in glob.glob(f'{data_dir}/2021.1.00172.L_*_001_of_001.tar')]
print(f"Downloaded uids = {downloaded_uids}")
weblog_names = [os.path.basename(x) for x in glob.glob('/orange/adamginsburg/web/secure/ACES/weblogs/humanreadable/*')]


sb_status = {}

# loop through oids, not uids: the SB names are _not_ unique, but the UIDs are
for new_oid in unique_oids:
    new_sb_issuename = uids_to_sbs[new_oid]
    new_sb = new_sb_issuename.split(" ")[0]

    need_update = [] # empty list = False
    matches = results.loc[new_oid]
    new_uid = matches['obs_id'][0]
    delivered = '3000' not in matches['obs_release_date'][0]

    uuid = new_oid.replace("/","_").replace(":","_")
    downloaded = uuid in downloaded_uids

    mous = matches['member_ous_uid'][0].replace(":","_").replace("/","_")
    gous = matches['group_ous_uid'][0].replace(":","_").replace("/","_")
    calibrated_dir = f'/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.{gous}/member.{mous}/calibrated'
    if os.path.exists(calibrated_dir):
        mses = (glob.glob(f'{calibrated_dir}/*.ms')
                + glob.glob(f'{calibrated_dir}/*.ms.split.cal'))
    pipeline_run = os.path.exists(calibrated_dir) and len(mses) > 0
    pipeline_links = [x.replace("/orange/adamginsburg/web/secure/", "https://data.rc.ufl.edu/secure/adamginsburg/")
            for x in glob.glob(f"/orange/adamginsburg/web/secure/ACES/weblogs-reimaging/member.{mous}/pipeline*/html/t1-4.html")]

    product_dir = f'/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.{gous}/member.{mous}/product'
    product_filenames = (glob.glob(f"{product_dir}/*Sgr_A_star_sci.spw*.cube.I.pbcor.fits") +
            glob.glob(f"{product_dir}/*Sgr_A_star_sci.spw*.mfs.I.pbcor.fits") +
            glob.glob(f"{product_dir}/*Sgr_A_star_sci.spw*.cont.I.*.fits") +
            glob.glob(f"{product_dir}/*Sgr_A_star_sci.spw*.cube.I.sd.fits")
            )

    reproc_product_dir = f'/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.{gous}/member.{mous}/calibrated/working/'
    reproc_product_filenames = (glob.glob(f"{reproc_product_dir}/*Sgr_A_star_sci.spw*.cont.I.iter1.image.tt0.pbcor.fits") +
            glob.glob(f"{reproc_product_dir}/*Sgr_A_star_sci.spw*.mfs.I.iter1.image.pbcor.fits") +
            glob.glob(f"{reproc_product_dir}/*Sgr_A_star_sci.spw*.cube.I.iter1.image.pbcor.fits"))

    # https://g-76492b.55ba.08cc.data.globus.org/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X114/product/member.uid___A001_X15a0_X114.J1427-4206_bp.spw18.mfs.I.pbcor.fits
    product_links = [f" - [{os.path.basename(fn)}](https://g-76492b.55ba.08cc.data.globus.org/{fn[25:]})" for fn in product_filenames]
    product_link_text = "\n".join(product_links)

    reproc_product_links = [f" - [{os.path.basename(fn)}](https://g-76492b.55ba.08cc.data.globus.org/{fn[25:]})" for fn in reproc_product_filenames]
    reproc_product_link_text = "\n".join(reproc_product_links)

    weblog_url = f'https://data.rc.ufl.edu/secure/adamginsburg/ACES/weblogs/humanreadable/{new_sb.strip().replace(" ","_")}/html/'
    print(f"Operating on sb={new_sb}, oid={new_oid}, dl={downloaded}, delivered={delivered}, url={weblog_url}, weblognames={new_sb in weblog_names}.  pipeline_run={pipeline_run}")
    sb_status[new_sb] = {'downloaded': downloaded,
                         'delivered': delivered,
                         'weblog_url': weblog_url}
    array = sb_re.search(new_sb).groups()[1]

    spw_lines = """
  * [ ] SPW16 H13CN
  * [ ] SPW18 H13CO+/SiO
  * [ ] SPW20 HCO+
  * [ ] SPW22 HNCO
  * [ ] SPW24 Cont 1
  * [ ] SPW26 Cont 2
""" if array in ('7M', 'TP') else """
  * [ ] SPW25 H13CN
  * [ ] SPW27 H13CO+/SiO
  * [ ] SPW29 HCO+
  * [ ] SPW31 HNCO
  * [ ] SPW33 Cont 1
  * [ ] SPW35 Cont 2
"""

    if new_oid in new_oids:
        issuebody = (f"""
{new_sb}
[{new_uid}](https://almascience.org/aq/?result_view=observation&mous={new_uid})

* [x] Observations completed?
* [{'x' if delivered else ' '}] Delivered?
* [{'x' if downloaded else ' '}] Downloaded? (specify where)
  * [{'x' if downloaded else ' '}] hipergator""" + ("" if array == 'TP' else f"""
  * [{'x' if pipeline_run else ' '}] hipergator pipeline run
""") +
f"""
* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked
* [ ] [Weblog]({weblog_url}) Quality Assessment?
* [ ] Imaging: Continuum
* [ ] Imaging: Lines
{spw_lines}

## Product Links:

{product_link_text}

## Reprocessed Product Links:

{reproc_product_link_text}
""".replace("\r",""))


        print(f"Posting new issue for {new_sb} -> {new_sb_issuename}")

        title=f"Execution Block ID {new_uid} {new_sb_issuename}"

        labels = ['EB', array]
        if delivered:
            labels.append('Delivered')

        new_issue = api.issues.create(title=title,
                                      body=issuebody,
                                      labels=labels)
    else:
        #print(f"Issue exists: Possibly updating existing issue {new_sb_issuename}")
        issue = sbs_to_issues[new_sb_issuename]
        body = issue.body
        labels = [lb.name for lb in issue.labels]

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
                #don't overwrite!  this would un-check!
                # (probably this was added earlier when some versions were missing the link URL, but it's not now)
                # lines[ii] = f"* [ ] [Weblog]({weblog_url}) Quality Assessment?"
            elif 'unpacked' in line and '[Weblog]' in line:
                lines[ii] = f"* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked"

        # we never want to insert at 0, so it's OK if 0 evaluates to False
        if insert_weblog_at:
            lines.insert(insert_weblog_at, f"* [{'x' if new_sb in weblog_names else ' '}] [Weblog]({weblog_url}) unpacked")

        if 'hipergator' not in lines[insert_hipergator_at]:
            need_update.append("Downloaded")
            lines.insert(insert_hipergator_at, f"  * [{'x' if downloaded else ' '}] hipergator")

        pipeline_linenumber = insert_hipergator_at+1
        pipeline_line_text = f"  * [{'x' if pipeline_run else ' '}] hipergator pipeline run"
        if 'hipergator pipeline run' not in lines[pipeline_linenumber] and array != "TP":
            need_update.append("Pipelined")
            lines.insert(pipeline_linenumber, pipeline_line_text)
        elif 'hipergator pipeline run' in lines[pipeline_linenumber] and array != "TP" and lines[pipeline_linenumber] != pipeline_line_text:
            need_update.append("Pipelined")
            lines[pipeline_linenumber] = pipeline_line_text

        pipeline_links_linenumber = pipeline_linenumber+1
        # check if the link is included at all
        pipeline_links_done = [link in body for link in pipeline_links]
        if len(pipeline_links) > sum(pipeline_links_done):
            need_update.append("New pipeline run links")
            for link,done in zip(pipeline_links, pipeline_links_done):
                if not done:
                    assert 'https' in link
                    pipenumber = link.split("/")[-3].split("-")[1]
                    plink_text = f"    * [Pipeline Run {pipenumber}]({link})"
                    lines.insert(pipeline_links_linenumber, plink_text)


        issuebody = "\n".join(lines)

        if delivered and 'Delivered' not in labels:
            labels.append('Delivered')
            need_update.append('Delivered')

        if re.sub('\s', '', issue.body) != re.sub('\s', '', issuebody):
            need_update.append("Generic: something changed")


        if product_link_text:
            productlinks = f"""

## Product Links:

{product_link_text}
""".replace("\r","")

            if '## Product Links:' not in issue.body:
                need_update.append("New product links - none were present")
                issebody += productlinks
            elif issue.body.strip().endswith("## Product Links:"):
                need_update.append("New product links - end was blank")
                issuebody = issuebody.strip().split("## Product Links:")[0] + productlinks
            elif "## Product Links:\n\n\n\n\n##" in issue.body:
                need_update.append("Update product links")
                issuebody = issuebody.replace("## Product Links:\n\n\n\n\n", productlinks)


        if reproc_product_link_text:
            reproductlinks = f"""

## Reprocessed Product Links:

{reproc_product_link_text}
""".replace("\r","")

            if '## Reprocessed Product Links:' not in issue.body:
                need_update.append("New reproduct links")
                issuebody += reproductlinks
            elif "## Reprocessed Product Links:\n\n\n" in issue.body:
                need_update.append("Update reproduct links")
                issuebody = issuebody.replace("## Reprocessed Product Links:\n\n\n", reproductlinks)

        if need_update:
            print(f"Updating issue for {new_sb} -> {new_sb_issuename}.  need_update={need_update}")
            if False:
                print('\n'.join(difflib.ndiff(issuebody.split("\n"),
                                              issue.body.split("\n"))
                             ))


            if not dryrun:
                api.issues.update(issue_number=issue.number,
                                  title=issue.title,
                                  body=issuebody,
                                  labels=labels)

            # use this to break raise Exception("Completed a run; check it")


paged_issues = paged(api.issues.list_for_repo, state='all')
issues = [x for page in paged_issues for x in page]
assert len(issues) > 30

# should only ever be 1 project, so pagination not needed
projects = api.projects.list_for_repo()

columns = api.projects.list_columns(projects[0].id) #api(projects[0].columns_url)
coldict = {column.name: column for column in columns}
#cards = [api(col.cards_url) for col in columns]
cards = [x for col in columns for x in all_flat(api.projects.list_cards, column_id=col.id)]

issue_urls_in_cards = [card.content_url for card in cards]
issue_urls = [issue.url for issue in issues]

for issue in issues:
    if issue.state == 'closed':
        # skip closed issues
        continue

    print(".", end="")
    if 'EB' in [label.name for label in issue.labels]:
        completed = '[x] Observations completed' in issue.body
        delivered = '[x] Delivered' in issue.body

        if issue.url not in issue_urls_in_cards:
            # need to add it
            if completed and not delivered:
                col = coldict['Completed but not delivered/downloaded']
            elif completed and delivered:
                col = coldict['Delivered Execution Blocks']
            else:
                continue

            print(f"Adding issue {issue.title} to the {'Completed' if completed else ''}{'Delivered' if delivered else ''} column")
            if not dryrun:
                api(path=f'/projects/columns/{col.id}/cards', verb='POST',
                    data={'content_id': issue.id,
                          'column_id': col.id,
                          'content_type': 'Issue'
                         },
                         headers={"Accept": "application/vnd.github.v3+json"}
                      )
        else:
            # check if issue is categorized right

            completed_not_delivered = api(coldict['Completed but not delivered/downloaded'].cards_url)
            completed_and_delivered = api(coldict['Delivered Execution Blocks'].cards_url)
            other = [api(coldict[key].cards_url) for key in coldict if key not in ['Completed but not delivered/downloaded', 'Delivered Execution Blocks']]

            if completed and delivered and issue.url not in [card.content_url for ccard in other for card in ccard]:
                if issue.url not in [card.content_url for card in completed_and_delivered]:
                    # move issue to completed_and_delivered
                    current_card = [card for card in completed_not_delivered if card.content_url == issue.url][0]
                    print(f"MOVING issue {issue.title}")

                    if not dryrun:
                        # remove it from current location
                        api(path=f'/projects/columns/cards/{current_card.id}',
                            verb='DELETE', headers={"Accept":
                                                    "application/vnd.github.v3+json"})

                        # add it to new location
                        CADid = coldict['Delivered Execution Blocks'].id
                        api(path=f'/projects/columns/{CADid}/cards', verb='POST',
                            data={'content_id': issue.id,
                                  'column_id': CADid,
                                  'content_type': 'Issue'
                                 },
                                 headers={"Accept": "application/vnd.github.v3+json"}
                              )


import shutil
import glob
import numpy as np
import tarfile
from astroquery.alma import Alma
import os

alma = Alma()
alma.cache_location = Alma.cache_location = '.'

# try the ESO servers?
#alma.archive_url = 'https://almascience.eso.org'

username = input("Username: ")
alma.login(username)

results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=False, cache=False)

existing_tarballs = glob.glob("2021.1.00172.L/weblog_tarballs/*weblog.tgz")
mouses = results['obs_id']
mouses_filtered = [x for x in mouses
                   if not any([x[6:].replace("/","_") in y
                               for y in existing_tarballs])]

print("Found {0} new files out of {1}".format(len(mouses_filtered), len(mouses)))

files = alma.stage_data(mouses_filtered)

product_mask = np.array(['asdm' not in row['URL'] for row in files])

product = files[product_mask]
print("Found {0} total product files to download".format(product))
weblog_mask = np.array(['weblog.tgz' in row['URL'] for row in product], dtype='bool')
weblog_files = product[weblog_mask]
print("Found {0} weblogs to download".format(len(weblog_files)))
weblog_fns = [x.split("/")[-1] for x in weblog_files['URL']]
existing_weblog_fns = [x.split("/")[-1] for x in existing_tarballs]
weblog_urls_to_download = [row['URL'] for row,fn in zip(weblog_files, weblog_fns) if fn not in existing_weblog_fns]
print("Found {0} *new* weblogs to download".format(len(weblog_urls_to_download)))

for fn in weblog_urls_to_download:
    if 'tgz' not in fn:
        raise ValueError

weblog_tarballs = alma.download_files(weblog_urls_to_download)

if not os.path.exists('2021.1.00172.L'):
    os.mkdir('2021.1.00172.L')
if not os.path.exists('2021.1.00172.L/weblog_tarballs'):
    os.mkdir('2021.1.00172.L/weblog_tarballs')

#weblogs = weblogs_band3+weblogs_band6
for logfile in weblog_tarballs:
    print(logfile)
    with tarfile.open(logfile) as tf:
        tf.extractall('2021.1.00172.L')

for dirpath, dirnames, filenames in os.walk('.'):
    for fn in filenames:
        if "weblog.tgz" in fn:
            shutil.move(os.path.join(dirpath, fn),
                        os.path.join('2021.1.00172.L/weblog_tarballs', fn))

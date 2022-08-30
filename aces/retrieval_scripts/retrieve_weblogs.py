import shutil
import glob
import numpy as np
import tarfile
from astroquery.alma import Alma
import os
import six
import sys

from aces import conf

basepath = conf.basepath


def main():
    if len(sys.argv) > 1:
        username = sys.argv[1]
        print(f"Using username {username}")
    else:
        username = six.moves.input("Username: ")

    alma = Alma()
    alma.cache_location = Alma.cache_location = '.'

    # try the ESO servers?
    alma.archive_url = 'https://almascience.eso.org'
    alma.dataarchive_url = 'https://almascience.eso.org'

    alma.login(username)
    print(f"Logged in as {username}.  Performing weblog query.", flush=True)

    results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=None)

    release_dates = results['obs_release_date']
    # remove all 3000- dates
    ok_release_dates = np.array([int(x[0]) < 3 for x in release_dates])
    results = results[ok_release_dates]

    existing_tarballs = glob.glob(f"{basepath}/data/2021.1.00172.L/weblog_tarballs/*weblog.tgz")
    mouses = results['member_ous_uid']
    mouses_filtered = [x for x in mouses
                       if not any([x[6:].replace("/", "_") in y
                                   for y in existing_tarballs])]

    print("Found {0} new files out of {1}".format(len(mouses_filtered), len(mouses)))
    if len(mouses_filtered) > 0:

        # files = alma.stage_data(mouses_filtered)
        files = alma.get_data_info(mouses_filtered, expand_tarfiles=True)

        print(f"Found {len(files)} files using get_data_info")

        product_mask = np.array(['asdm' not in row['access_url'] for row in files])
        print(f"There are {product_mask.sum()} products")

        product = files[product_mask]
        print("Found {0} total product files to download".format(product))
        weblog_mask = np.array(['weblog.tgz' in row['access_url'] for row in product], dtype='bool')
        weblog_files = product[weblog_mask]
        print("Found {0} weblogs to download".format(len(weblog_files)))
        weblog_fns = [x.split("/")[-1] for x in weblog_files['access_url']]
        existing_weblog_fns = [x.split("/")[-1] for x in existing_tarballs]
        weblog_urls_to_download = [row['access_url'] for row, fn in zip(weblog_files, weblog_fns) if fn not in existing_weblog_fns]
        print("Found {0} *new* weblogs to download".format(len(weblog_urls_to_download)))

        for fn in weblog_urls_to_download:
            if 'tgz' not in fn:
                raise ValueError

        weblog_tarballs = alma.download_files(weblog_urls_to_download)

        if not os.path.exists(f'{basepath}/data/2021.1.00172.L'):
            os.mkdir(f'{basepath}/data/2021.1.00172.L')
        if not os.path.exists(f'{basepath}/data/2021.1.00172.L/weblog_tarballs'):
            os.mkdir(f'{basepath}/data/2021.1.00172.L/weblog_tarballs')
        if not os.path.exists(f'{basepath}/data/2021.1.00172.L/weblogs'):
            os.mkdir(f'{basepath}/data/2021.1.00172.L/weblogs')

        # weblogs = weblogs_band3+weblogs_band6
        for logfile in weblog_tarballs:
            print(logfile)
            with tarfile.open(logfile) as tf:
                tf.extractall(f'{basepath}/data/2021.1.00172.L/weblogs')

        for dirpath, dirnames, filenames in os.walk('.'):
            for fn in filenames:
                if "weblog.tgz" in fn:
                    shutil.move(os.path.join(dirpath, fn),
                                os.path.join(f'{basepath}/data/2021.1.00172.L/weblog_tarballs', fn))

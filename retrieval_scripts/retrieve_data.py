import numpy as np
from astroquery.alma import Alma
import six
import os
import sys

if len(sys.argv) > 1:
    username = sys.argv[1]
    print(f"Using username {username}")
else:
    username = six.moves.input("Username: ")

alma = Alma()
alma.archive_url = 'https://almascience.eso.org'
alma.dataarchive_url = 'https://almascience.eso.org'
alma.cache_location = Alma.cache_location = '.'
alma.login(username)

results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=False, cache=False)

obsids = np.unique(results['obs_id'])

data = alma.retrieve_data_from_uid(obsids)


# Optional: %run retrieve_data <username> True will extract files
if len(sys.argv) > 2:
    extract = bool(sys.argv[2])
else:
    extract = False

if extract:
    import tarfile
    for fn in data:
        if fn.endswith('.tar'):
            with tarfile.TarFile(fn) as tf:
                for member in tf.getmembers():
                    if not os.path.exists(member.name):
                        tf.extract(member)

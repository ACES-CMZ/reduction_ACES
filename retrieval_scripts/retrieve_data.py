import numpy as np
from astroquery.alma import Alma
import six
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

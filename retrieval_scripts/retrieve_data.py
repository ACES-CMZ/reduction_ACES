from astroquery.alma import Alma
import six

alma = Alma()
alma.cache_location = Alma.cache_location = '.'
username = six.moves.input("Username: ")
alma.login(username)

results = alma.query(payload=dict(project_code='2021.1.00172.L'), public=False, cache=False)

obsids = np.unique(results['obs_id'])

data = alma.retrieve_data_from_uid(obsids)

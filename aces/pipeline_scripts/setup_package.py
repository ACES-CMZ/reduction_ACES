import os


def get_package_data():
    paths = [os.path.join('data', '*.json'), ]
    return {'aces.pipeline_scripts': paths}
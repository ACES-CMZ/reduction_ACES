import os


def get_package_data():
    paths = [os.path.join('.', '*.json'), ]
    return {'aces.pipeline_scripts': paths}
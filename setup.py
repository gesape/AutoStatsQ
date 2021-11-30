# import sys

# from setuptools import setup, find_packages
from setuptools import setup

version = '0.1'

setup(
    name='autostatsq',
    version=version,
    author='Gesa Petersen',
    author_email='gesap@gfz-potsdam.de',
    python_requires='!=3.0.*, !=3.1.*, !=3.2.*, <4',
    install_requires=[],
    packages=['autostatsq', 'autostatsq.generate_report'],
    package_dir={'autostatsq': 'src', 'autostatsq.generate_report': 'src/generate_report'},
    entry_points={
        'console_scripts': [
            'autostatsq = autostatsq.network_control:main',
        ]
    },
    include_package_data=True,
    package_data={'': ['generate_report/index_template.html', 
                       'generate_report/figures/*',
                       'generate_report/reveal/',
                       'generate_report/reveal/*',
                       'generate_report/reveal/*/*',
                       'generate_report/reveal/*/*/*',
                       'generate_report/reveal/*/*/*/*',
                       'generate_report/theme/*',
                       'generate_report/theme/*/*',
                       'generate_report/theme/*/*/*',
                       'generate_report/theme/*/*/*/*' ]}
)

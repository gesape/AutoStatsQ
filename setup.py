from setuptools import setup

version = '0.1.1'

dependencies = [
    "numpy",
    "pyrocko",
    "grond",
    "matplotlib"
]

setup(
    name='autostatsq',
    version=version,
    author='Gesa Petersen',
    author_email='gesap@gfz-potsdam.de',
    python_requires='!=3.0.*, !=3.1.*, !=3.2.*, <4',
    install_requires=dependencies,
    packages=['autostatsq'],
    package_dir={'autostatsq': 'src'},
    entry_points={
        'console_scripts': [
            'autostatsq = autostatsq.network_control:main',
        ]
    }
)

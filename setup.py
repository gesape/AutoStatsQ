# import sys

# from setuptools import setup, find_packages
from setuptools import setup

version = "0.1"

setup(
    name="autostatsq",
    version=version,
    author="Gesa Petersen",
    author_email="gesap@gfz-potsdam.de",
    python_requires="!=3.0.*, !=3.1.*, !=3.2.*, <4",
    install_requires=[],
    packages=["autostatsq"],
    package_dir={"autostatsq": "src"},
    entry_points={
        "console_scripts": [
            "autostatsq = autostatsq.network_control:main",
        ]
    },
)

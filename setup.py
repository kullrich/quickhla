#!/usr/bin/python
# -*- coding: UTF-8 -*-

import sys
import os
from setuptools import setup

if sys.version_info < (3, 8):
    sys.exit('\nERROR: quickhla requires python 3.8 or greater')

__version__ = open(os.path.join('quickhla', 'VERSION')).read().strip()


def main():
    setup(
        name="quickhla",
        version=__version__,
        description="classifies HLA genes",
        long_description=open("README.md").read(),
        long_description_content_type="text/markdown",
        author="Kristian Ullrich",
        author_email="ullrich@evolbio.mpg.de",
        packages=["quickhla"],
        keywords=["hla", "classifier"],
        url="https://gitlab.gwdg.de/kristian.ullrich/quickhla",
        license="MIT",
        classifiers=[
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3.8",
        ],
        install_requires=['argparse',
                          'numpy',
                          'biopython'],
        include_package_data=True,
        zip_safe=False,
        entry_points={
            "console_scripts": [
                "quickhla=quickhla.quickhla:main",
            ]
        },
    )


if __name__ == "__main__":
    main()

#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import msnpy


def main():

    if sys.version_info[0] != 2 and sys.version_info[1] <= 7:
        sys.exit("Python-2.7.8 is required ")

    setuptools.setup(name="msnpy",
        version=msnpy.__version__,
        description="Python package for data processing of direct-infusion mass spectrometry-based metabolomics and lipidomics data",
        long_description=open('README.rst').read(),
        author="Ralf Weber",
        author_email="r.j.weber@bham.ac.uk",
        url="https://github.com/computational-metabolomics/msnpy",
        license="GPLv3",
        platforms=['Windows, UNIX'],
        keywords=['Metabolomics', 'Lipidomics', 'Mass spectrometry', 'Data Processing', 'Direct-Infusion Mass Spectrometry'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        install_requires=open('requirements.txt').read().splitlines(),
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 2",
          "Programming Language :: Python :: 2.7",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
        ],
        entry_points={
         'console_scripts': [
             'msnpy = msnpy.__main__:main'
         ]
        }
    )


if __name__ == "__main__":
    main()

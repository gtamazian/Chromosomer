#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from setuptools import find_packages, setup

import chromosomer

setup(
    name="chromosomer",
    version=chromosomer.__version__,
    description="Reference-assisted chromosome assembly tool",
    author="Gaik Tamazian",
    author_email="gaik.tamazian@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(exclude=["doc", "tests*"]),
    install_requires=["pyfaidx", "future", "bioformats"],
    entry_points={"console_scripts": ["chromosomer = chromosomer.cli:chromosomer"]},
)

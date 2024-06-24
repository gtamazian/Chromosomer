#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import chromosomer.cli
import os
import shutil
import sys
import tempfile
import unittest


class TestFragmentSimulator(unittest.TestCase):
    def setUp(self):
        self.__output_dir = tempfile.NamedTemporaryFile().name
        os.mkdir(self.__output_dir)

    def test_simulator(self):
        sys.argv = [
            "",
            "simulator",
            "-g",
            "10",
            "-p",
            "5",
            "--prefix",
            "test_",
            "20",
            "100",
            "5",
            self.__output_dir,
        ]
        chromosomer.cli.chromosomer()

    def tearDown(self):
        if os.path.isdir(self.__output_dir):
            shutil.rmtree(self.__output_dir)

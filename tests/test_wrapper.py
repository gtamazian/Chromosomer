#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import tempfile
import unittest

from chromosomer.shims import FastaWriter, RandomSequence
from chromosomer.wrapper.blast import BlastN, MakeBlastDb

path = os.path.dirname(__file__)
os.chdir(path)


class TestWrapperMakeBlastDb(unittest.TestCase):
    def setUp(self):
        self.__fasta = tempfile.mkstemp()[1]
        self.__seq_length = 100
        self.__seq_number = 10
        self.__dbname = os.path.join(os.path.split(self.__fasta)[0], "test")

    def tearDown(self):
        temp_files = glob.glob("{}*".format(self.__dbname))
        for i in temp_files:
            os.unlink(i)

    def test_launch(self):
        """
        The the makeblastdb launching routine.
        """
        with FastaWriter(self.__fasta) as fasta_writer:
            for i in range(self.__seq_number):
                seq_generator = RandomSequence(self.__seq_length)
                fasta_writer.write("seq{}".format(i + 1), seq_generator.get())

        wrapper = MakeBlastDb(self.__fasta, self.__dbname)
        wrapper.launch()


class TestWrapperBlastN(unittest.TestCase):
    def setUp(self):
        # create a FASTA file of random sequences and a BLAST
        # database from it
        self.__fasta = tempfile.mkstemp()[1]
        self.__output = tempfile.mkstemp()[1]
        self.__seq_length = 100
        self.__seq_number = 10

        with FastaWriter(self.__fasta) as fasta_writer:
            for i in range(self.__seq_number):
                seq_generator = RandomSequence(self.__seq_length)
                fasta_writer.write("seq{}".format(i + 1), seq_generator.get())

        wrapper = MakeBlastDb(self.__fasta)
        wrapper.launch()

    def tearDown(self):
        temp_files = glob.glob("{}*".format(self.__fasta))
        for i in temp_files:
            os.unlink(i)
        if os.path.isfile(self.__output):
            os.unlink(self.__output)

    def test_launch(self):
        """
        Test the blastn testing routine.
        """
        wrapper = BlastN(self.__fasta, self.__fasta, self.__output)
        wrapper.set("-outfmt", 6)
        wrapper.get("-outfmt")
        wrapper.launch()

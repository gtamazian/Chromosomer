#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import logging
import os
import random
import tempfile
import unittest

import pyfaidx

from chromosomer.fragment import (
    AlignmentToMap,
    AlignmentToMapError,
    Map,
    MapError,
    SeqLengths,
    Simulator,
)
from chromosomer.shims import BedReader, BlastTab, FastaWriter, RandomSequence
from chromosomer.wrapper.blast import BlastN, MakeBlastDb

path = os.path.dirname(__file__)
os.chdir(path)


class TestFragmentMap(unittest.TestCase):
    def setUp(self):
        self.__test_line = os.path.join("data", "fragment_map", "fragment_map_line.txt")
        self.__incorrect_file_dir = os.path.join(
            "data", "fragment_map", "incorrect_input"
        )
        self.__output_dir = os.path.join("data", "fragment_map")
        self.__output_file = tempfile.NamedTemporaryFile().name
        self.__incorrect_files = os.listdir(self.__incorrect_file_dir)
        # silence the logging message
        logging.disable(logging.ERROR)

    def test_add_record(self):
        """
        Check if fragment records are added correctly.
        """
        fragment_map = Map()
        new_record = Map.Record(
            fr_name="fragment1",
            fr_length=180,
            fr_start=0,
            fr_end=180,
            fr_strand="+",
            ref_chr="chr1",
            ref_start=5000,
            ref_end=5180,
        )
        fragment_map.add_record(new_record)

    def test_read(self):
        """
        Test the Map reading routine.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        fragment = next(fragment_map.fragments("chr1"))
        self.assertEqual(fragment.fr_name, "fragment1")
        self.assertEqual(fragment.fr_length, 180)
        self.assertEqual(fragment.fr_start, 0)
        self.assertEqual(fragment.fr_end, 180)
        self.assertEqual(fragment.fr_strand, "+")
        self.assertEqual(fragment.ref_chr, "chr1")
        self.assertEqual(fragment.ref_start, 5000)
        self.assertEqual(fragment.ref_end, 5180)

        # check for incorrect input files
        for i in self.__incorrect_files:
            with self.assertRaises(MapError):
                fragment_map.read(os.path.join(self.__incorrect_file_dir, i))

    def test_chromosomes(self):
        """
        Test the Map chromosomes iterator.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        chromosomes = list(fragment_map.chromosomes())
        self.assertEqual(chromosomes, ["chr1"])

    def test_fragments(self):
        """
        Test the Map fragments iterator.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        fragments = list(fragment_map.fragments("chr1"))
        self.assertEqual(len(fragments), 1)
        self.assertIsInstance(fragments[0], Map.Record)

        # check if the missing chromosome is processed correctly
        with self.assertRaises(MapError):
            list(fragment_map.fragments("chrN"))

    def test_summary(self):
        """
        Test the Map summary routine.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        self.assertIsInstance(fragment_map.summary(), dict)

    def test_convert2bed(self):
        """
        Test the BED conversion routine.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        fragment_map.convert2bed(self.__output_file)
        # try to read the produced BED file
        with open(self.__output_file) as bed_file:
            reader = BedReader(bed_file)
            for _ in reader.records():
                pass

    def test_write(self):
        """
        Test the Map writing routine.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)

        output_filename = os.path.join(
            "data", "fragment_map", "fragment_map_output.txt"
        )
        fragment_map.write(output_filename)

        with open(output_filename) as output_file:
            with open(self.__test_line) as original_file:
                for x, y in zip(original_file, output_file):
                    self.assertEqual(x, y)

        os.unlink(output_filename)

    def test_assemble(self):
        """
        Test the assemble routine.
        """
        # first, we form fragment and chromosome sequences
        fragments = {}
        fragment_pattern = ["AC", "AG", "CT", "CG", "AT"]
        for i, pattern in enumerate(fragment_pattern):
            fragments["fragment{}".format(i + 1)] = pattern * 5
        # a negative number indicated reverse orientation of a fragment
        chromosome_content = {"chr1": [1, -2, 3], "chr2": [-4, 5]}
        # get chromosome sequences
        chromosomes = {}
        complement = str.maketrans("ATCGatcgNnXx", "TAGCtagcNnXx")
        gap_size = 10
        for i, chromosome_fragments in chromosome_content.items():
            chromosomes[i] = []
            for j in chromosome_fragments:
                fr_seq = fragments["fragment{}".format(abs(j))]
                if j < 0:
                    chromosomes[i].append(fr_seq[::-1].translate(complement))
                else:
                    chromosomes[i].append(fr_seq)
                chromosomes[i].append("N" * gap_size)
            chromosomes[i] = "".join(chromosomes[i])
        # contruct a fragment __map
        fragment_map = Map()
        for i, chromosome_fragments in chromosome_content.items():
            current_start = 0
            for j in chromosome_fragments:
                fr_name = "fragment{}".format(abs(j))
                fr_length = 10
                fr_start = 0
                fr_end = fr_length
                fr_strand = "+" if j > 0 else "-"
                ref_chr = i
                ref_start = current_start
                ref_end = current_start + fr_length
                fragment_map.add_record(
                    Map.Record(
                        fr_name,
                        fr_length,
                        fr_start,
                        fr_end,
                        fr_strand,
                        ref_chr,
                        ref_start,
                        ref_end,
                    )
                )
                current_start += fr_length
                # add the gap
                fr_name = "GAP"
                fr_length = gap_size
                fr_start = 0
                fr_end = gap_size
                fr_strand = "+"
                ref_chr = i
                ref_start = current_start
                ref_end = current_start + fr_end
                fragment_map.add_record(
                    Map.Record(
                        fr_name,
                        fr_length,
                        fr_start,
                        fr_end,
                        fr_strand,
                        ref_chr,
                        ref_start,
                        ref_end,
                    )
                )
                current_start += fr_length

        output_chromosomes = os.path.join(self.__output_dir, "temp_chromosomes.txt")
        output_fragments = os.path.join(self.__output_dir, "temp_fragments.txt")

        # write the fragment sequences to a FASTA file
        with FastaWriter(output_fragments) as writer:
            for i, j in fragments.items():
                writer.write(i, j)

        fragment_map.assemble(output_fragments, output_chromosomes)

        # read fragments from the written FASTA file and compare them
        # to the original ones
        assembled_chromosomes = pyfaidx.Fasta(output_chromosomes)
        for i, seq in chromosomes.items():
            self.assertEqual(seq, assembled_chromosomes[i][:].seq)

        # try to use the fragment absent in the FASTA file of
        # fragment sequences
        fragment_map.add_record(
            Map.Record(
                fr_name="missing_fragment",
                fr_length=0,
                fr_start=0,
                fr_end=0,
                fr_strand="+",
                ref_chr="chr3",
                ref_start=0,
                ref_end=0,
            )
        )
        with self.assertRaises(MapError):
            fragment_map.assemble(output_fragments, output_chromosomes)

        os.unlink(output_chromosomes)
        os.unlink(output_chromosomes + ".fai")
        os.unlink(output_fragments)
        os.unlink(output_fragments + ".fai")

    def tearDown(self):
        if os.path.isfile(self.__output_file):
            os.unlink(self.__output_file)


class TestFragmentLength(unittest.TestCase):
    def setUp(self):
        self.__fragment_number = 10
        self.__fragment_length = 10
        self.__fasta_temp = tempfile.mkstemp()[1]

        # create a FASTA file of random sequences
        with FastaWriter(self.__fasta_temp) as fasta_writer:
            seq_generator = RandomSequence(self.__fragment_length)
            for i in range(self.__fragment_length):
                fasta_writer.write("seq{}".format(i + 1), seq_generator.get())

    def test_lengths(self):
        """
        Test the lengths method.
        """
        x = SeqLengths(self.__fasta_temp)
        lengths = x.lengths()
        for i in lengths.values():
            self.assertEqual(i, self.__fragment_length)

    def tearDown(self):
        os.unlink(self.__fasta_temp)


class TestFragmentSimulator(unittest.TestCase):
    def setUp(self):
        self.__fragment_number = 10
        self.__chromosome_number = 2
        self.__fragment_length = 10
        self.__unplaced_number = 2
        self.__gap_size = 5

        self.__simulator = Simulator(
            self.__fragment_length,
            self.__fragment_number,
            self.__chromosome_number,
            self.__unplaced_number,
            self.__gap_size,
        )

    def test_write(self):
        """
        The the writing method of the fragment simulator.
        """
        self.__fragments = tempfile.mkstemp()[1]
        self.__chromosomes = tempfile.mkstemp()[1]
        self.__map = tempfile.mkstemp()[1]

        self.__simulator.write(self.__map, self.__fragments, self.__chromosomes)

        # check if the correct number of fragment and chromosome
        # sequences was written
        fragment_fasta = pyfaidx.Fasta(self.__fragments)
        self.assertEqual(
            len(list(fragment_fasta.keys())),
            self.__fragment_number + self.__unplaced_number,
        )
        chromosome_fasta = pyfaidx.Fasta(self.__chromosomes)
        self.assertEqual(len(list(chromosome_fasta.keys())), self.__chromosome_number)

        # check if a correct fragment map was written
        test_map = Map()
        test_map.read(self.__map)

        os.unlink(self.__fragments)
        os.unlink(self.__fragments + ".fai")
        os.unlink(self.__chromosomes)
        os.unlink(self.__chromosomes + ".fai")
        os.unlink(self.__map)


class TestFragmentAlignmentToMap(unittest.TestCase):
    def setUp(self):
        # simulate fragments and chromosomes
        self.__fragment_number = 10
        self.__chromosome_number = 2
        self.__fragment_length = 100
        self.__unplaced_number = 5
        self.__gap_size = 5

        self.__simulator = Simulator(
            self.__fragment_length,
            self.__fragment_number,
            self.__chromosome_number,
            self.__unplaced_number,
            self.__gap_size,
        )

        # create the corresponding files
        self.__map_file = tempfile.mkstemp()[1]
        self.__fragment_file = tempfile.mkstemp()[1]
        self.__chromosome_file = tempfile.mkstemp()[1]

        self.__simulator.write(
            self.__map_file, self.__fragment_file, self.__chromosome_file
        )

        # create the chromosome sequence database and align the
        # fragments to it
        makeblastdb_wrapper = MakeBlastDb(self.__chromosome_file)
        makeblastdb_wrapper.launch()

        self.__alignment_file = tempfile.mkstemp()[1]
        blastn_wrapper = BlastN(
            self.__fragment_file, self.__chromosome_file, self.__alignment_file
        )
        blastn_wrapper.set("-outfmt", 6)
        blastn_wrapper.set("-dust", "no")
        blastn_wrapper.launch()

    def test_blast(self):
        """
        Test the blast method which utilizes BLASTN alignments to
        construct a fragment map.
        """
        fragment_lengths = SeqLengths(self.__fragment_file)
        map_creator = AlignmentToMap(self.__gap_size, fragment_lengths.lengths())
        with open(self.__alignment_file) as alignment_file:
            blast_alignments = BlastTab(alignment_file)
            new_map = map_creator.blast(blast_alignments, 1.2)[0]
            orig_map = Map()
            orig_map.read(self.__map_file)

            # compare the obtained fragment map with the original one
            for chromosome in orig_map.chromosomes():
                for orig, new in zip(
                    orig_map.fragments(chromosome), new_map.fragments(chromosome)
                ):
                    self.assertEqual(orig, new)

            # now test againt the situation when a fragment which length
            # is missing is added to the alignments
        with open(self.__alignment_file) as alignment_file:
            blast_alignments = BlastTab(alignment_file)
            incomplete_lengths = fragment_lengths.lengths()
            del incomplete_lengths[sorted(incomplete_lengths.keys())[0]]
            map_creator = AlignmentToMap(
                self.__gap_size, incomplete_lengths, min_fragment_length=50
            )
            with self.assertRaises(AlignmentToMapError):
                map_creator.blast(blast_alignments, 1.2)

    def test_shrink_gaps(self):
        """
        Test the gap shrinkage function.
        """
        test_map = Map()
        end = 0
        for i in range(10):
            for name in ("fr_{}".format(i), "GAP"):
                fr_length = random.randrange(100)
                new_record = Map.Record(
                    name,
                    fr_length,
                    0,
                    fr_length,
                    random.choice(("+", "-")),
                    "ref_1",
                    end,
                    end + fr_length,
                )
                end += fr_length
                test_map.add_record(new_record)

        gap_size = 50
        test_map.shrink_gaps(gap_size)

        # check if gap sizes are of the specified value
        for i in test_map.chromosomes():
            for j in test_map.fragments(i):
                if j.fr_name == "GAP":
                    self.assertEqual(j.fr_length, gap_size)
                    self.assertEqual(j.fr_end - j.fr_start, gap_size)
                    self.assertEqual(j.ref_end - j.ref_start, gap_size)
            # check that fragments are adjacent to each other
            for j, k in zip(
                list(test_map.fragments(i))[:-1], list(test_map.fragments(i))[1:]
            ):
                self.assertEqual(j.ref_end, k.ref_start)

    def tearDown(self):
        for i in (
            self.__map_file,
            self.__fragment_file,
            self.__chromosome_file,
            self.__alignment_file,
        ):
            if os.path.isfile(i):
                os.unlink(i)
        for i in glob.glob("{}*".format(self.__chromosome_file)):
            os.unlink(i)

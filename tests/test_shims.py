import logging
import os
import tempfile
import unittest

import vcfpy as vcf

from chromosomer.shims import (
    BedReader,
    BedRecord,
    BedWriter,
    Gff3Reader,
    Gff3Record,
    Gff3Writer,
)
from chromosomer.transfer import BedTransfer, Gff3Transfer, VcfTransfer


class TestTransfer(unittest.TestCase):
    def setUp(self):
        self.__test_line = os.path.join("data", "fragment_map", "fragment_map_line.txt")
        self.__test_bed = os.path.join("data", "shim_test", "simple.bed")
        self.__test_vcf = os.path.join("data", "shim_test", "simple.vcf")
        self.__test_blast_tab = os.path.join("data", "shim_test", "simple.blasttab3")
        self.__test_gff3 = os.path.join("data", "shim_test", "simple.gff3")
        self.__output_file = tempfile.NamedTemporaryFile().name
        logging.disable(logging.ERROR)

    def test_transfer_gff3(self):
        transferrer = Gff3Transfer(self.__test_line)
        total_features = 0
        transferred_features = 0
        with open(self.__test_gff3, "r") as gff_file:
            with Gff3Writer(self.__output_file) as writer:
                for feature in Gff3Reader(gff_file).records():
                    total_features += 1
                    tx_feat = transferrer.feature(feature)
                    if tx_feat is not None:
                        transferred_features += 1
                        writer.write(tx_feat)
        with open(self.__output_file, "r") as f:
            assert 1 == sum([1 for _ in f])

    def test_transfer_bed(self):
        total_features = 0
        transferred_features = 0
        transferrer = BedTransfer(self.__test_line)
        with open(self.__test_bed, "r") as bed_file:
            with BedWriter(self.__output_file) as writer:
                for feature in BedReader(bed_file).records():
                    total_features += 1
                    tx_feat = transferrer.feature(feature)
                    if tx_feat is not None:
                        transferred_features += 1
                        writer.write(tx_feat)
        with open(self.__output_file, "r") as f:
            assert 0 == sum([1 for _ in f])

    def test_transfer_vcf(self):
        total_features = 0
        transferred_features = 0
        transferrer = VcfTransfer(self.__test_line)
        with open(self.__output_file, "w") as vcf_out:
            with open(self.__test_vcf) as vcf_in:
                reader = vcf.Reader(vcf_in)
                writer = vcf.Writer(vcf_out, reader.header)
                for variant in reader:
                    total_features += 1
                    tx_variant = transferrer.feature(variant)
                    if tx_variant is not None:
                        transferred_features += 1
                        writer.write_record(tx_variant)
        with open(self.__output_file, "r") as f:
            assert 18 == sum([1 for _ in f])

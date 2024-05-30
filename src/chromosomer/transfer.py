#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


from chromosomer.fragment import Map
from chromosomer.shims import BedRecord, Gff3Record


class Transfer(object):
    """
    Implements transfering routines for abstract data.
    """

    def __init__(self, fragment_map):
        """
        Create a Transfer object.

        :param fragment_map: a fragment map for feature transfer
        :type fragment_map: Map
        """
        self.__fragment_map = Map()
        self.__fragment_map.read(fragment_map)

    def find_fragment(self, fragment):
        """
        Given a fragment name, return its record from the fragment
        map. If the specified fragment is absent in the map, return
        None.

        :param fragment: a fragment name
        :type fragment: str
        :return: a fragment map record corresponding to the specified
            fragment
        :rtype Map.Record
        """
        for chromosome in self.__fragment_map.chromosomes():
            for record in self.__fragment_map.fragments(chromosome):
                if record.fr_name == fragment:
                    return record

        return None

    def coordinate(self, fragment, pos):
        """
        Given a position on a fragment, return the corresponding
        coordinates on the assembled chromosomes according to the
        fragment map specified when the object was created.

        :param fragment: a fragment name
        :param pos: a position on a fragment (zero-based)
        :type fragment: str
        :type pos: int
        :return: a tuple of the chromosome name and a position on it
        :rtype: tuple
        """
        fr_record = self.find_fragment(fragment)
        if fr_record is None:
            # the fragment is absent in the assembly, skip the feature
            return None

        chrom = fr_record.ref_chr
        if fr_record.fr_strand == "+":
            chrom_pos = fr_record.ref_start + pos
        else:
            chrom_pos = fr_record.ref_end - pos

        return chrom, chrom_pos


class BedTransfer(Transfer):
    """
    Implements transferring routines for files in the BED format.
    """

    def feature(self, bed_record: BedRecord):
        """
        Given a record from a BED file, return the transferred one.

        :param bed_record: a line from a BED file
        :type bed_record: BedRecord
        :return: a transferred feature
        :rtype: BedRecord
        """
        fr_record = self.find_fragment(bed_record.seq)
        if fr_record is None:
            # the fragment is absent in the assembly, skip it
            return None

        chrom = fr_record.ref_chr
        start = self.coordinate(bed_record.seq, bed_record.start)[1]
        end = self.coordinate(bed_record.seq, bed_record.end)[1]

        # determine the transferred feature strand
        if bed_record.strand is not None:
            feature_strand = -1 if bed_record.strand == "-" else 1
            fragment_strand = -1 if bed_record.strand == "-" else 1
            if feature_strand * fragment_strand == -1:
                strand = "-"
            else:
                strand = "+"
        else:
            strand = None

        transferred_record = list(bed_record)
        transferred_record[0] = chrom
        transferred_record[1] = min(start, end)
        transferred_record[2] = max(start, end)
        transferred_record[5] = strand

        return BedRecord(*transferred_record)


class Gff3Transfer(Transfer):
    """
    Implements transferring routines for files in the GFF3 format.
    """

    def feature(self, gff3_record: Gff3Record):
        """
        Given a record from a GFF3 file, return the transferred one.

        :param gff3_record: a line from a GFF3 file
        :type gff3_record: Gff3Record
        :return: a transferred feature
        :rtype: BedRecord
        """
        fr_record = self.find_fragment(gff3_record.seqid)
        if fr_record is None:
            # the fragment is absent in the assembly, skip it
            return None

        chrom = fr_record.ref_chr
        start = self.coordinate(gff3_record.seqid, gff3_record.start - 1)[1]
        end = self.coordinate(gff3_record.seqid, gff3_record.end)[1]

        # determine the transferred feature strand
        if gff3_record.strand != ".":
            feature_strand = -1 if gff3_record.strand == "-" else 1
            fragment_strand = -1 if gff3_record.strand == "-" else 1
            if feature_strand * fragment_strand == -1:
                strand = "-"
            else:
                strand = "+"
        else:
            strand = "."

        transferred_record = list(gff3_record)
        transferred_record[0] = chrom
        transferred_record[3] = min(start, end) + 1
        transferred_record[4] = max(start, end)
        transferred_record[6] = strand

        return Gff3Record(*transferred_record)


class VcfTransfer(Transfer):
    """
    Implements transferring routines for files in the VCF format.
    """

    def feature(self, vcf_record):
        """
        Given a record from a VCF file, return the transferred one.

        :param vcf_record: a variant record
        :type vcf_record: vcf.model._Record
        :return: a transferred variant
        :rtype: vcf.model._Record
        """
        fr_record = self.find_fragment(vcf_record.CHROM)
        if fr_record is None:
            # the fragment is absent in the assembly, skip it
            return None

        chrom = fr_record.ref_chr
        pos = self.coordinate(vcf_record.CHROM, vcf_record.POS)[1]

        vcf_record.CHROM = chrom
        vcf_record.POS = pos

        return vcf_record

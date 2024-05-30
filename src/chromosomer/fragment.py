#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015-2016 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import random
from collections import defaultdict, namedtuple
from operator import attrgetter

import pyfaidx

from chromosomer.exception import AlignmentToMapError, MapError
# Bioformats no longer supported
# from bioformats.blast import BlastTab as _BlastTab
from chromosomer.shims import BlastTab
from chromosomer.shims import FastaWriter as Writer
from chromosomer.shims import RandomSequence

logging.basicConfig()
logger = logging.getLogger(__name__)


class Map(object):
    """
    The class implements routines related to creation, reading and
    writing a fragment map that describes how genome fragments are
    situated on reference genome chromosomes.
    """

    numeric_values = (1, 2, 3, 6, 7)

    record_names = (
        "fr_name",
        "fr_length",
        "fr_start",
        "fr_end",
        "fr_strand",
        "ref_chr",
        "ref_start",
        "ref_end",
    )

    Record = namedtuple("Record", record_names)

    def __init__(self):
        """
        Initializes a Map object.
        """
        self.__fragments = defaultdict(list)
        self.__block_adding = False

    def add_record(self, new_record):
        """
        Given a new fragment record, add it to the fragment map.

        :param new_record: a record to be added to the map
        :type new_record: Map.Record
        """
        self.__fragments[new_record.ref_chr].append(new_record)

    def read(self, filename):
        """
        Read a fragment map from the specified file. The file records
        are added to the records in the map.

        :param filename: a name of a file to read a fragment map from
        :type: str
        """
        lineno = 0
        i = 0
        with open(filename) as input_map_file:
            for line in input_map_file:
                lineno += 1
                line_parts = line.split("\t", 8)
                if len(line_parts) < 8:
                    logger.error("line %d: the incorrect number of " "columns", lineno)
                    raise MapError
                for i in self.numeric_values:
                    try:
                        line_parts[i] = int(line_parts[i])
                    except ValueError:
                        logger.error(
                            "line %d: the incorrect numeric " "value %s",
                            lineno,
                            line_parts[i],
                        )
                        raise MapError
                new_record = Map.Record(*line_parts)
                self.add_record(new_record)
                i += 1
        logger.debug(
            "map of %d fragments was successfully read from " "%s", i, filename
        )

    @property
    def records(self):
        return self.__fragments

    def chromosomes(self):
        """
        Return an iterator to the chromosomes the fragment map
        describes.

        :return: an iterator to iterate through the fragment map
            chromosomes
        """
        sorted_chromosomes = sorted(self.__fragments.keys())
        for i in sorted_chromosomes:
            yield i

    def fragments(self, chromosome):
        """
        Return an iterator to fragments of the specified chromosome
        describes by the map. If the chromosome is absent, None is
        returned.

        :param chromosome: a chromosome which fragments are to be
            iterated
        :type chromosome: str
        :return: an iterator to iterate through the chromosome's
            fragments
        """
        if chromosome not in self.__fragments:
            logging.error("%s missing in the fragment map", chromosome)
            raise MapError

        sorted_fragments = sorted(
            self.__fragments[chromosome], key=attrgetter("ref_start")
        )

        for i in sorted_fragments:
            yield i

    def write(self, filename):
        """
        Write the fragment map to the specified file.

        :param filename: a name of a file to write the fragment map to
        :type filename: str
        """
        template = "\t".join(["{}"] * len(self.record_names)) + "\n"
        i = 0
        with open(filename, "w") as output_map_file:
            for chromosome in self.chromosomes():
                for fragment in self.fragments(chromosome):
                    new_line = template.format(*fragment)
                    output_map_file.write(new_line)
                    i += 1
        logger.debug(
            "the map of %d fragments successfully written " "to %s", i, filename
        )

    def assemble(self, fragment_filename, output_filename, save_soft_mask=False):
        """
        Assemble chromosome sequences from fragments.

        :param fragment_filename: a name of a FASTA file of fragment
            sequences
        :param output_filename: a name of the output FASTA file of
            the assembled chromosomes
        :param save_soft_mask: save soft-masking in sequences being
            assembled or not
        :type fragment_filename: str
        :type output_filename: str
        :type save_soft_mask: bool
        """
        logger.debug("assembling chromosomes...")
        logger.debug("FASTA of fragments: %s", fragment_filename)
        logger.debug("FASTA of chromosomes: %s", output_filename)
        logger.debug("saving soft mask: %r", save_soft_mask)
        num_fragments = 0
        num_chromosomes = 0

        fragment_fasta = pyfaidx.Fasta(fragment_filename)
        complement = str.maketrans("ATCGatcgNnXx", "TAGCtagcNnXx")
        with Writer(output_filename) as chromosome_writer:
            for chromosome in self.chromosomes():
                seq = []
                for record in self.fragments(chromosome):
                    if record.fr_name == "GAP":
                        record_seq = "N" * (record.fr_end - record.fr_start)
                    else:
                        if record.fr_name not in fragment_fasta:
                            logger.error(
                                "the fragment %s sequence " "missing", record.fr_name
                            )
                            raise MapError
                        record_seq = fragment_fasta[record.fr_name][
                            record.fr_start : record.fr_end
                        ].seq
                        # convert the sequence to non-unicode
                        record_seq = str(record_seq)
                        if not save_soft_mask:
                            record_seq = record_seq.upper()
                        # if the fragment orientation is reverse, then
                        # the reverse complement of the fragment
                        # sequence is written
                        if record.fr_strand == "-":
                            record_seq = record_seq[::-1].translate(complement)
                    seq.append(record_seq)
                    num_fragments += 1
                chromosome_writer.write(chromosome, "".join(seq))
                num_chromosomes += 1

        logger.debug(
            "%d fragments assembled to %d chromosomes", num_fragments, num_chromosomes
        )

    def shrink_gaps(self, gap_size):
        """
        Shrink gaps inserted into the map to the specified size.

        :param gap_size: a required gap size
        :type gap_size: int
        """
        logger.debug("shrinking gaps to %d bp", gap_size)
        total_shift = 0

        # process each chromosome separately
        for chrom in list(self.__fragments.keys()):
            if len(self.__fragments[chrom]) > 1:
                shifts = []
                # iterate through gaps
                for i in self.__fragments[chrom]:
                    if i.fr_name == "GAP":
                        len_diff = i.fr_length - gap_size
                        shifts.append(len_diff)
                        total_shift += len_diff
                    else:
                        shifts.append(0)
                # calculate absolute shifts for chromosome fragments
                accumulated_shift = 0
                for i in range(len(shifts)):
                    accumulated_shift += shifts[i]
                    shifts[i] = accumulated_shift
                for i in range(len(self.__fragments[chrom])):
                    fragment = list(self.__fragments[chrom][i])
                    if fragment[0] == "GAP":
                        # change gap length and end position
                        fragment[1] = gap_size
                        fragment[3] = gap_size
                        fragment[6] -= shifts[i - 1]
                        fragment[7] = fragment[6] + gap_size
                    else:
                        fragment[6] -= shifts[i]
                        fragment[7] -= shifts[i]
                    self.__fragments[chrom][i] = Map.Record(*fragment)

        logger.debug("in total, gaps shrinked by %d bp", total_shift)

    def summary(self):
        """
        Return a summary on the fragment map.

        :return: a dictionary of tuples each describing one assembled
            chromosome
        :rtype: dict
        """
        summary = {}
        for chromosome in self.chromosomes():
            gaps = 0
            chr_fragments = list(self.fragments(chromosome))
            for fragment in chr_fragments:
                if fragment.fr_name == "GAP":
                    gaps += fragment.fr_length
            chr_length = chr_fragments[-1].ref_end
            fr_num = len(chr_fragments)
            summary[chromosome] = (fr_num, gaps, chr_length)
        return summary

    def convert2bed(self, bed_filename):
        """
        Given a name of the output BED3 file, write the fragment map
        to it in the BED format.

        :param bed_filename: a name of the output BED file of the
            fragment map
        :type bed_filename: str
        """
        template = "\t".join(["{}"] * (len(Map.record_names) + 1)) + "\n"
        with open(bed_filename, "w") as bed_file:
            for chromosome in sorted(self.chromosomes()):
                for fragment in self.fragments(chromosome):
                    bed_file.write(
                        template.format(
                            fragment.ref_chr,
                            fragment.ref_start,
                            fragment.ref_end,
                            fragment.fr_name,
                            1000,
                            fragment.fr_strand,
                            fragment.fr_start,
                            fragment.fr_end,
                            fragment.fr_length,
                        )
                    )


class AlignmentToMap(object):
    """
    The class implements routines to create a fragment map from a set
    of alignments between fragments to be assembled and reference
    chromosomes.
    """

    Anchor = namedtuple(
        "Anchor",
        (
            "fragment",
            "fr_start",
            "fr_end",
            "fr_strand",
            "ref_chr",
            "ref_start",
            "ref_end",
        ),
    )

    def __init__(
        self, gap_size, fragment_lengths, min_fragment_length=None, centromeres=None
    ):
        """
        Create a converter object to create fragment maps from
        alignmets between reference chromosomes and fragments to be
        assembled.

        :param gap_size: a size of a gap between fragments
        :param fragment_lengths: a dictionary of fragment lengths
            which keys are their names and values are their lengths
        :param min_fragment_length: the minimal length of a fragment
            to be included in the map
        :param centromeres: a dictionary of reference chromosome
            centromere locations
        :type gap_size: int
        :type fragment_lengths: dict
        :type min_fragment_length: int
        :type centromeres: dict
        """
        self.__gap_size = gap_size
        self.__fragment_lengths = fragment_lengths
        self.__min_fragment_length = min_fragment_length
        self.__centromeres = centromeres

        self.__anchors = {}
        self.__unlocalized = []
        self.__unplaced = []
        self.__fragment_map = Map()

    def blast(self, blast_alignments: BlastTab, bitscore_ratio_threshold: float):
        """
        Create a fragment map from BLAST blast_alignments between
        fragments and reference chromosomes.

        :param blast_alignments: BLAST blast_alignments
        :param bitscore_ratio_threshold: the minimal ratio of two
            greatest fragment alignment bit scores to consider the
            fragment placed to a reference
        :return: a tuple containing the fragment map constructed from
            the provided BLAST alignments, the list of unlocalized
            fragments and a list of unplaced fragments
        :rtype: tuple
        """
        self.__anchors = {}
        self.__unlocalized = []
        self.__unplaced = []

        temp_anchors = defaultdict(list)

        for alignment in blast_alignments.alignments():
            if self.__min_fragment_length is not None:
                # check if the fragment length is equal or greater
                # than the threshold value
                try:
                    if (
                        self.__fragment_lengths[alignment.query]
                        < self.__min_fragment_length
                    ):
                        # skip the alignment
                        continue
                except KeyError:
                    logger.error("the fragment %s length is missing", alignment.query)
                    raise AlignmentToMapError

            # consider the centromeres if required
            if (
                self.__centromeres is not None
                and alignment.subject in self.__centromeres
            ):
                # the chromosome a fragment was aligned to has a
                # centromere, so we determine which arm the alignment
                # refers to and modify the chromosome name by adding
                # '_1' or '_2' to it
                if (
                    min(alignment.s_start, alignment.s_end)
                    < self.__centromeres[alignment.subject].start
                ):
                    arm_prefix = "_1"
                else:
                    arm_prefix = "_2"
                new_alignment = list(alignment)
                new_alignment[1] += arm_prefix
                alignment = BlastTab.Alignment(*new_alignment)

            temp_anchors[alignment.query].append(alignment)
            # check if there is more than 2 alignments for the
            # fragment; if there is, then leave two fragments with
            # the greatest bit-score values
            if len(temp_anchors[alignment.query]) > 2:
                temp_anchors[alignment.query] = sorted(
                    temp_anchors[alignment.query],
                    key=attrgetter("bit_score"),
                    reverse=True,
                )
                temp_anchors[alignment.query] = temp_anchors[alignment.query][0:2]

        for fragment, alignments in temp_anchors.items():
            if len(alignments) > 1:
                # check if the ratio of the alignment bit scores is
                # greater than the required threshold to consider a
                # fragment places
                if (
                    alignments[0].bit_score / alignments[1].bit_score
                    > bitscore_ratio_threshold
                ):
                    self.__anchors[fragment] = AlignmentToMap.Anchor(
                        fragment=alignments[0].query,
                        fr_start=alignments[0].q_start - 1,
                        fr_end=alignments[0].q_end,
                        fr_strand="+"
                        if alignments[0].s_start < alignments[0].s_end
                        else "-",
                        ref_chr=alignments[0].subject,
                        ref_start=min(alignments[0].s_start, alignments[0].s_end) - 1,
                        ref_end=max(alignments[0].s_start, alignments[0].s_end),
                    )
                elif alignments[0].subject == alignments[1].subject:
                    # the fragment is considered unlocalized
                    self.__unlocalized.append((fragment, alignments[0].subject))
                else:
                    # the fragment is considered unplaced
                    self.__unplaced.append(fragment)
            else:
                # there is a single alignment, use it as an anchor
                self.__anchors[fragment] = AlignmentToMap.Anchor(
                    fragment=alignments[0].query,
                    fr_start=alignments[0].q_start - 1,
                    fr_end=alignments[0].q_end,
                    fr_strand="+"
                    if alignments[0].s_start < alignments[0].s_end
                    else "-",
                    ref_chr=alignments[0].subject,
                    ref_start=min(alignments[0].s_start, alignments[0].s_end) - 1,
                    ref_end=max(alignments[0].s_start, alignments[0].s_end),
                )

        # get total lengths of mapped, unlocalized and unplaced
        # fragments
        total_mapped = total_unlocalized = total_unplaced = 0
        for i in self.__anchors.values():
            total_mapped += self.__fragment_lengths[i.fragment]
        for i in self.__unlocalized:
            total_unlocalized += self.__fragment_lengths[i[0]]
        for i in self.__unplaced:
            total_unplaced += self.__fragment_lengths[i]

        logger.info(
            "%d mapped fragments of total length %d bp",
            len(self.__anchors),
            total_mapped,
        )
        logger.info(
            "%d unlocalized fragments of total length %d bp",
            len(self.__unlocalized),
            total_unlocalized,
        )
        logger.info(
            "%d unplaced fragments of total length %d bp",
            len(self.__unplaced),
            total_unplaced,
        )

        self.__anchor_fragments()

        return (self.__fragment_map, self.__unlocalized, self.__unplaced)

    def __anchor_fragments(self):
        """
        Build a fragment map from anchors.
        """
        # first, we split anchors by reference genome chromosomes
        chr_anchors = defaultdict(list)
        for anchor in self.__anchors.values():
            chr_anchors[anchor.ref_chr].append(anchor)

        # second, we sort the anchors by their position on the
        # chromosomes
        for chr_name in chr_anchors.keys():
            chr_anchors[chr_name] = sorted(
                chr_anchors[chr_name], key=attrgetter("ref_start")
            )

        # now we form a fragment map from the anchors
        total_inserted_gaps = 0
        self.__fragment_map = Map()
        for chr_name in chr_anchors.keys():
            previous_end = 0
            for anchor in chr_anchors[chr_name]:
                try:
                    fragment_length = self.__fragment_lengths[anchor.fragment]
                except KeyError:
                    logger.error("the fragment %s length is missing", anchor.fragment)
                    raise AlignmentToMapError

                # determine the fragment's start and end positions
                if anchor.fr_strand == "+":
                    ref_start = anchor.ref_start - anchor.fr_start
                    ref_end = ref_start + fragment_length
                else:
                    ref_end = anchor.ref_end + anchor.fr_start
                    ref_start = ref_end - fragment_length

                new_record = Map.Record(
                    fr_name=anchor.fragment,
                    fr_length=fragment_length,
                    fr_start=0,
                    fr_end=fragment_length,
                    fr_strand=anchor.fr_strand,
                    ref_chr=anchor.ref_chr,
                    ref_start=previous_end,
                    ref_end=previous_end + ref_end - ref_start,
                )
                self.__fragment_map.add_record(new_record)
                previous_end += ref_end - ref_start

                # add a gap
                new_gap = Map.Record(
                    fr_name="GAP",
                    fr_length=self.__gap_size,
                    fr_start=0,
                    fr_end=self.__gap_size,
                    fr_strand="+",
                    ref_chr=anchor.ref_chr,
                    ref_start=previous_end,
                    ref_end=previous_end + self.__gap_size,
                )
                previous_end += self.__gap_size
                self.__fragment_map.add_record(new_gap)
                total_inserted_gaps += self.__gap_size

        logger.info("%d chromosomes formed", len(chr_anchors))
        logger.info("%d bp of gaps inserted", total_inserted_gaps)


class Simulator(object):
    """
    The class describes routines to simulate genome fragments and
    chromosomes that are composed from them.
    """

    def __init__(
        self,
        fragment_length,
        fragment_number,
        chromosome_number,
        unplaced_number,
        gap_size,
    ):
        """
        Create a fragment simulator object.

        :param fragment_length: the length of a fragment
        :param fragment_number: the number of fragments constituting
            the chromosomes
        :param chromosome_number: the number of chromosomes
        :param unplaced_number: the number of fragments not included
            in the chromosomes
        :param gap_size: the length of gaps between fragments in
            chromosomes
        :type fragment_length: int
        :type fragment_number: int
        :type chromosome_number: int
        :type unplaced_number: int
        :type gap_size: int
        """
        self.__fragment_length = fragment_length
        self.__fragment_number = fragment_number
        self.__chromosome_number = chromosome_number
        self.__gap_size = gap_size

        # create fragment sequences
        self.__fragments = {}
        seq_generator = RandomSequence(self.__fragment_length)
        for i in range(self.__fragment_number):
            fr_name = "fragment{}".format(i + 1)
            self.__fragments[fr_name] = seq_generator.get()

        self.__map = Map()
        self.__create_map()
        self.__assemble_chromosomes()

        # add unplaced fragments
        for i in range(
            self.__fragment_number, self.__fragment_number + unplaced_number
        ):
            fr_name = "fragment{}".format(i + 1)
            self.__fragments[fr_name] = seq_generator.get()

    def __create_map(self):
        """
        Assign fragments to chromosomes randomly and create a
        fragment map.
        """
        fragment_positions = [0] * self.__chromosome_number
        for i in range(self.__fragment_number):
            chr_num = random.randrange(self.__chromosome_number)
            fr_strand = random.choice(("+", "-"))
            self.__map.add_record(
                Map.Record(
                    fr_name="fragment{}".format(i + 1),
                    fr_length=self.__fragment_length,
                    fr_start=0,
                    fr_end=self.__fragment_length,
                    fr_strand=fr_strand,
                    ref_chr="chr{}".format(chr_num + 1),
                    ref_start=fragment_positions[chr_num],
                    ref_end=fragment_positions[chr_num] + self.__fragment_length,
                )
            )
            fragment_positions[chr_num] += self.__fragment_length
            self.__map.add_record(
                Map.Record(
                    fr_name="GAP",
                    fr_length=self.__gap_size,
                    fr_start=0,
                    fr_end=self.__gap_size,
                    fr_strand="+",
                    ref_chr="chr{}".format(chr_num + 1),
                    ref_start=fragment_positions[chr_num],
                    ref_end=fragment_positions[chr_num] + self.__gap_size,
                )
            )
            fragment_positions[chr_num] += self.__gap_size

    def __assemble_chromosomes(self):
        """
        Get chromosome sequences from fragments using the constructed
        fragment map.
        """
        complement = str.maketrans("ATCGatcgNnXx", "TAGCtagcNnXx")
        chromosomes = defaultdict(list)
        for i in self.__map.chromosomes():
            for fr in self.__map.fragments(i):
                if fr.fr_name == "GAP":
                    temp_fragment = "N" * fr.fr_length
                else:
                    temp_fragment = self.__fragments[fr.fr_name]
                    if fr.fr_strand == "-":
                        temp_fragment = temp_fragment[::-1].translate(complement)
                chromosomes[i].append(temp_fragment)
            chromosomes[i] = "".join(chromosomes[i])

        self.__chromosomes = chromosomes

    def write(self, map_file, fragment_file, chromosome_file):
        """
        Write the produced data - a fragment map, a FASTA file of
        fragments and a FASTA file of chromosomes - to the specified
        files.

        :param map_file: a name of a file to write the fragment map to
        :param fragment_file: a name of a file to write fragment
            sequences to
        :param chromosome_file: a name of a file to write chromosome
            sequences to
        :type map_file: str
        :type fragment_file: str
        :type chromosome_file: str
        """
        self.__map.write(map_file)
        with Writer(fragment_file) as fragment_fasta:
            for name, seq in self.__fragments.items():
                fragment_fasta.write(name, seq)
        with Writer(chromosome_file) as chromosome_fasta:
            for name, seq in self.__chromosomes.items():
                chromosome_fasta.write(name, seq)

        logger.debug(
            "a simulated map of %d fragments written to %s",
            len(self.__map.records),
            map_file,
        )
        logger.debug(
            "%d simulated fragments written to %s", len(self.__fragments), fragment_file
        )
        logger.debug(
            "%d simulated chromosomes written to %s",
            len(self.__chromosomes),
            chromosome_file,
        )


class SeqLengths(object):
    """
    The class implements routines to handle fragment sequence lengths.
    """

    def __init__(self, filename):
        """
        Create a SeqLengths object to handle sequence lengths of the
        specified FASTA file.

        :param filename: a name of a FASTA file with sequences which
            lengths are to be derived
        :type filename: str
        """
        self.__filename = filename
        self.__lengths = {}

    def lengths(self):
        """
        Return a dictionary of sequence lengths.

        :return: a dictionary which keys are sequence names and
            values are their lengths
        :rtype: dict
        """
        total_length = 0
        if not self.__lengths:
            reader = pyfaidx.Fasta(self.__filename)
            for seq in list(reader.keys()):
                self.__lengths[seq] = len(reader[seq])
                total_length += len(reader[seq])

        logger.debug(
            "%d sequences analyzed with the total length of " "%d bp",
            len(self.__lengths),
            total_length,
        )

        return self.__lengths


def agp2map(agp_filename, map_filename):
    """
    Given a name of an AGP file, convert it to the fragment map format.

    :param agp_filename: a name of an AGP file
    :param map_filename: a name of an output fragment map file
    :type agp_filename: str
    :type map_filename: str
    """
    with open(agp_filename) as agp_file:
        with open(map_filename, "w") as map_file:
            for line in agp_file:
                if line.startswith("#"):
                    continue
                line = line.rstrip()
                line = line.split(None, 8)
                if line[4] == "N":
                    output = (
                        "GAP",
                        int(line[5]),
                        0,
                        int(line[5]),
                        "+",
                        line[0],
                        int(line[1]) - 1,
                        int(line[2]),
                    )
                else:
                    frag_len = int(line[7]) - int(line[6]) + 1
                    output = (
                        line[5],
                        frag_len,
                        0,
                        frag_len,
                        line[8],
                        line[0],
                        int(line[1]) - 1,
                        int(line[2]),
                    )
                map_file.write("\t".join(map(str, output)) + "\n")

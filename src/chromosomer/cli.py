#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import csv
import logging
import os
from os.path import splitext

import vcfpy as vcf

from chromosomer.fragment import (
    AlignmentToMap,
    Map,
    SeqLengths,
    Simulator,
    agp2map,
    logger,
)
from chromosomer.shims import BedReader, BedWriter, BlastTab, Gff3Reader, Gff3Writer
from chromosomer.transfer import BedTransfer, Gff3Transfer, VcfTransfer


def read_fragment_lengths(filename):
    """
    Given a name of a file with fragment lengths, read them to a
    dictionary.

    :param filename: a name of a file with fragment lengths
    :type filename: str
    :return: a dictionary which keys are fragment sequence names and
        values are their lengths
    :rtype: dict
    """
    result = dict()
    with open(filename) as length_file:
        length_reader = csv.reader(length_file, delimiter="\t")
        for fragment, length in length_reader:
            result[fragment] = int(length)
    return result


def chromosomer():
    """
    The main function that is run if Chromosomer was launched. It
    defines a command-line parser which processed arguments passed to
    the program.
    """
    parser = argparse.ArgumentParser(
        description="Reference-assisted chromosome assembly tool."
    )
    subparsers = parser.add_subparsers(dest="command")

    parser.add_argument("-v", "--version", action="version", version="%(prog)s 0.1.4")

    parser.add_argument(
        "-d", "--debug", action="store_true", help="show debugging messages"
    )

    # Parser for the 'chromosomer assemble' part that produces a FASTA
    # file of assembled chromosomes from the specified fragment map.
    assemble_parser = subparsers.add_parser(
        "assemble",
        help="get sequences of assembled chromosomes",
        description="Get the FASTA file of assembled chromosomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # required arguments for the 'assemble' routine
    assemble_parser.add_argument("map", help="a fragment map file")
    assemble_parser.add_argument(
        "fragment_fasta", help="a FASTA file of fragment " "sequences to be assembled"
    )
    assemble_parser.add_argument(
        "output_fasta",
        help="the output FASTA file of the " "assembled chromosome sequences",
    )

    # optinal arguments for the 'assemble' routine
    assemble_parser.add_argument(
        "-s",
        "--save_soft_mask",
        action="store_true",
        help="keep soft masking from the " "original fragment sequences",
    )

    # Parser for the 'chromosomer fragmentmap' part that
    # produces a map of fragment positions on reference
    # chromosomes from BLAST alignments of the fragments to the
    # chromosomes.
    fragmentmap_parser = subparsers.add_parser(
        "fragmentmap",
        description="Construct a fragment map from fragment "
        "alignments to reference chromosomes.",
        help="construct a fragment map from alignments",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # required arguments for the 'fragmentmap' routine
    fragmentmap_parser.add_argument(
        "alignment_file",
        help="a BLAST tabular file of fragment alignments to " "reference chromosomes",
    )
    fragmentmap_parser.add_argument(
        "gap_size", type=int, help="a size of a gap inserted between mapped fragments"
    )
    fragmentmap_parser.add_argument(
        "fragment_lengths",
        help="a file containing lengths of fragment sequences; it can "
        "be obtained using the 'chromosomer fastalength' tool",
    )
    fragmentmap_parser.add_argument(
        "output_map", help="an output fragment map file name"
    )

    # optional arguments for the 'fragmentmap' routine
    fragmentmap_parser.add_argument(
        "-r",
        "--ratio_threshold",
        type=float,
        default=1.2,
        help="the least ratio of two greatest fragment alignment "
        "scores to determine the fragment placed to a reference "
        "genome",
    )

    fragmentmap_parser.add_argument(
        "-s",
        "--shrink_gaps",
        action="store_true",
        help="shrink large interfragment gaps to the specified size",
    )

    # Parser for the 'chromosomer fragmentmapstat' part that reports
    # statistics on a fragment map
    fragmentmapstat_parser = subparsers.add_parser(
        "fragmentmapstat",
        description="Show statistics on a fragment map.",
        help="show fragment map statistics",
    )

    # required arguments for the 'fragmentmapstat' routine
    fragmentmapstat_parser.add_argument("map", help="a fragment map file")
    fragmentmapstat_parser.add_argument(
        "output", help="an output file of " "fragment map statistics"
    )

    # Parser for the 'chromosomer fragmentmapbed' part that converts
    # a fragement map to the BED format
    fragmentmapbed_parser = subparsers.add_parser(
        "fragmentmapbed",
        description="Convert a fragment map to the BED format.",
        help="convert a fragment map to the BED format",
    )

    # required arguments for the 'fragmentmapbed' routine
    fragmentmapbed_parser.add_argument("map", help="a fragment map file")
    fragmentmapbed_parser.add_argument(
        "output", help="an output BED file " "representing the " "fragment map"
    )

    # Parser for the 'chromosomer transfer' part that transfers
    # genome feature annotation from fragments to their assembly
    transfer_parser = subparsers.add_parser(
        "transfer",
        description="Transfer annotated genomic features from "
        "fragments to their assembly.",
        help="transfer annotated features from fragments to " "chromosomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # required arguments for the 'transfer' routine
    transfer_parser.add_argument("map", help="a fragment map file")
    transfer_parser.add_argument(
        "annotation", help="a file of annotated genome " "features"
    )
    transfer_parser.add_argument(
        "output", help="an output file of the " "transfered annotation"
    )

    # optional arguments for the 'transfer' routine
    transfer_parser.add_argument(
        "-f",
        "--format",
        default="bed",
        choices=["bed", "gff3", "vcf"],
        help="the format of a file of " "annotated features (bed, " "gff3 or vcf)",
    )

    # Parser for the 'chromosomer fastalength' part that calculates
    # lengths of sequences in the given FASTA file.
    fastalength_parser = subparsers.add_parser(
        "fastalength",
        description="Get lengths of sequences in the specified FASTA "
        "file (required to build a fragment map).",
        help="get lengths of sequences from a FASTA file",
    )

    # required arguments for the 'fastalength' routine
    fastalength_parser.add_argument(
        "fasta", help="a FASTA file which sequence " "lengths are to be obtained"
    )
    fastalength_parser.add_argument(
        "output", help="an output file of sequence " "lengths"
    )

    # Parser for the 'chromosomer simulator' routine
    simulator_parser = subparsers.add_parser(
        "simulator",
        description="Simulate fragments and test assembly for " "testing purposes.",
        help="fragment simulator for testing purposes",
    )

    # required arguments for the 'simulator' routine
    simulator_parser.add_argument(
        "fr_num", type=int, help="the number of " "chromosome fragments"
    )
    simulator_parser.add_argument("fr_len", type=int, help="the length of fragments")
    simulator_parser.add_argument("chr_num", type=int, help="the number of chromosomes")
    simulator_parser.add_argument("output_dir", help="the directory for output files")
    simulator_parser.add_argument(
        "-g",
        "--gap_size",
        type=int,
        default=2000,
        help="the size of gaps between " "fragments on a chromosome",
    )
    simulator_parser.add_argument(
        "-p", "--unplaced", type=int, help="the number of unplaced " "fragments"
    )
    simulator_parser.add_argument(
        "--prefix", default="", help="the prefix for output file " "names"
    )

    # Parser for the 'chromosomer agp2map' routine
    agp2map_parser = subparsers.add_parser(
        "agp2map",
        description="Convert an AGP file to the fragment map format.",
        help="convert an AGP file to a fragment map",
    )

    # required arguments for the 'agp2map' routine
    agp2map_parser.add_argument("agp_file", help="an AGP file")
    agp2map_parser.add_argument(
        "output_file", help="the output " "fragment map " "file"
    )

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logger.propagate = False
        formatter = logging.Formatter(
            "%(asctime)-15s - %(message)s", "%Y-%m-%d %H:%M:%S"
        )
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        logging.basicConfig()
        cli_logger = logging.getLogger(__name__)
        cli_logger.propagate = False
        cli_logger.addHandler(ch)
        cli_logger.setLevel(logging.INFO)

    if args.command == "assemble":
        fragment_map = Map()
        fragment_map.read(args.map)
        fragment_map.assemble(
            args.fragment_fasta, args.output_fasta, args.save_soft_mask
        )
    elif args.command == "fragmentmap":
        fragment_lengths = read_fragment_lengths(args.fragment_lengths)
        map_creator = AlignmentToMap(args.gap_size, fragment_lengths)
        with open(args.alignment_file) as alignment_file:
            ## TODO: Replace BlastTab
            alignments = BlastTab(alignment_file)
            fragment_map, unlocalized, unplaced = map_creator.blast(
                alignments, args.ratio_threshold
            )
            if args.shrink_gaps:
                fragment_map.shrink_gaps(args.gap_size)
            fragment_map.write(args.output_map)
            # write unlocalized and unplaced fragments
            with open(
                splitext(args.output_map)[0] + "_unlocalized.txt", "w"
            ) as unlocalized_file:
                for i in unlocalized:
                    unlocalized_file.write("{}\t{}\n".format(*i))
            with open(
                splitext(args.output_map)[0] + "_unplaced.txt", "w"
            ) as unplaced_file:
                for i in unplaced:
                    unplaced_file.write("{}\n".format(i))
    elif args.command == "transfer":
        total_count = transferred_count = 0
        if args.format == "bed":
            transferrer = BedTransfer(args.map)
            with open(args.annotation) as input_file:
                with BedWriter(args.output) as output_file:
                    for feature in BedReader(input_file).records():
                        total_count += 1
                        transferred_feature = transferrer.feature(feature)
                        if transferred_feature is not None:
                            transferred_count += 1
                            output_file.write(transferred_feature)
        elif args.format == "gff3":
            transferrer = Gff3Transfer(args.map)
            with open(args.annotation) as input_file:
                with Gff3Writer(args.output) as output_file:
                    for feature in Gff3Reader(input_file).records():
                        total_count += 1
                        transferred_feature = transferrer.feature(feature)
                        if transferred_feature is not None:
                            transferred_count += 1
                            output_file.write(transferred_feature)
        elif args.format == "vcf":
            transferrer = VcfTransfer(args.map)
            reader = vcf.Reader(open(args.annotation))
            writer = vcf.Writer(open(args.output, "w"), reader.header)
            for variant in reader:
                total_count += 1
                transferred_feature = transferrer.feature(variant)
                if transferred_feature is not None:
                    transferred_count += 1
                    writer.write_record(transferred_feature)
            writer.close()

        logger.info("%d features transferred", transferred_count)
        logger.info("%d features skipped", total_count - transferred_count)
    elif args.command == "fastalength":
        seq_lengths = SeqLengths(args.fasta)
        with open(args.output, "wt") as length_file:
            length_writer = csv.writer(length_file, delimiter="\t")
            for header, length in seq_lengths.lengths().items():
                length_writer.writerow(
                    (
                        header,
                        length,
                    )
                )
    elif args.command == "simulator":
        fr_simulator = Simulator(
            args.fr_len, args.fr_num, args.chr_num, args.unplaced, args.gap_size
        )
        map_file = os.path.join(args.output_dir, args.prefix + "map.txt")
        chr_file = os.path.join(args.output_dir, args.prefix + "chromosomes.fa")
        fr_file = os.path.join(args.output_dir, args.prefix + "fragments.fa")
        fr_simulator.write(map_file, fr_file, chr_file)
    elif args.command == "fragmentmapstat":
        fragment_map = Map()
        fragment_map.read(args.map)
        summary = fragment_map.summary()
        template = "\t".join(["{}"] * 4) + "\n"
        with open(args.output, "w") as output_file:
            for chromosome in sorted(summary.keys()):
                output_file.write(template.format(chromosome, *summary[chromosome]))
    elif args.command == "fragmentmapbed":
        fragment_map = Map()
        fragment_map.read(args.map)
        fragment_map.convert2bed(args.output)
    elif args.command == "agp2map":
        agp2map(args.agp_file, args.output_file)

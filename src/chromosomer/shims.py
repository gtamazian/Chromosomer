"""
Collection of classes to handle biological data (BED, BlastTab, and GFF3 files)
    Taken from biformats.blast.BlastTab, no longer available on pypi
    see original at <https://github.com/yqwu1983/bioformats/blob/dev/bioformats/blast.py>
"""

import csv
import logging
import random
from collections import OrderedDict, namedtuple

from future.utils import iteritems


class BlastTab:
    """
    Taken from biformats.blast.BlastTab, no longer available on pypi
    see original at <https://github.com/yqwu1983/bioformats/blob/dev/bioformats/blast.py>
    """

    blast_field_names = (
        "query",
        "subject",
        "identity",
        "length",
        "mismatches",
        "gap_openings",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "e_value",
        "bit_score",
    )

    Alignment = namedtuple("Alignment", blast_field_names)

    def __init__(self, alignments_file):
        self.__line_parts = []
        self.__reader = csv.reader(alignments_file, delimiter="\t")

    def alignments(self):
        """
        Iterate through alignments in the file the object was created
        from.

        :return: an alignment for the file the object was created from
        :rtype: Blast.Alignment
        """
        for self.__line_parts in self.__reader:
            if not self.__line_parts[0].startswith("#"):
                # if the line starts with '#', then it is a
                # comment and we skip it
                yield self.__parse_blast_line()

    def __parse_blast_line(self):
        """
        Parse the current line from the BLAST tabular file.

        :return: an alignment from the file the object was created from
        :rtype: Blast.Alignment
        """
        line_parts = self.__line_parts

        # check if the line contains the proper number of columns
        if len(line_parts) != 12:
            logging.error(
                "line %d: the incorrect number of columns" % self.__reader.line_num
            )
            raise Exception

        # convert numeric values of identity, e-value and bit score
        # to float numbers
        for i in (2, 10, 11):
            try:
                line_parts[i] = float(line_parts[i])
            except ValueError:
                logging.error(
                    "line %d: the incorrect numerical value %s"
                    % self.__reader.line_num,
                    line_parts[i],
                )
                raise Exception

        # convert numeric values of alignment length, the number of
        # mismatches, the number of gap openings and query and subject
        # coordinates
        for i in range(3, 10):
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logging.error(
                    "line %d: the incorrect integer value %s" % self.__reader.line_num,
                    line_parts[i],
                )
                raise Exception

        return BlastTab.Alignment(*line_parts)


def is_itemrgb(bed_field):
    color_vals = bed_field.split(",", 2)
    return len(color_vals) >= 3 and all(
        [
            x.isnumeric() and "." not in x and 0 <= int(x) and int(x) <= 255
            for x in color_vals
        ]
    )


def is_block_starts(bed_field):
    values = bed_field.split(",")
    if (
        not all([x.isnumeric() and "." not in x and int(x) > 0 for x in values])
        or int(values[0]) != 0
    ):
        return False
    prev_value = 0
    for v in values[1:]:
        if prev_value > int(v):
            return False
        prev_value = int(v)
    return True


bed_field_check = (
    lambda _: True,
    lambda x: x.isnumeric() and "." not in x and int(x) >= 0,  # coordinate
    lambda x: x.isnumeric() and "." not in x and int(x) >= 0,
    lambda _: True,
    lambda x: x.isnumeric()
    and "." not in x
    and 0 <= int(x)
    and int(x) <= 1000,  # Score
    lambda x: x in ("+,-"),  # Strand
    lambda x: x.isnumeric() and "." not in x and int(x) >= 0,
    lambda x: x.isnumeric() and "." not in x and int(x) >= 0,
    is_itemrgb,
    lambda x: x.isnumeric() and "." not in x and int(x) > 0,
    lambda x: all(
        [y.isnumeric() and "." not in y and int(y) > 0 for y in x.split(",")]
    ),
    is_block_starts,
)


class BedReader(object):
    """
    This class implements a parser to read data from a file in the
    BED format.
    """

    def __init__(self, handle):
        """
        Given a handle of a file, create a BED reader object to read
        data from it.

        :param handle: a handle of a BED file
        """
        self.__reader = csv.reader(handle, delimiter="\t")
        self.__line_parts = []
        self.__bed_col = 12  # the number of BED columns
        self.__aux_col = 0  # the number of auxiliary columns

    def records(self, check_order=False):
        """
        Iterate through records in the BED file the object was
        created from.

        :param check_order: check if BED records are sorted; if they
            are not, then raise the exception
        :type check_order: bool
        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        prev_seq = ""
        prev_start = -1
        for self.__line_parts in self.__reader:
            new_record = self.__parse_bed_line()
            if check_order and (
                (prev_seq > new_record.seq) or (prev_start > new_record.start)
            ):
                raise Exception(
                    "line %d: BED record order violated" % self.__reader.line_num
                )
            prev_seq = new_record.seq
            prev_start = new_record.start
            yield new_record

    def __get_bed_format(self):
        """
        Determine the number of BED columns and the number of extra
        columns in a BED file.

        :return: a tuple of two numbers: the number of BED columns
            and the number of extra columns
        :rtype: tuple
        """
        i = 0
        for i, value in enumerate(self.__line_parts):
            if i < 12:
                if not bed_field_check[i](value):
                    i -= 1
                    break
            else:
                return 12, len(self.__line_parts) - 12
        return i + 1, len(self.__line_parts) - (i + 1)

    @property
    def bed_columns(self):
        """
        Get the number of BED columns in a BED file being read.

        :return: the number of BED columns
        :rtype: int
        """
        return self.__bed_col

    @property
    def aux_columns(self):
        """
        Get the number of auxiliary columns in a BED file being read.

        :return: the number of auxiliary BED columns
        :rtype: int
        """
        return self.__aux_col

    def __parse_bed_line(self):
        """
        Parse the current line from the BED file.

        :return: a record from the BED file the object was created from
        :rtype: Record
        """
        bed_col, aux_col = self.__get_bed_format()
        self.__bed_col = min(self.__bed_col, bed_col)
        self.__aux_col = max(self.__aux_col, aux_col)
        if self.__bed_col == 7:
            # thickStart and thickEnd columns must be present together
            self.__bed_col = 6
            self.__aux_col += 1
        elif 10 <= self.bed_columns < 12:
            # blockCount, bloclSizes and blockStarts columns must be
            # present together
            self.__aux_col += self.__bed_col - 9
            self.__bed_col = 9

        if self.__bed_col < 3:
            # The first three columns of a BED file are mandatory.
            raise Exception(
                f"BED line {self.__reader.line_num} has fewer than 3 columns ({self.__bed_col}'{self.__line_parts}')"
            )

        # convert numeric values: start, end, score, thick_start,
        # thick_end, block_num, blocks_sizes and block_starts; the
        # given BED line may contain the lesser number of columns,
        # so we adjust the tuple of numeric value positions
        bed_numeric_fields = (1, 2, 4, 6, 7, 9)
        line_numeric_pos = [x for x in bed_numeric_fields if x < self.__bed_col]
        for i in line_numeric_pos:
            self.__line_parts[i] = int(self.__line_parts[i])

        # form the tuple to be returned as a result
        bed_parts = self.__line_parts[: self.__bed_col]
        bed_parts += [None] * (12 - self.__bed_col)
        aux_parts = self.__line_parts[self.__bed_col :]
        result = BedRecord(*bed_parts, extra=aux_parts)

        return result


class BedWriter(object):
    """
    The class implements writing to a file in the BED format.
    """

    def __init__(self, filename):
        """
        Given a name of a file, create a BED writer object to write
        data to it.

        :param filename: a name of a file to write BED records to
        :type filename: str
        """
        self.__filename = filename

    def __enter__(self):
        self.__output = open(self.__filename, "w")
        return self

    def write(self, bed_record):
        """
        Given a BED record, write it to the file specified when the
        object was creted.

        :param bed_record: a BED record to be written to the file
        :type: Record
        """
        bed_record = [x for x in bed_record if x is not None]
        num_fields = len(bed_record)
        # check if the last column contains any values
        if not bed_record[-1]:
            num_fields -= 1
        else:
            num_fields += len(bed_record[-1]) - 1
            bed_record = bed_record[:-1] + bed_record[-1]
        template = "\t".join(["{}"] * num_fields) + "\n"
        self.__output.write(template.format(*bed_record))

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()


bed_columns = (
    "seq",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "color",
    "block_num",
    "block_sizes",
    "block_starts",
    "extra",
)
BedRecord = namedtuple("BedRecord", bed_columns)

gff3_columns = (
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
)
Gff3Record = namedtuple("Gff3Record", gff3_columns)


class FastaWriter:
    def __init__(self, output_filename: str):
        self.__context_flag = False
        self.__output_filename = output_filename

    def __enter__(self):
        self.__context_flag = True
        self.__output_file = open(self.__output_filename, "w")
        return self

    def write(self, name: str, sequence: str):
        if not self.__context_flag:
            raise Exception("Cannot write FASTA outside of 'with ... as' acontext")
        SPLIT = 40
        self.__output_file.write(">" + name)
        self.__output_file.write("\n")
        self.__output_file.write(
            "\n".join([sequence[i : i + SPLIT] for i in range(0, len(sequence), SPLIT)])
        )
        self.__output_file.write("\n")

    def __exit__(self, exc_type, exc_value, exc_tb):
        self.__output_file.close()
        if exc_type is not None:
            raise exc_value


class Gff3Reader(object):
    """
    This class implements a parser to read data from a file in the
    GFF3 format.
    """

    def __init__(self, handle):
        """
        Given a handle of a file, create a GFF3 reader object to read
        data from it.

        :param handle: a handle of a GFF3 file
        """
        self.__reader = csv.reader(handle, delimiter="\t")
        self.__line_parts = []

    def records(self, check_order=False):
        """
        Iterate through records in the GFF3 file the object was
        created from.

        :param check_order: check if GFF3 records are sorted; if they
            are not, then raise the exception
        :type check_order: bool
        :return: a record from the GFF3 file the object was created
        from
        :rtype: Record
        """
        # check the first line of the file
        prev_seq = ""
        prev_start = -1
        for self.__line_parts in self.__reader:
            new_record = self.__parse_gff3_line()
            if check_order and (
                (prev_seq > new_record.seqid)
                or ((prev_seq == new_record.seqid) and (prev_start > new_record.start))
            ):
                raise Exception(
                    "line %d: GFF3 record order violated" % self.__reader.line_num
                )
            prev_seq = new_record.seqid
            prev_start = new_record.start
            yield new_record

    def __parse_gff3_line(self):
        """
        Parse the current line from the GFF3 file.

        :return: a record from the GFF3 file the object was created
            from
        :rtype: Record
        """
        if len(self.__line_parts) < 8:
            raise Exception(
                "line %d: the incorrect number of columns - %d"
                % (self.__reader.line_num, len(self.__line_parts)),
            )

        if len(self.__line_parts) > 9:
            # in the attributes column, some values may contain tab
            # characters that leads to multiple fields; so we
            # concatenate these fields to a single field
            self.__line_parts[8] = "\t".join(self.__line_parts[8:])
            self.__line_parts = self.__line_parts[:9]

        # convert numeric values: start, end, score (if present) and
        # phase (if present)
        line_int_pos = [3, 4]
        if self.__line_parts[7] != ".":
            line_int_pos += [7]
        for i in line_int_pos:
            try:
                self.__line_parts[i] = int(self.__line_parts[i])
            except ValueError:
                raise Exception(
                    "line %d: the incorrect numeric value %s"
                    % (self.__reader.line_num, self.__line_parts[i]),
                )
        if self.__line_parts[5] != ".":
            # if the score is specified, try to convert it to a float
            # number value
            try:
                self.__line_parts[5] = float(self.__line_parts[5])
            except ValueError:
                raise Exception(
                    "line %d: the incorrect float number value %s"
                    % (self.__reader.line_num, self.__line_parts[5]),
                )

        if len(self.__line_parts) == 9:
            # parse the attributes
            attributes = OrderedDict()
            for x in self.__line_parts[8].strip(";").split(";"):
                try:
                    tag, value = x.split("=", 2)
                except ValueError:
                    raise Exception(
                        "line %d: the incorrect GFF3 attribute %s"
                        % (self.__reader.line_num, x),
                    )
                if not tag or not value:
                    raise Exception(
                        "line %d: the incorrect GFF3 attribute %s"
                        % (self.__reader.line_num, x),
                    )
                attributes[tag] = value
            self.__line_parts[8] = attributes
        else:
            self.__line_parts += [None]

        return Gff3Record(*self.__line_parts)


class Gff3Writer(object):
    """
    The class implements writing to a file in the GFF3 format.
    """

    def __init__(self, filename):
        """
        Given a name of a file, create a GFF3 writer object to write
        data to it.

        :param filename: a name of a file to write GFF3 records to
        :type filename: str
        """
        self.__filename = filename

    def __enter__(self):
        self.__output = open(self.__filename, "w")
        self.__output.write("##gff-version 3\n")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()

    def write(self, gff3_record):
        """
        Given a GFF3 record, write it to the file specified when the
        object was created.

        :param gff3_record: a GFF3 record to be written to the file
        :type gff3_record: Record
        """
        gff3_record = list(gff3_record)
        # prepare the attributes line
        if gff3_record[-1]:
            attributes = ["{}={}".format(x, y) for x, y in iteritems(gff3_record[-1])]
            gff3_record[-1] = ";".join(attributes)
        else:
            gff3_record = gff3_record[:-1]

        template = "\t".join(["{}"] * len(gff3_record)) + "\n"
        self.__output.write(template.format(*gff3_record))


class RandomSequence:
    def __init__(self, length):
        self.__length = length

    def get(self):
        nucs = "ATCG"
        return "".join([random.choice(nucs) for _ in range(self.__length)])

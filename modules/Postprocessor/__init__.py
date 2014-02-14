from collections import defaultdict
from formats.Bed import BedRecord


class Postprocessor():

    def __init__(self, config):
        self.config = config

        self.pirna = set()
        self._parse_pirna_records()

        self.mirna = set()
        self._parse_mirna_records()

        self.read_counts = defaultdict(int)

    def _parse_pirna_records(self):
        for record in open(self.config.pirna_records, "rU").xreadlines():
            _record = record.strip().split(",")

            chromosome = _record[0]
            start = int(_record[1])
            end = int(_record[2])
            strand = _record[3]

            self.pirna.add((chromosome, start, end, strand))

    def _parse_mirna_records(self):
        for record in open(self.config.mirna_records, "rU").xreadlines():
            _record = record.strip().split(",")

            chromosome = _record[0]
            start = int(_record[1])
            end = int(_record[2])
            strand = _record[3]

            self.mirna.add((chromosome, start, end, strand))

    def _count_read_usage(self, alignments):
        for alignment in alignments:
            self.read_counts[BedRecord(alignment).name] += 1

    def _intersects_pirna_record(self, bed_record):
        for chromosome, start, end, strand in self.pirna:
            if (
                bed_record.chrom == chromosome and
                bed_record.chrom_start >= start and
                bed_record.chrom_end <= end and
                bed_record.strand == strand
            ):
                return True

        return False

    def _intersects_mirna_record(self, bed_record):
        for chromosome, start, end, strand in self.mirna:
            if (
                bed_record.chrom == chromosome and
                bed_record.chrom_start >= start-4 and
                bed_record.chrom_end <= end+4 and
                bed_record.strand == strand
            ):
                return True

        return False

    def _is_valid_alignment(self, bed_record):
        return(
            self.read_counts[bed_record.name] == 1 and
            not self._intersects_pirna_record(bed_record) and
            not self._intersects_mirna_record(bed_record)
        )

    def process_alignments(self, alignments):
        print "Post-Processing alignments..."

        self._count_read_usage(alignments)

        return set(
            alignment for alignment in alignments if
            self._is_valid_alignment(BedRecord(alignment))
        )
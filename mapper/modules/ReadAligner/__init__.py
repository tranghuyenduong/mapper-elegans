import re
import subprocess

from collections import defaultdict


class AlignmentRecord():

    def __init__(self, record):
        _record = record.strip().split("\t")

        self.read = _record[0]
        self.read_count = int(re.search("(?<=-)[\d]+", self.read).group(0))
        self.strand = _record[1]
        self.ref = _record[2]
        self.ref_public_name = None
        self.chromosome = None
        self.start_coord = int(_record[3])
        self.seq = _record[4]
        self.end_coord = self.start_coord + len(self.seq)

    def __str__(self):
        return "%s\t%i\t%i\t%s\t%i\t%s" % (
            self.chromosome,
            self.start_coord,
            self.end_coord,
            self.ref_public_name,
            self.read_count,
            self.seq
        )


class RecordFormatter():

    def __init__(self, transcript_bed, transcript_parents, public_names, exon_coords):
        self.transcript_bed = defaultdict(tuple)
        self._parse_transcript_bed(transcript_bed)

        self.transcript_parents = defaultdict(str)
        self._parse_transcript_parents(transcript_parents)

        self.public_names = defaultdict(str)
        self._parse_public_names(public_names)

        self.exon_coords = defaultdict(list)
        self._parse_exon_coords(exon_coords)

    def _parse_transcript_bed(self, transcript_bed):
        for line in open(transcript_bed, "rU").xreadlines():
            chromosome, start, end, transcript_id, _, strand = line.strip().split("\t")

            if strand == "+":
                self.transcript_bed[transcript_id] = chromosome, int(start), 1
            else:
                self.transcript_bed[transcript_id] = chromosome, int(end), -1

    def _parse_transcript_parents(self, transcript_parents):
        for line in open(transcript_parents, "rU").xreadlines():
            transcript_id, name, gene_id = line.strip().split(",")

            self.transcript_parents[transcript_id] = gene_id
            self.transcript_parents[name] = gene_id

    def _parse_public_names(self, public_names):
        for line in open(public_names, "rU").xreadlines():
            gene_id, public_name = line.strip().split(",")

            self.public_names[gene_id] = public_name

    def _parse_exon_coords(self, exon_coords):
        for line in open(exon_coords, "rU").xreadlines():
            transcript_id, coords = line.strip().split("\t")

            self.exon_coords[transcript_id] = map(int, coords.split())

    def format_record(self, alignment_record, sync_coords):
        if sync_coords:
            alignment_record.start_coord = self.exon_coords[alignment_record
            .ref][alignment_record.start_coord]
            alignment_record.end_coord = self.exon_coords[alignment_record
            .ref][alignment_record.end_coord]

        alignment_record.chromosome = self.transcript_bed[alignment_record.ref][0]
        alignment_record.ref_public_name = self.public_names[self
        .transcript_parents[alignment_record.ref]]
        alignment_record.start_coord = self.transcript_bed[alignment_record.ref][1] + alignment_record.start_coord*self.transcript_bed[alignment_record.ref][2]
        alignment_record.end_coord = self.transcript_bed[alignment_record.ref][1] + alignment_record.end_coord*self.transcript_bed[alignment_record.ref][2]


class ReadAligner():

    def __init__(self, config, reads):
        self.config = config
        self.reads = reads

        self.formatter = RecordFormatter(
            self.config.transcript_bed,
            self.config.transcript_parents,
            self.config.public_names,
            self.config.exon_coords
        )

    def align_to(self, ref, sync_coords=False):
        align = subprocess.Popen(
            [
                "bowtie",
                "-f",
                "-v 0",
                "--all",
                "--best",
                "--strata",
                "--nofw",
                ref,
                self.reads
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        while True:
            line = align.stdout.readline()

            if line:
                alignment_record = AlignmentRecord(line)
                self.formatter.format_record(alignment_record, sync_coords)

                yield str(alignment_record)
            else:
                break

        print align.communicate()[1]

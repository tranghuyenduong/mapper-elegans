import subprocess

from collections import defaultdict
from formats.Alignment import AlignmentRecord


class RecordFormatter():

    def __init__(self, transcript_bed, transcript_parents, exon_coords):
        self.complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

        self.transcript_bed = defaultdict(tuple)
        self._parse_transcript_bed(transcript_bed)

        self.transcript_parents = defaultdict(str)
        self._parse_transcript_parents(transcript_parents)

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

    def _parse_exon_coords(self, exon_coords):
        for line in open(exon_coords, "rU").xreadlines():
            transcript_id, coords = line.strip().split("\t")

            self.exon_coords[transcript_id] = map(int, coords.split())

    def format_record(self, alignment_record, coding_transcript):
        if not coding_transcript:
            return

        chromosome, start, strand = self.transcript_bed[alignment_record.ref]

        alignment_record.start_coord = \
            start + self.exon_coords[alignment_record.ref][alignment_record.start_coord]*strand
        alignment_record.end_coord = \
            start + self.exon_coords[alignment_record.ref][alignment_record.end_coord]*strand

        alignment_record.ref = chromosome

        if strand == -1:
            alignment_record.strand = "+"
            alignment_record.start_coord, alignment_record.end_coord = \
                alignment_record.end_coord, alignment_record.start_coord
            alignment_record.seq = "".join(self.complement[base] for base in\
                    reversed(alignment_record.seq))


class ReadAligner():

    def __init__(self, config, reads):
        self.config = config
        self.reads = reads

        self.formatter = RecordFormatter(
            self.config.transcript_bed,
            self.config.transcript_parents,
            self.config.exon_coords
        )

    def align_to(self, ref, coding_transcript=False):
        bt_params = [
            "bowtie",
            "-f",
            "-v 0",
            "--all",
            "--best",
            "--strata",
            ref,
            self.reads
        ]
        if coding_transcript:
            bt_params.insert(6, "--nofw")

        align = subprocess.Popen(
            bt_params,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        std_out, std_err = align.communicate()

        for result in std_out.splitlines():
            alignment_record = AlignmentRecord(result)
            self.formatter.format_record(alignment_record, coding_transcript)

            yield str(alignment_record)

        print std_err
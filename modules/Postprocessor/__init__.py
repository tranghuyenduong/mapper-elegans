import subprocess

from collections import defaultdict
from formats.Bed import BedRecord


class Postprocessor():

    def __init__(self, config):
        self.config = config
        self.read_loci = defaultdict(int)
        self.pirna_mirna_reads = defaultdict(bool)

    def _count_mapped_loci(self, alignments):
        for alignment in alignments:
            self.read_loci[alignment.name] += 1

    def _correct_read_counts(self, alignments):
        for alignment in alignments:
            alignment.score = alignment.score / self.read_loci[alignment.name]

    def _find_pirna_mirna_reads(self, alignments):
        intersect = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-c",
                "-f",
                "1.0",
                "-s",
                "-a",
                "stdin",
                "-b",
                self.config.pirna_mirna_records
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdin = "\n".join(str(i) for i in alignments)
        for result in intersect.communicate(input=stdin)[0].splitlines():
            _, _, _, name, _, _, intersections = result.split("\t")

            if int(intersections):
                self.pirna_mirna_reads[name] = True

    def process_alignments(self, alignments):
        print "Post-Processing alignments..."

        _alignments = [BedRecord(a) for a in alignments]
        self._count_mapped_loci(_alignments)
        self._correct_read_counts(_alignments)
        self._find_pirna_mirna_reads(_alignments)

        return set(
            str(a) for a in _alignments
            if not self.pirna_mirna_reads[a.name]
        )

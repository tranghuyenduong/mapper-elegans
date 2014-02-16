import subprocess

from collections import defaultdict
from formats.Bed import BedRecord


class Postprocessor():

    def __init__(self, config):
        self.config = config
        self.read_scores = defaultdict(int)

    def _calc_read_scores(self, alignments):
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

        for result in intersect.communicate(input="\n".join(alignments))[0].splitlines():
            _, _, _, name, score, _, intersections = result.split("\t")
            self.read_scores[name] += int(score) + int(intersections)

    def process_alignments(self, alignments):
        print "Post-Processing alignments..."

        self._calc_read_scores(alignments)

        return set(
            alignment for alignment in alignments if
            self.read_scores[BedRecord(alignment).name] == 1
        )
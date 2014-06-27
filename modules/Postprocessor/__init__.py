import subprocess

from collections import defaultdict


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
            mapped_loci = self.read_loci[alignment.name]

            if mapped_loci > 1:
                alignment.is_multi_mapped = True

            alignment.score = alignment.score / mapped_loci

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

        self._count_mapped_loci(alignments)
        self._correct_read_counts(alignments)
        self._find_pirna_mirna_reads(alignments)

        return set(a for a in alignments if not self.pirna_mirna_reads[a.name])

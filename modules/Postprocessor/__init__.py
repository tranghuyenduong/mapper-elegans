import subprocess

from collections import defaultdict
from formats.Intersection import IntersectionRecord


class Postprocessor():

    def __init__(self, config):
        self.config = config
        self.read_loci = defaultdict(int)

    def _count_mapped_loci(self, pp_map):
        for a in pp_map:
            self.read_loci[a.name] += 1

    def _correct_read_counts(self, pp_map):
        for a in pp_map:
            a.mapped_loci = self.read_loci[a.name]
            a.score = a.score / a.mapped_loci

    def _find_pirna_mirna_reads(self, pp_map):
        intersect = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-wao",
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

        stdin = "\n".join(str(i) for i in pp_map)
        for result in intersect.communicate(input=stdin)[0].splitlines():
            ir = IntersectionRecord(*result.strip().split())

            if ir.s_name != ".":
                pp_map[(
                    ir.q_chrom,
                    ir.q_chrom_start,
                    ir.q_chrom_end,
                    ir.q_strand
                )].append(ir.s_name)

    def process_alignments(self, alignments):
        print "Post-Processing alignments..."

        pp_map = {a: a.pirnas_mirnas for a in alignments}

        self._count_mapped_loci(pp_map)
        self._correct_read_counts(pp_map)
        self._find_pirna_mirna_reads(pp_map)

        return set(a for a in pp_map)

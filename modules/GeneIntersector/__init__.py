import subprocess

from collections import defaultdict
from formats.Intersection import IntersectionRecord


class GeneIntersector():

    def __init__(self, config):
        self.config = config
        self.gene_intersects = defaultdict(set)

    def _find_gene_intersections(self, alignments):
        intersect = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-wo",
                "-f",
                "1.0",
                "-S",
                "-a",
                "stdin",
                "-b",
                self.config.genes
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        for result in intersect.communicate(input="\n".join(alignments))[0].splitlines():
            intersect_record = IntersectionRecord(result)

            self.gene_intersects[intersect_record.q_name].add(
                "%s\t%i\t%i\t%s\t%i\t%s" % (
                    intersect_record.q_chrom,
                    intersect_record.q_chrom_start,
                    intersect_record.q_chrom_end,
                    intersect_record.s_name,
                    intersect_record.q_score,
                    intersect_record.s_strand
                )
            )

    def find_gene_intersections(self, alignments):
        print "Extracting gene intersections..."

        self._find_gene_intersections(alignments)

        return set(alignment.pop() for alignment in self.gene_intersects.values() if len(alignment) == 1)
import subprocess

from collections import defaultdict
from formats.Intersection import IntersectionRecord


class GeneIntersector():

    def __init__(self, config):
        self.config = config

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
                alignments,
                "-b",
                self.config.genes
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        while True:
            line = intersect.stdout.readline()

            if line:
                yield IntersectionRecord(line)
            else:
                break

        print intersect.communicate()[1]

    def find_gene_intersections(self, alignments):
        read_alignments = defaultdict(set)

        for intersect_record in self._find_gene_intersections(alignments):
            read_alignments[intersect_record.q_name].add(
                "%s\t%i\t%i\t%s\t%i\t%s\t%s" % (
                    intersect_record.q_chrom,
                    intersect_record.q_chrom_start,
                    intersect_record.q_chrom_end,
                    intersect_record.s_name,
                    intersect_record.q_score,
                    intersect_record.q_name,
                    intersect_record.s_strand
                )
            )

        for alignment in read_alignments.values():
            if len(alignment) > 1:
                continue

            yield alignment.pop()
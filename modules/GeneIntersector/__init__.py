import subprocess

from formats.Intersection import IntersectionRecord


class GeneIntersector():

    def __init__(self, config):
        self.config = config

    def _find_gene_intersections(self, gi_map):

        # Setup call to bedtools
        bedtools_call = [
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
        ]
        print("Calling Bedtools with the following")
        print(bedtools_call)

        # Makes call to bedtools
        intersect = subprocess.Popen(
            bedtools_call,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        stdin = "\n".join(str(a) for a in gi_map)
        # Trang's note: need to change this to make python3 happy
        stdin = str.encode(stdin)
        for result in intersect.communicate(input=stdin)[0].splitlines():
            # Trang's note:
            # original
            # ir = IntersectionRecord(*result.strip().split())
            
            # convert back to string to make python 3 happy when defining match dictionary
            ir = IntersectionRecord(*result.strip().decode("utf-8").split())

            gi_map[(
                ir.q_chrom,
                ir.q_chrom_start,
                ir.q_chrom_end,
                ir.q_strand
            )].append(ir.s_name)

    def _correct_read_counts(self, gi_map):
        for a in gi_map:
            if a.genes:
                a.score = a.score / len(a.genes)

    def find_gene_intersections(self, alignments):
        print("Extracting gene intersections...")

        gi_map = {a: a.genes for a in alignments}
        self._find_gene_intersections(gi_map)
        self._correct_read_counts(gi_map)

        return set(a for a in gi_map) if self.config.annotate_only \
                else set(a for a in gi_map if a.genes)

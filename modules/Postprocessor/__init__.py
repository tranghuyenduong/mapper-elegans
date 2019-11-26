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

        bedtools_call = [
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
        ]
        #pirna_mirna_records = p("refs/pirna_mirna") in settings.py
        print "Calling Bedtools with the following"
        print bedtools_call

        intersect = subprocess.Popen(
            bedtools_call,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        records_encountered = 0
        records_saved = 0
        stdin = "\n".join(str(i) for i in pp_map)
        for result in intersect.communicate(input=stdin)[0].splitlines():
            ir = IntersectionRecord(*result.strip().split())

            records_encountered += 1
            if ir.s_name != ".":
                pp_map[(
                    ir.q_chrom,
                    ir.q_chrom_start,
                    ir.q_chrom_end,
                    ir.q_strand
                )].append(ir.s_name)
                records_saved += 1

        print "Saved %d of %d records" % (records_saved, records_encountered)

    def process_alignments(self, alignments):
        print "Post-Processing alignments..."

        pp_map = {a: a.pirnas_mirnas for a in alignments}

        # Normalize the score based on the number of reads
        self._count_mapped_loci(pp_map)
        self._correct_read_counts(pp_map)

        readlocifile = open("postprocess_readloci","w")
        for a in pp_map:
            readlocifile.write("%s -- %d | %d\n" % (a.name, self.read_loci[a.name], a.score))

        self._find_pirna_mirna_reads(pp_map)

        return set(a for a in pp_map)

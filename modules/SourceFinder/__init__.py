import subprocess

from formats.Intersection import IntersectionRecord


class SourceFinder(object):

    def __init__(self, config):
        self.config = config

        self.genome_alignments = None
        self.cds_alignments = None

        # Initialize special cases containers as empty sets
        self.intron_exon_alignments = set()
        self.intron_or_intergenic_alignments = set()
        self.exon_exon_alignments = set()
        self.exon_alignments = set()

    def classify_all_alignments(self, genome_alignments, cds_alignments):
        self.genome_alignments = genome_alignments
        self.cds_alignments = cds_alignments

        self._classify_all_alignments()

    @property
    def all_alignments(self):
        return self.genome_alignments|self.cds_alignments

    @property
    def _intron_or_intergenic_or_intron_exon_alignments(self):
        return self.genome_alignments - cds_alignments

    def _classify_all_alignments(self):
        self._classify_intron_exon_alignments()
        self._classify_intron_or_intergenic_alignments()
        self._classify_exon_exon_alignments()
        self._classify_exon_alignments()

    def _classify_intron_exon_alignments(self):
        iiie = self._intron_or_intergenic_or_intron_exon_alignments
        temp = {a: a.source for a in iiie}

        intersect = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-wo",
                "-S",
                "-a",
                "stdin",
                "-b",
                self.config.exons
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdin = "\n".join(str(a) for a in temp)
        for result in intersect.communicate(input=stdin)[0].splitlines():
            ir = IntersectionRecord(*result.strip().split())

            temp[(
                ir.q_chrom,
                ir.q_chrom_start,
                ir.q_chrom_end,
                ir.q_strand
            )] = 'Intron-Exon'

        self.intron_exon_alignments = set(a for a in iiie
                                          if a.source == 'Intron-Exon')

    def _classify_intron_or_intergenic_alignments(self):
        self.intron_or_intergenic_alignments = \
            self._intron_or_intergenic_or_intron_exon_alignments - \
                self.intron_exon_alignments

        for a in self.intron_or_intergenic_alignments:
            a.source = 'Intron'

    def _classify_exon_exon_alignments(self):
        self.exon_exon_alignments = \
            self.cds_alignments - self.genome_alignments

        for a in self.exon_exon_alignments:
            a.source = 'Exon-Exon'

    def _classify_exon_alignments(self):
        self.exon_alignments = \
            self.cds_alignments - self.exon_exon_alignments

        for a in self.exon_alignments:
            a.source = 'Exon'


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
        print("Classifying alignments...")

        self.genome_alignments = genome_alignments
        self.cds_alignments = cds_alignments

        self._classify_all_alignments()

    @property
    def all_alignments(self):
        return (self.intron_exon_alignments |
                self.intron_or_intergenic_alignments |
                self.exon_exon_alignments |
                self.exon_alignments)

    @property
    def _intron_or_intergenic_or_intron_exon_alignments(self):
        return self.genome_alignments - self.cds_alignments

    def _classify_all_alignments(self):
        self._classify_intron_exon_alignments()
        self._classify_intron_or_intergenic_alignments()
        self._classify_exon_exon_alignments()
        self._classify_exon_alignments()

    def _classify_intron_exon_alignments(self):
        iiie = self._intron_or_intergenic_or_intron_exon_alignments
        
        # Create a dictionary out of the set iiie
        temp = {a: a for a in iiie}

        bedtools_call = [
            "bedtools",
            "intersect",
            "-wo",
            "-S",
            "-a",
            "stdin",
            "-b",
            self.config.exons
        ]
        print("Running bedtools with the following call:")
        print(bedtools_call)

        intersect = subprocess.Popen(
            bedtools_call,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdin = "\n".join(str(a) for a in temp)

        # convert to bytes to make python 3 happy
        stdin = str.encode(stdin)

        for result in intersect.communicate(input=stdin)[0].splitlines():
            # original
            # ir = IntersectionRecord(*result.strip().split())
            
            # convert back to string to make python 3 happy when defining match dictionary
            ir = IntersectionRecord(*result.strip().decode("utf-8").split())

            # Pulls out a matching AlignmentRecord from temp dictionary
            # Remember temp dictionary is created from iiie, which is a 
            # set of alignents that are not coding (i.e. iiie = genome - coding)
            match = temp[(
                ir.q_chrom,
                ir.q_chrom_start,
                ir.q_chrom_end,
                ir.q_strand
            )]

            match.source = 'Intron-Exon'
            self.intron_exon_alignments.add(match)

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


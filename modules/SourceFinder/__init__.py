class SourceFinder(object):

    def __init__(self, genome_alignments, cds_alignments):
        self.genome_alignments = genome_alignments
        self.cds_alignments = cds_alignments

    def all_alignments(self):
        return self.genome_alignments|self.cds_alignments

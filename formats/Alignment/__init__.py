from formats.Bed import BedRecord


class AlignmentRecord(BedRecord):

    def __init__(self, chrom, chrom_start, chrom_end, name, score, strand):
        BedRecord.__init__(self, chrom, chrom_start, chrom_end, name, score,
            strand)

        self.is_multi_mapped = False
        self.genes = []

    @property
    def _unique_identifier(self):
        return (
            self.chrom,
            self.chrom_start,
            self.chrom_end,
            self.strand
        )

    def __hash__(self):
        return hash(self._unique_identifier)

    def __eq__(self, x):
        return self._unique_identifier == x

    def __ne__(self, x):
        return self._unique_identifier != x

    @property
    def is_multi_gene(self):
        return len(self.genes) > 1

    def summary(self):
        for gene in self.genes:
            yield "%s\t%i\t%i\t%s\t%i\t%s\t%s\t%s\n" % (
                self.chrom,
                self.chrom_start,
                self.chrom_end,
                gene,
                self.score,
                self.strand,
                self.is_multi_mapped,
                self.is_multi_gene
            )

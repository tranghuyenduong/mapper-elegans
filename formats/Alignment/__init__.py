from formats.Bed import BedRecord


class AlignmentRecord(BedRecord):

    def __init__(self, chrom, chrom_start, chrom_end, name, score, strand):
        BedRecord.__init__(self, chrom, chrom_start, chrom_end, name, score,
            strand)

        self.mapped_loci = 0
        self.pirnas_mirnas = []
        self.genes = []

        # Populated by the SourceFinder module
        self.source = None

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
    def _is_multi_mapped(self):
        return self.mapped_loci > 1

    @property
    def _is_multi_gene(self):
        return len(self.genes) > 1

    @property
    def _is_pirna_mirna_read(self):
        return len(self.pirnas_mirnas) > 0

    @property
    def summary(self):
        return "%s\t%i\t%i\t%s\t%i\t%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\n" % (
            self.chrom,
            self.chrom_start,
            self.chrom_end,
            self.name,
            self.score,
            self.strand,
            self.source,
            self.mapped_loci,
            self._is_multi_mapped,
            ",".join(self.genes),
            self._is_multi_gene,
            ",".join(self.pirnas_mirnas),
            self._is_pirna_mirna_read)

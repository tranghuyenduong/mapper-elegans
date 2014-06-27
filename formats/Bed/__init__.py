class BedRecord():

    def __init__(self, chrom, chrom_start, chrom_end, name, score, strand):
        self.chrom = chrom
        self.chrom_start = int(chrom_start)
        self.chrom_end = int(chrom_end)
        self.name = name

        try:
            self.score = int(score)
        except ValueError:
            self.score = score

        self.strand = strand

    def __str__(self):
        return "%s\t%i\t%i\t%s\t%i\t%s" % (
            self.chrom,
            self.chrom_start,
            self.chrom_end,
            self.name,
            self.score,
            self.strand
        )


def parse(handle):
    for record in handle.xreadlines():
        yield BedRecord(*record.strip().split("\t"))

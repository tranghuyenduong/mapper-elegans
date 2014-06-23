class BedRecord():

    def __init__(self, record):
        _record = record.strip().split("\t")

        self.chrom = _record[0]
        self.chrom_start = int(_record[1])
        self.chrom_end = int(_record[2])
        self.name = _record[3]

        try:
            self.score = int(_record[4])
        except ValueError:
            self.score = _record[4]

        self.strand = _record[5]

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
        yield BedRecord(record)

class BedRecord():

    def __init__(self, record):
        self.string = record.strip()
        _record = self.string.split("\t")

        self.chrom = _record[0]
        self.chrom_start = int(_record[1])
        self.chrom_end = int(_record[2])
        self.name = _record[3]
        self.score = _record[4]
        self.strand = _record[5]

    def __str__(self):
        return self.string


def parse(handle):
    for record in handle.xreadlines():
        yield BedRecord(record)

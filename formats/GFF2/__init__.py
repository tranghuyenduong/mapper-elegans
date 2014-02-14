from collections import defaultdict


class GFF2Record():

    def __init__(self, record):
        self.string = record.strip()
        _record = self.string.split("\t")

        self.seqid = _record[0]
        self.source = _record[1]
        self.type = _record[2]
        self.start = int(_record[3])
        self.end = int(_record[4])
        self.score = _record[5]
        self.strand = _record[6]
        self.phase = _record[7]

        try:
            _record[8]
        except IndexError:
            self.attr = None
        else:
            self.attr = _record[8]

    def __str__(self):
        return self.string


def parse(handle):
    for record in handle.xreadlines():
        if record.startswith("##"):
            continue

        yield GFF2Record(record)

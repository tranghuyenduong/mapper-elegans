from collections import defaultdict


class GFF3Record():

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
            self.attr = defaultdict(list)
            self._parse_attr(_record[8])

    def __str__(self):
        return self.string

    def _parse_attr(self, attributes):
        try:
            for attr in attributes.split(";"):
                charIndex = attr.index("=")
                tag = attr[0:charIndex]
                values = attr[charIndex+1:]
                #tag, values = attr.split("=")

                for value in values.split(","):
                    self.attr[tag].append(value)
        except ValueError:
            print "This attribute could not be parsed: %s" % attributes


def parse(handle):
    for record in handle.xreadlines():
        if record.startswith("##"):
            continue

        yield GFF3Record(record)

import re

class AlignmentRecord():

    def __init__(self, record):
        _record = record.strip().split("\t")

        self.read = _record[0]
        self.read_count = int(re.search("(?<=-)[\d]+", self.read).group(0))
        self.strand = _record[1]
        self.ref = _record[2]
        self.start_coord = int(_record[3])
        self.seq = _record[4]
        self.end_coord = self.start_coord + len(self.seq)

    def __str__(self):
        return "%s\t%i\t%i\t%s\t%i\t%s" % (
            self.ref,
            self.start_coord,
            self.end_coord,
            self.read,
            self.read_count,
            self.strand
        )

class AlignmentRecord():

    def __init__(self, record):
        _record = record.strip().split("\t")

        self.chromosome = _record[0]
        self.start = int(_record[1])
        self.end = int(_record[2])
        self.geneid = _record[3]
        self.read_count = int(_record[4])
        self.readid = _record[5]

def parse(handle):
    for record in handle.xreadlines():
        yield AlignmentRecord(record)

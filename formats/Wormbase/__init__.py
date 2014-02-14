class WormbaseRecord():

    def __init__(self, record):
        self.string = record.strip()
        _record = self.string.split(",")

        self.taxid = int(_record[0])
        self.geneid = _record[1]
        self._public_name = _record[2]
        self.seqid = _record[3]
        self.status = _record[4]

    def __str__(self):
        return self.string

    @property
    def public_name(self):
        return self._public_name if self._public_name else self.seqid


def parse(handle):
    for record in handle.xreadlines():
        if record.startswith("#") :
            continue

        yield WormbaseRecord(record)

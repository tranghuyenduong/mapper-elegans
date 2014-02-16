class IntersectionRecord():

    def __init__(self, record):
        _record = record.strip().split("\t")

        self.q_chrom = _record[0]
        self.q_chrom_start = int(_record[1])
        self.q_chrom_end = int(_record[2])
        self.q_name = _record[3]

        try:
            self.q_score = int(_record[4])
        except ValueError:
            self.q_score = _record[4]

        self.q_strand = _record[5]

        self.s_chrom = _record[6]
        self.s_chrom_start = int(_record[7])
        self.s_chrom_end = int(_record[8])
        self.s_name = _record[9]

        try:
            self.s_score = int(_record[10])
        except ValueError:
            self.s_score = _record[10]

        self.s_strand = _record[11]

        self.overlap = int(_record[12])
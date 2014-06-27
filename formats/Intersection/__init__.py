class IntersectionRecord():

    def __init__(self, q_chrom, q_chrom_start, q_chrom_end, q_name, q_score,
        q_strand, s_chrom, s_chrom_start, s_chrom_end, s_name, s_score,
        s_strand, overlap):
        self.q_chrom = q_chrom
        self.q_chrom_start = int(q_chrom_start)
        self.q_chrom_end = int(q_chrom_end)
        self.q_name = q_name

        try:
            self.q_score = int(q_score)
        except ValueError:
            self.q_score = q_score

        self.q_strand = q_strand

        self.s_chrom = s_chrom
        self.s_chrom_start = int(s_chrom_start)
        self.s_chrom_end = int(s_chrom_end)
        self.s_name = s_name

        try:
            self.s_score = int(s_score)
        except ValueError:
            self.s_score = s_score

        self.s_strand = s_strand

        self.overlap = int(overlap)

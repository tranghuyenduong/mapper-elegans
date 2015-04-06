import re

from formats.Alignment import AlignmentRecord

class BowtieRecord():

    def __init__(self, read, strand, ref, start_coord, seq, read_qual,
        count_ceil):
        self.read = read
        self.read_count = int(re.search('.*-(\d+)$', self.read).group(1))
        self.strand = strand
        self.ref = ref
        self.start_coord = int(start_coord)
        self.seq = seq
        self.end_coord = self.start_coord + len(self.seq)

    def to_alignment(self):
        return AlignmentRecord(*[
            self.ref,
            self.start_coord,
            self.end_coord,
            self.read,
            self.read_count,
            self.strand
        ])

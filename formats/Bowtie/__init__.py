import re

from formats.Alignment import AlignmentRecord

class BowtieRecord():

    #
    # Trang, Michael - The reason it wasn't running is the following line in ReadAligner.align_to
    #
    #     bt_record = BowtieRecord(*result.split("\t"))
    #
    # expects that the number of columns in the result line matches exactly the number of constructor parameters.
    #
    # Long story short, if you end up using -S format, you'll probably need to change the constructor parameters below to match.
    #
    # For now, I've reverted back Sanjeev's original contructor parameters that matches up with bowtie called WITHOUT the -S flag.
    #
    # def __init__(self, read, strand, ref, start_coord, mapping_qual, cigar, mate_ref, mate_offset, infer_insert_size, seq,
    #     read_qual, xa, md=None, nm=None):
    #     self.read = read
    #     self.read_count = int(re.search('.*-(\d+)$', self.read).group(1))
    #
    #     self.original_strand = strand
    #     strand = int(strand)
    #     if(strand & 0x10 != 0):
    #         self.strand = '-'
    #     else:
    #         self.strand = '+'
    #
    #     self.original_ref = ref
    #     self.ref = ref
    #     self.original_start_coord = start_coord
    #     self.start_coord = int(start_coord)
    #     self.mapping_qual = int(mapping_qual)
    #     self.cigar = cigar
    #     self.mate_ref = mate_ref
    #     self.mate_offset = mate_offset
    #     self.infer_insert_size = infer_insert_size
    #     self.seq = seq
    #     self.read_qual = read_qual
    #     self.end_coord = self.start_coord + len(self.seq)
    #     self.original_seq = seq
    #     self.xa = xa
    #     self.md = md
    #     self.nm = nm

    def __init__(self, read, strand, ref, start_coord, seq, read_qual, count_ceil):
        self.read = read
        self.strand = strand
        self.original_ref = ref
        self.ref = ref
        self.original_start_coord = start_coord
        self.start_coord = int(start_coord)
        self.original_seq = seq
        self.seq = seq
        self.read_qual = read_qual
        self.end_coord = self.start_coord + len(self.seq)

        # Works for both `1-23-2` format and just `23` format
        regexResult = re.search('.*-(\d+)$', self.read)
        if regexResult != None:
            self.read_count = int(regexResult.group(1))
        else :
            self.read_count = int(self.read)

    def to_alignment(self):
        return AlignmentRecord(*[
            self.ref,
            self.start_coord,
            self.end_coord,
            self.read,
            self.read_count,
            self.strand,
            self.original_seq,
            self
        ])
    #
    # Trang, Michael - This is formatted for the -S option. Keeping around in case you need it.
    #
    # @property
    # def summary(self):
    #     return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
    #         self.read,
    #         self.original_strand,
    #         self.original_ref,
    #         self.original_start_coord,
    #         self.mapping_qual,
    #         self.cigar,
    #         self.mate_ref,
    #         self.mate_offset,
    #         self.infer_insert_size,
    #         self.original_seq,
    #         self.read_qual
    #         )

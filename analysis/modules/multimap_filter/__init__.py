from collections import defaultdict
from modules import alignment_record


class MultimapFilter():

    def __init__(self, alignments):
        self.alignments = alignments

        self.reads = defaultdict(int)
        self._count_read_usage()

    def _count_read_usage(self):
        for alignment in alignment_record.parse(open(self.alignments, "rU")):
            self.reads[alignment.readid] += 1

    def get_unique_mapped(self):
        for alignment in alignment_record.parse(open(self.alignments, "rU")):
            if self.reads[alignment.readid] == 1:
                yield alignment

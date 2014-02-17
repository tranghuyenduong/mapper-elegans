from collections import defaultdict
from formats import Bed
from math import ceil


class ReadCounter():

    def __init__(self, genes, bin_size):
        self.bin_size = bin_size
        self.bins = defaultdict(list)
        self.counter = {}
        self.samples = []

        self._create_bins(genes)

    def _get_bins(self, start, end):
        bin_count = int(ceil(float(end-start) / float(self.bin_size)))
        adj_bin_size = int(ceil(float(end-start) / float(bin_count)))

        for idx in xrange(start, end, adj_bin_size):
            yield idx

    def _create_bins(self, genes):
        for record in Bed.parse(open(genes, "rU")):
            for bin in self._get_bins(record.chrom_start, record.chrom_end):
                self.counter[record.name, bin] = defaultdict(int)
                self.bins[record.name].append(bin)

    def _find_bin(self, geneid, coord):
        for bin in reversed(self.bins[geneid]):
            if bin <= coord:
                return geneid, bin

    def init_sample(self, sampleid):
        self.samples.append(sampleid)

    def log_alignment(self, sampleid, alignment):
        self.counter[
            self._find_bin(
                alignment.name,
                alignment.chrom_start if alignment.strand == "+" else alignment.chrom_end
            )
        ][sampleid] += alignment.score

    def write_counter(self, output):
        counter = sorted(self.counter.iteritems(), key=lambda s: (s[0][0], s[0][1]))

        self.samples.sort()

        with open(output, "w") as output_handle:
            output_handle.write("Gene:Bin,%s\n" % ",".join(self.samples))

            for (geneid, bin), count_dict in counter:
                output_handle.write("%s:%i,%s\n" % (
                    geneid,
                    bin,
                    ",".join(str(count_dict[sample]) for sample in self.samples)
                ))

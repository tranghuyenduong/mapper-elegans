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
            yield idx, idx + adj_bin_size - 1

    def _create_bins(self, genes):
        print "Using BED file from [%s]" % genes
        for record in Bed.parse(open(genes, "rU")):
            for start, end in self._get_bins(record.chrom_start,
                                             record.chrom_end):
                self.counter[record.name, start, end] = defaultdict(int)
                self.bins[record.name].append((start, end))

    def _find_bin(self, geneid, coord):
        for start, end in reversed(self.bins[geneid]):
            if start <= coord:
                return geneid, start, end

    def init_sample(self, sampleid):
        self.samples.append(sampleid)

    def log_alignment(self, sampleid, gene, coord, score):
        try:
            self.counter[self._find_bin(gene, coord)][sampleid] += score
        except KeyError:
            print sampleid
            print gene
            print coord
            print score
            print self._find_bin(gene, coord)
            print self.bins[gene]
            quit(1)

    def write_counter(self, output):
        counter = sorted(self.counter.iteritems(),
                         key=lambda s: (s[0][0], s[0][1]))

        self.samples.sort()

        with open(output, "w") as output_handle:
            output_handle.write("Gene:Bin,%s\n" % ",".join(self.samples))

            for (geneid, start, end), count_dict in counter:
                output_handle.write("%s:%i-%i,%s\n" % (
                    geneid,
                    start,
                    end,
                    ",".join(str(count_dict[sample]) \
                             for sample in self.samples)
                ))

    def write_bam(self, output):
        counter = sorted(self.counter.iteritems(),
                         key=lambda s: (s[0][0], s[0][1]))

        self.samples.sort()
        #to be continued!
        with open(output+".bam", "w") as output_bam:
            output_bam.write("Gene:Bin,%s\n" % ",".join(self.samples))

            for (geneid, start, end), count_dict in counter:
                output_bam.write("%s:%i-%i,%s\n" % (
                    geneid,
                    start,
                    end,
                    ",".join(str(count_dict[sample]) \
                             for sample in self.samples)
                ))

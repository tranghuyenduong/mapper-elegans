import os
import subprocess

from Bio import SeqIO
from collections import defaultdict
from modules import is_existing_file


class Preprocessor():

    def __init__(self, config):
        self.config = config

    def processed_reads(self):
        clipped = self._clip_adapter(self.config.raw_input)
        collapsed = self._collapse_reads(clipped)

        print self._read_size_dist(collapsed)
        return collapsed

    def _clip_adapter(self, input):
        print "Clipping adapter sequences..."

        output = "%s_clipped.fastq" % os.path.splitext(input)[0]

        if self.config.force_preprocess or not is_existing_file(output):
            clip_seq = subprocess.Popen(
                [
                    "fastx_clipper",
                    "-a",
                    self.config.adapter_seq,
                    "-c",
                    "-M",
                    str(self.config.min_overlap),
                    "-l",
                    str(self.config.min_seq_len),
                    "-v",
                    "-i",
                    input,
                    "-o",
                    output
                ]
            )
            clip_seq.wait()

        print "Adapters clipped!\n"
        return output

    def _collapse_reads(self, input):
        print "Collapsing reads..."

        output = "%s_collapsed.fa" % os.path.splitext(input)[0]

        if self.config.force_preprocess or not is_existing_file(output):
            total_reads = 0
            too_long_reads = 0
            copied_reads = 0
            unique_reads = 0

            reads = defaultdict(int)

            for read in SeqIO.parse(open(input, "rU"), "fastq"):
                total_reads += 1

                if len(read) > self.config.max_seq_len:
                    too_long_reads += 1
                    continue

                reads[str(read.seq)] += 1

            with open(output, "w") as output_handle:
                for idx, (read, read_count) in enumerate(reads.iteritems()):
                    output_handle.write(">%i-%i\n%s\n" % (
                        idx,
                        read_count,
                        read
                    ))

                    copied_reads += read_count
                    unique_reads += 1

            print "Max. Length: %i\nInput: %i reads.\nOutput: %i reads (%i " \
                  "unique).\ndiscarded %i too-long reads." % (
                self.config.max_seq_len,
                total_reads,
                copied_reads,
                unique_reads,
                too_long_reads
            )

        print "Reads collapsed!\n"
        return output

    @staticmethod
    def _read_size_dist(reads):
        read_sizes = defaultdict(int)

        for record in SeqIO.parse(open(reads, "rU"), "fasta"):
            record_len = len(record.seq)
            read_sizes[record_len] += 1

        s = "Read Size Distribution:\n"

        template = "{0:>5}{1:>10}\n"

        s += template.format("Size", "Count")
        for size, count in read_sizes.iteritems():
            s += template.format(size, count)

        return s

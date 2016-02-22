import re
import subprocess
import time

from Bio import SeqIO
from collections import defaultdict
from modules import is_existing_file, new_or_existing_dir, base_file_name
from os import path


class Preprocessor():

    def __init__(self, config):
        self.config = config
        self.tmp = new_or_existing_dir(self.config.tmp_dir)

    def processed_reads(self, reads, barcode):
        trimmed = self._trim_barcode(reads)
        clipped = self._clip_adapter(trimmed, barcode)
        collapsed = self._collapse_reads(clipped)

        print self._read_size_dist(collapsed)
        return collapsed

    def _trim_barcode(self, input):
        if not self.config.trim_barcode:
            return input

        print "Trimming barcodes from 5\' end..."

        output = path.join(self.tmp, base_file_name(input) + "_trimmed.fq")

        if self.config.force_preprocess or not is_existing_file(output):
            trim_seq = subprocess.Popen(
                [
                    "fastx_trimmer",
                    "-f",
                    str(self.config.barcode_len + 1),
                    "-Q33",
                    "-i",
                    input,
                    "-o",
                    output
                ]
            )
            trim_seq.wait()

        print "Barcodes trimmed...\n"
        return output

    def _clip_adapter(self, input, barcode):
        start_time = time.time()
        print "Clipping adapter sequences..."

        output = path.join(self.tmp, base_file_name(input) + "_clipped.fq")

        if self.config.force_preprocess or not is_existing_file(output):
            clip_seq = subprocess.Popen(
                [
                    "fastx_clipper",
                    "-a",
                    self.config.template % barcode,
                    "-c",
                    "-M",
                    str(self.config.min_overlap),
                    "-l",
                    str(self.config.min_seq_len),
                    "-v",
                    "-Q33",
                    "-i",
                    input,
                    "-o",
                    output
                ]
            )
            clip_seq.wait()

        elapsed_seconds = (time.time() - start_time)
        print "Adapters clipped in {0:.0f} seconds\n".format(elapsed_seconds)
        return output

    def _collapse_reads(self, input):
        print "Collapsing reads..."

        output = path.join(self.tmp, base_file_name(input) + "_collapsed.fq")

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
                    output_handle.write(">%i-%i-%i\n%s\n" % (
                        idx,
                        len(read),
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
            record_count = int(re.search(
                r".*-(\d+)$",
                record.name
                ).group(1))
            read_sizes[record_len] += record_count

        s = "Read Size Distribution:\n"

        template = "{0:>5}{1:>10}\n"

        s += template.format("Size", "Count")
        for size, count in read_sizes.iteritems():
            s += template.format(size, count)

        return s

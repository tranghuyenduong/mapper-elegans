import os
import subprocess

from Bio import SeqIO
from collections import defaultdict
from modules import is_existing_file, expr_cutoff


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

            min_expr_cutoff = expr_cutoff(total_reads, self.config.min_expr)

            with open(output, "w") as output_handle:
                for idx, (read, read_count) in enumerate(reads.iteritems()):
                    if read_count >= min_expr_cutoff:
                        output_handle.write(">%i-%i\n%s\n" % (
                            idx,
                            read_count,
                            read
                        ))

                        copied_reads += read_count
                        unique_reads += 1

            print "Total Reads: %i\nDiscarded %i too long reads.\nCopied " \
                  "Reads: %i (Unique: %i)\nWritten to: %s." % (
                total_reads,
                too_long_reads,
                copied_reads,
                unique_reads,
                output
            )

        print "Reads collapsed!\n"
        return output

    def _filter_multi_mapped(self, input):
        print "Filtering out multi-mapped reads..."

        input_prefix = os.path.splitext(input)[0]

        output = "%s_nomultreads.fa" % input_prefix

        if self.config.force_preprocess or not is_existing_file(output):
            temp = "%s_multmapreads.fa" % input_prefix

            align = subprocess.Popen(
                [
                    "bowtie",
                    "-f",
                    "-v 0",
                    "--all",
                    "--best",
                    "--strata",
                    "-m",
                    str(self.config.multi_map_limit),
                    "--max",
                    temp,
                    self.config.genome_ref,
                    input
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            print align.communicate()[1]

            multi_mapped_reads = defaultdict(bool)

            for read in SeqIO.parse(open(temp, "rU"), "fasta"):
                multi_mapped_reads[read.name] = True

            total_reads = 0
            copied_reads = 0

            with open(output, "w") as output_handle:
                for record in SeqIO.parse(open(input, "rU"), "fasta"):
                    total_reads += 1

                    if not multi_mapped_reads[record.name]:
                        SeqIO.write(record, output_handle, "fasta")

                        copied_reads += 1

            print "%i of %i total reads have been written to %s." % (
                copied_reads,
                total_reads,
                output
            )

        print "Multi-mapped reads filtered!\n"
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

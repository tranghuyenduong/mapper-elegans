import re
import os
import subprocess
import multiprocessing
import time

from multiprocessing import Process
from Bio import SeqIO
from collections import defaultdict
from modules import is_existing_file, new_or_existing_dir, base_file_name
from os import path

# Specific to M1 apple
# PULLING FUNCTION TO TOP LEVEL TO AVOID PICKLING ERROR when running map.py and multiprocessing code starts
# Error we were seeing: AttributeError: Can't pickle local object 'Preprocessor._clip_adapter.<locals>.fastx_clipper'
# 

# Problem: fastx_clipper is not multithreaded and SUPER slow on large fastq files.
# Solution: break up fastq file into smaller chunks, process each chunk on a separate core.
def fastx_clipper(chunk_input, chunk_output, template, min_overlap, min_seq_len):
    clip_seq_cmd = [
            "fastx_clipper",
            "-a",
            #clip off the adapter sequence (or linker sequence)
            template,
            #Sanjeev's set up for sarah's data
            #self.config.template % barcode,
            # discard reads without adaptor (comment out only if the reads are clipped)
            "-c",
            # discard reads have less than min_overlap
            "-M",
            str(min_overlap),
            # discard reads less than min_seq_len
            "-l",
            str(min_seq_len),
            "-v",
            "-Q33",
            "-i",
            chunk_input,
            "-o",
            chunk_output
        ]
    print("Clipping with the following parameters.")
    print("%% %s" % (' '.join(str(p) for p in clip_seq_cmd)))

    clip_seq = subprocess.Popen(clip_seq_cmd)
    clip_seq.wait()

    # Clean up
    os.remove(chunk_input)

class Preprocessor():

    def __init__(self, config):
        self.config = config
        self.tmp = new_or_existing_dir(self.config.tmp_dir)

    def processed_reads(self, reads, barcode):

        trimmed = self._trim_barcode(reads)
        clipped = self._clip_adapter(trimmed, barcode)
        collapsed = self._collapse_reads(clipped)
        return collapsed

    def _trim_barcode(self, input):
        if not self.config.trim_barcode:
            return input

        print("Trimming barcodes from 5\' end...")

        output = path.join(self.tmp, base_file_name(input) + "_trimmed.fq")

        if self.config.force_preprocess or not is_existing_file(output):

            trim_seq_cmd = [
                    "fastx_trimmer",
                    # first base pair to keep, the 5th in this case
                    "-f",
                    str(self.config.barcode_len + 1),
                    "-Q33",
                    "-i",
                    input,
                    "-o",
                    output
                ]

            print("Trimming with the following parameters.")
            print("%% %s" % (' '.join(str(p) for p in trim_seq_cmd)))

            trim_seq = subprocess.Popen(trim_seq_cmd)
            trim_seq.wait()

        print("Barcodes trimmed...\n")
        return output

    def _clip_adapter(self, input, barcode):
        start_time = time.time()

        output = path.join(self.tmp, base_file_name(input) + "_clipped.fq")

        if self.config.force_preprocess or not is_existing_file(output):
            print("===== BARCODE={0} TEMPLATE={1}\n".format(barcode, self.config.template))

            # Figure out chunk size for each core to handle
            LINES_PER_RECORD = 4
            number_of_lines = sum(1 for line in open(input))
            number_of_records = number_of_lines / LINES_PER_RECORD
            number_of_cores = multiprocessing.cpu_count()
            number_of_records_per_core = number_of_records / number_of_cores

            print("Clipping adapter sequences across {0} cores...\n".format(number_of_cores))

            # Prepare for multicore processing
            processes = []
            chunk_output_paths = []

            # Break up fastq into smaller chunks
            chunk_index = 0
            chunk_input_handle = None
            with open(input) as input_handle:
                for index, line in enumerate(input_handle):
                    # Lazy create chunked fastq file
                    if chunk_input_handle is None:
                        chunk_input_path = path.join(self.tmp, base_file_name(input) + "_{0}_chunk.fq".format(chunk_index))
                        chunk_input_handle = open(chunk_input_path, "w")

                    # Write line to chunk file to be processed later
                    chunk_input_handle.write(line)

                    # If this is the last line of the last record for a chunk, start processing the file.
                    record_index = index / LINES_PER_RECORD
                    is_last_line_in_record = (index + 1) % LINES_PER_RECORD == 0
                    is_last_record_in_chunk = (record_index + 1) % number_of_records_per_core == 0
                    is_last_chunk = chunk_index == number_of_cores - 1
                    if is_last_line_in_record and is_last_record_in_chunk and not is_last_chunk:
                        print("Starting fastx_clipper on core {0}".format(chunk_index + 1))

                        # Finish writing the chunk
                        chunk_input_handle.close()

                        # Figure out where to put the output file
                        chunk_output_path = path.join(self.tmp, base_file_name(input) + "_{0}_clipped.fq".format(chunk_index))
                        chunk_output_paths.append(chunk_output_path)

                        # Kick off fastx_clipping
                        process = Process(target=fastx_clipper, args=(chunk_input_path, chunk_output_path, self.config.template, self.config.min_overlap, self.config.min_seq_len))
                        process.start()
                        processes.append(process)

                        # Reset state for the next chunk
                        chunk_input_handle = None
                        chunk_index += 1

            if chunk_input_handle is not None:
                print("Starting fastx_clipper on core {0}".format(chunk_index + 1))

                # Finish writing the chunk
                chunk_input_handle.close()

                # Figure out where to put the output file
                chunk_output_path = path.join(self.tmp, base_file_name(input) + "_{0}_clipped.fq".format(chunk_index))
                chunk_output_paths.append(chunk_output_path)

                # Kick off processing
                print("FASTX CLIPPER: input={0} output={1}".format(chunk_input_path, chunk_output_path))
                process = Process(target=fastx_clipper, args=(chunk_input_path, chunk_output_path, self.config.template, self.config.min_overlap, self.config.min_seq_len))
                process.start()
                processes.append(process)

            print("\nWaiting for {0} processes to finish, this could take a while...\n".format(number_of_cores))

            # Wait for all processes to finish
            for process in processes:
                process.join()

            print("Clipping complete, combining results into a single output file...\n")

            # Build back up a single combined output file
            combined_output_handle = open(output, "w")
            for chunk_output_path in chunk_output_paths:
                with open(chunk_output_path) as chunk_output_handle:
                    for line in chunk_output_handle:
                        combined_output_handle.write(line)

                # Clean up
                os.remove(chunk_output_path)

            # And we're done
            combined_output_handle.close()

            elapsed_seconds = (time.time() - start_time)
            print("Adapters clipped in {0:.0f} seconds\n".format(elapsed_seconds))
        else:
            print("Using previously-clipped adapter sequence at {0}\n".format(output))

        return output

    def _collapse_reads(self, input):
        print("Collapsing reads...")

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

            print("Max. Length: %i\nInput: %i reads.\nOutput: %i reads (%i " \
                  "unique).\ndiscarded %i too-long reads." % (
                self.config.max_seq_len,
                total_reads,
                copied_reads,
                unique_reads,
                too_long_reads
            ))

            #print(self._read_size_dist(collapsed))
            print("self._read_size_dist(collapsed)")
            
        print("Reads collapsed!\n")
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

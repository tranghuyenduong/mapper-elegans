import os

# Alignments input and output
INPUT_SUFFIX = "_alignments.txt"
root_path = "/Volumes/seagate/Project/mapper-elegans"

# BED file filtering
def build_bed_filter_table(filter_table, filename):
    with open(filename, "rU") as input_handle:
        for line in input_handle:
            components = line.split("\t")
            chromosome = components[0]
            start = int(components[1])
            end = int(components[2])

            # Lazy init
            if chromosome not in filter_table:
                filter_table[chromosome] = []

            # Array of (start, end) tuples
            chromosome_tuples = filter_table[chromosome]
            chromosome_tuples.append((start, end))

# Combined filter table for r_rna, sn_rna, sno_rna, and t_rna entries
filter_table = {}

# Build filter table for fast processing
R_RNA_BED_FILE = "/Volumes/seagate/Project/mapper-elegans/data/refs/r_rna.bed"
build_bed_filter_table(filter_table, R_RNA_BED_FILE)
SN_RNA_BED_FILE = "/Volumes/seagate/Project/mapper-elegans/data/refs/sn_rna.bed"
build_bed_filter_table(filter_table, SN_RNA_BED_FILE)
SNO_RNA_BED_FILE = "/Volumes/seagate/Project/mapper-elegans/data/refs/sno_rna.bed"
build_bed_filter_table(filter_table, SNO_RNA_BED_FILE)
T_RNA_BED_FILE = "/Volumes/seagate/Project/mapper-elegans/data/refs/t_rna.bed"
build_bed_filter_table(filter_table, T_RNA_BED_FILE)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# Filter helper
def filtered_alignments(alignments_path):
    lines_to_keep = []
    lines_to_remove = []
    with open(alignments_path, "rU") as input_handle:
        line_count = 0
        increment_count = 0
        total_line_count = file_len(alignments_path)
        for line in input_handle:
            # Debug printing progress "5000 of 33689", "10000 of 33689", etc
            line_count += 1
            if line_count / 10000 != increment_count:
                increment_count += 1
                print "{0} of {1}".format(line_count, total_line_count)

            # Parse
            components = line.split("\t")
            chromosome = components[0]
            start = int(components[1])
            end = int(components[2])

            # Filter
            filter_tuples = filter_table[chromosome]
            for filter_tuple in filter_tuples:
                should_keep = True
                filter_start = filter_tuple[0] <= start <= filter_tuple[1]
                filter_end = filter_tuple[0] <= end <= filter_tuple[1]
                if filter_start or filter_end:
                    # print "filter_start={0} filter_end={1}\nalign_start={2} align_end={3}".format(filter_tuple[0], filter_tuple[1], start, end)
                    should_keep = False
                    break

            # Apply filter
            if should_keep:
                lines_to_keep.append(line)
            else:
                lines_to_remove.append(line)

    return (lines_to_keep, lines_to_remove)

# Find all alignment files
for (dirpath, dirnames, filenames) in os.walk(root_path):
    # Filter out everything but alignment files
    alignemnt_filenames = filter(lambda x:x.endswith(INPUT_SUFFIX), filenames)
    for alignment_filename in alignemnt_filenames:

        # Read the file
        alignment_path = "{0}/{1}".format(dirpath, alignment_filename)
        start_line_count = file_len(alignment_path)

        # Filter
        print "Start filtering {0}...".format(alignment_path)
        lines_to_keep, lines_to_remove = filtered_alignments(alignment_path)
        print "End filtering {0} start_lines={1} end_lines={2}".format(alignment_path, start_line_count, len(lines_to_keep))

        # Write lines we kept
        suffix = alignment_filename[-4:]
        output_filename = alignment_filename[:-4] + "_r_sn_sno_t" + suffix
        output_path = "{0}/{1}".format(dirpath, output_filename)
        with open(output_path, "w") as output_handle:
            for line in lines_to_keep:
                output_handle.write(line)

        # Write lines we filtered out
        output_filename = alignment_filename[:-4] + "_r_sn_sno_t_filter_out" + suffix
        output_path = "{0}/{1}".format(dirpath, output_filename)
        with open(output_path, "w") as output_handle:
            for line in lines_to_remove:
                output_handle.write(line)

        print "Wrote results to {0}".format(output_path)

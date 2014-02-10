def write_alignments(records, output):
    if not output:
        return

    print "Generating a sorted list of alignments for reference...\n"

    total_records = sorted(
        [x.strip().split("\t") for x in records],
        key=lambda s: (s[0], s[3], int(s[1]))
    )

    write_count = 0

    with open(output, "w") as output_handle:
        for record in total_records:
            output_handle.write("%s\n" % "\t".join(record))

            write_count += 1

    return write_count

def is_existing_file(fileName):
    try:
        with open(fileName):
            return True
    except IOError:
        return False

def expr_cutoff(reads_count, min_expr):
    return float(min_expr) * (float(reads_count) / float(1000000))
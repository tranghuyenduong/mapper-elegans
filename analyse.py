import settings
import time

from modules.ReadCounter import ReadCounter
from os import listdir, path


def main():
    config = settings.AnalysisConfig

    print "\nInitializing read counter..."

    rc = ReadCounter(config.genes, config.bin_size)

    for f in listdir(config.alignments_dir):
        # Bail if it's not an alignment file.
        # This guards against processing .DS_STORE files
        ALIGNMENTS_SUFFIX = "_alignments.txt"
        if not f.endswith(ALIGNMENTS_SUFFIX):
            continue

        sampleid = f.strip(ALIGNMENTS_SUFFIX)

        print "\nCounting reads for sample %s..." % sampleid

        rc.init_sample(sampleid)

        alignments = 0
        reads = 0
        for alignment in open(path.join(config.alignments_dir, f), "rU"):
            _, chrom_start, chrom_end, _, score, strand, source, _, \
                is_multi_mapped, gene, is_multi_gene, _, \
                is_pirna_mirna_read = alignment.strip().split("\t")

            if gene == '':
                continue

            if config.exclude_multi_mapped and is_multi_mapped == "True":
                continue

            if config.exclude_pirna_mirna and is_pirna_mirna_read == "True":
                continue

            if config.exclude_multi_gene and is_multi_gene == "True":
                continue

            rc.log_alignment(
                sampleid,
                gene,
                int(chrom_start if strand == "+" else chrom_end),
                int(score)
            )

            alignments += 1
            reads += int(score)

        print "%i alignments representing %i reads logged for sample %s." % (alignments, reads, sampleid)

    print "\nWriting read counts table to %s..." % config.output_file

    rc.write_counter(config.output_file)
    rc.write_bam(config.output_bam)


if __name__ == "__main__":
    start_time = time.time()

    main()

    print "\nScript completed in %i seconds." % (time.time() - start_time)

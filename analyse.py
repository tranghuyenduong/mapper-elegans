import settings
import time

from modules.ReadCounter import ReadCounter
from os import listdir, path


def main():
    config = settings.AnalysisConfig

    print "\nInitializing read counter..."

    rc = ReadCounter(config.genes, config.bin_size)

    for f in listdir(config.alignments_dir):
        sampleid = f.strip("_alignments.txt")

        print "\nCounting reads for sample %s..." % sampleid

        rc.init_sample(sampleid)

        alignments = 0
        reads = 0
        for alignment in open(path.join(config.alignments_dir, f), "rU"):
            _, chrom_start, chrom_end, _, score, strand, source, _, \
                is_multi_mapped, gene, is_multi_gene, _, \
                is_pirna_mirna_read = alignment.strip().split("\t")

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


if __name__ == "__main__":
    start_time = time.time()

    main()

    print "\nScript completed in %i seconds." % (time.time() - start_time)

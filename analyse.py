import settings
import time

from formats import Bed
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

        count = 0
        for alignment in Bed.parse(open(path.join(config.alignments_dir, f), "rU")):
            rc.log_alignment(sampleid, alignment)

            count += 1

        print "%i unique alignments logged for sample %s." % (count, sampleid)

    print "\nWriting read counts table to %s..." % config.output_file

    rc.write_counter(config.output_file)


if __name__ == "__main__":
    start_time = time.time()

    main()

    print "\nScript completed in %i seconds." % (time.time() - start_time)

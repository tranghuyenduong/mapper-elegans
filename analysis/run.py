import settings

from modules.multimap_filter import MultimapFilter
from modules.read_counter import ReadCounter
from os import listdir, path


def main():
    config = settings.Settings

    rc = ReadCounter(config.gene_boundaries, config.bin_size)

    for f in listdir(config.alignments_dir):
        sampleid = f.strip("_alignments.txt")

        rc.init_sample(sampleid)

        for alignment in MultimapFilter(path.join(config.alignments_dir, f)).get_unique_mapped():
            rc.log_alignment(sampleid, alignment)

    rc.write_counter(config.output_file)


if __name__ == "__main__":
    main()

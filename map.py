import os
import re
import settings
import time

from glob import glob
from subprocess import Popen

from modules import write_alignments, new_or_existing_dir
from modules.Preprocessor import Preprocessor
from modules.ReadAligner import ReadAligner
from modules.Postprocessor import Postprocessor
from modules.GeneIntersector import GeneIntersector
from modules.SourceFinder import SourceFinder

def main(reads, barcode, output, pre_process_config, bt_config,
        sf_config, post_process_config, gene_intersect_config):
    print "STEP 1: Pre-processing raw reads\n"
    print "Initializing Pre-processor...\n"
    pre_processor = Preprocessor(pre_process_config)
    processed_reads = pre_processor.processed_reads(reads, barcode)
    del pre_processor
    print "Pre-processing complete!\n"

    print "STEP 2: Generating formatted Bowtie alignment records\n"
    genome = set()
    cds = set()

    read_aligner = ReadAligner(
        bt_config,
        processed_reads
    )
    print "Aligning to genome...\n"
    for alignment in read_aligner.align_to(bt_config.genome_ref):
        genome.add(alignment)

    print "Aligning to coding transcripts...\n"
    for alignment in read_aligner.align_to(bt_config.coding_transcripts_ref,
                                           True):
        cds.add(alignment)
    del read_aligner
    print "Bowtie alignment records generated and formatted!\n"

    print "STEP 3: Classifying alignments...\n"
    source_finder = SourceFinder(sf_config)
    source_finder.classify_all_alignments(genome, cds)
    alignments = source_finder.all_alignments
    del source_finder
    print "Classification complete!\n"

    print "STEP 4: Post-processing alignments...\n"
    post_processor = Postprocessor(post_process_config)
    alignments = post_processor.process_alignments(alignments)
    del post_processor
    print "Post-processing complete!\n"

    print "STEP 5: Extracting gene intersections...\n"
    gene_intersector = GeneIntersector(gene_intersect_config)
    alignments = gene_intersector.find_gene_intersections(alignments)
    del gene_intersector
    print "Gene intersections extracted!\n"

    print "STEP 6: Writing alignments to file...\n"
    alignment_count = write_alignments(alignments, output)
    del alignments
    print "%i alignments written to %s." % (
        alignment_count,
        output
    )

def wrapper():
    mapper_config = settings.MapperConfig
    pre_process_config = settings.PreprocessConfig
    bt_config = settings.BowtieConfig
    sf_config = settings.SourceFinderConfig
    post_process_config = settings.PostprocessConfig
    gene_intersect_config = settings.GeneIntersectConfig

    for f in glob(mapper_config.reads_dir + "/*.fastq"):
        print "\nProcessing reads file %s...\n" % f

        root = re.search(".*/(.*)\.fastq", f).group(1).split("_")
        barcode = root[-1]
        output = "%s/%s_alignments.txt" % (
            new_or_existing_dir(mapper_config.alignments_dir),
            "_".join(root[:3]))

        main(f, barcode, output, pre_process_config, bt_config,
             sf_config, post_process_config, gene_intersect_config)

if __name__ == "__main__":
    start_time = time.time()

    wrapper()

    print "Script completed in %i seconds." % (time.time() - start_time)

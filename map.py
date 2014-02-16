import settings
import time

from modules import write_alignments
from modules.Preprocessor import Preprocessor
from modules.ReadAligner import ReadAligner
from modules.Postprocessor import Postprocessor
from modules.GeneIntersector import GeneIntersector


def wrapper():
    print "STEP 1: Pre-processing raw reads\n"
    print "Initializing Pre-processor...\n"
    pre_processor = Preprocessor(settings.PreprocessConfig)
    processed_reads = pre_processor.processed_reads()
    print "Pre-processing complete!\n"

    print "STEP 2: Generating formatted Bowtie alignment records\n"
    bt_config = settings.BowtieConfig
    alignments = set()
    read_aligner = ReadAligner(
        bt_config,
        processed_reads
    )
    print "Aligning to genome...\n"
    for alignment in read_aligner.align_to(bt_config.genome_ref):
        alignments.add(alignment)

    print "Aligning to coding transcripts...\n"
    for alignment in read_aligner.align_to(bt_config.coding_transcripts_ref, True):
        alignments.add(alignment)
    print "Bowtie alignment records generated and formatted!\n"

    print "STEP 3: Post-processing alignments...\n"
    print "Initializing Post-processor...\n"
    post_processor = Postprocessor(settings.PostprocessConfig)
    alignments = post_processor.process_alignments(alignments)
    print "Post-processing complete!\n"

    print "STEP 4: Extracting gene intersections..."
    output_config = settings.Output
    write_alignments(alignments, output_config.temp)
    gene_intersector = GeneIntersector(settings.GeneIntersectConfig)
    alignments = set(gene_intersector.find_gene_intersections(output_config.temp))
    print "Gene intersections extracted!\n"

    print "STEP 5: Writing alignments to file...\n"
    print "Generating a sorted list of alignments for reference...\n"
    alignment_count = write_alignments(alignments, output_config.alignments)
    print "%i alignments written to %s." % (
        alignment_count,
        output_config.alignments
    )


if __name__ == "__main__":
    start_time = time.time()

    wrapper()

    print "Script completed in %i seconds." % (time.time() - start_time)

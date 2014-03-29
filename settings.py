base_dir = "/home/sba/bioinf/noelle/"

class RefsConfig():

    genome = base_dir + "ws240/c_elegans.PRJNA13758.WS240.genomic.fa"
    gff2_records = base_dir + "ws240/c_elegans.PRJNA13758.WS240.annotations.gff2"
    gff3_records = base_dir + "ws240/c_elegans.PRJNA13758.WS240.annotations.gff3"
    gene_ids = base_dir + "ws240/c_elegans.PRJNA13758.WS240.geneIDs.txt"
    ref = base_dir + "refs"


class PreprocessConfig():

    raw_input = base_dir + "reads/WT_Odor_3.fastq"
    force_preprocess = False

    # Trim 5' barcode
    trim_barcode = False # Set to False when analysing ScriptMiner reads
    barcode_len = 4

    # For ScriptMiner Protocol
    template = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG"
    barcode_seq = "CAGATC"

    # Clip 3' adapter
    adapter_seq = template % barcode_seq
    min_overlap = 15
    min_seq_len = 17
    max_seq_len = 27


class BowtieConfig():

    genome_ref = base_dir + "refs/bt/genome"
    coding_transcripts_ref = base_dir + "refs/bt/coding_transcripts"
    transcript_bed = base_dir + "refs/primary_transcripts.bed"
    transcript_parents = base_dir + "refs/transcript_parents"
    exon_coords = base_dir + "refs/exon_coords"


class PostprocessConfig():

    pirna_mirna_records = base_dir + "refs/pirna_mirna"


class GeneIntersectConfig():

    genes = base_dir + "refs/genes"


class Output():

    alignments = base_dir + "alignments/WT_Odor_3_alignments.txt"


class AnalysisConfig():

    genes = base_dir + "refs/genes"
    bin_size = 400
    alignments_dir = base_dir + "alignments"
    output_file = base_dir + "data_400.csv"
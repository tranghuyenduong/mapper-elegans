class RefsConfig():

    genome = "/home/sba/noelle/ws240/c_elegans.PRJNA13758.WS240.genomic.fa"
    gff2_records = "/home/sba/noelle/ws240/c_elegans.PRJNA13758.WS240.annotations.gff2"
    gff3_records = "/home/sba/noelle/ws240/c_elegans.PRJNA13758.WS240.annotations.gff3"
    gene_ids = "/home/sba/noelle/ws240/c_elegans.PRJNA13758.WS240.geneIDs.txt"
    ref = "/home/sba/noelle/refs"


class PreprocessConfig():

    raw_input = "/home/sba/noelle/reads/WT_Odor_3.fastq"
    force_preprocess = False

    scriptminer = True

    # Trim 5' barcode (For non-scriptminer methods)
    barcode_len = 4

    # Clip 3' adapter
    barcode_seq = "CAGATC" # For scriptminer only
    scriptminer_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG"
    adapter_seq = scriptminer_adapter % barcode_seq
    min_overlap = 15
    min_seq_len = 17
    max_seq_len = 27


class BowtieConfig():

    genome_ref = "/home/sba/noelle/refs/bt/genome"
    coding_transcripts_ref = "/home/sba/noelle/refs/bt/coding_transcripts"
    transcript_bed = "/home/sba/noelle/refs/primary_transcripts.bed"
    transcript_parents = "/home/sba/noelle/refs/transcript_parents"
    exon_coords = "/home/sba/noelle/refs/exon_coords"


class PostprocessConfig():

    pirna_mirna_records = "/home/sba/noelle/refs/pirna_mirna"


class GeneIntersectConfig():

    genes = "/home/sba/noelle/refs/genes"


class Output():

    alignments = "/home/sba/noelle/alignments/WT_Odor_3_alignments.txt"


class AnalysisConfig():

    genes = "/home/sba/noelle/refs/genes"
    bin_size = 400
    alignments_dir = "/home/sba/noelle/alignments"
    output_file = "/home/sba/noelle/data_400.csv"

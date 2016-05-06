from os import path

BASE_DIR = "data"

def P(pathname):
    return path.join(BASE_DIR, pathname)

class RefsConfig():

    genome = P("ws253/c_elegans.PRJNA13758.WS253.genomic.fa")
    gff2_records = P("ws253/c_elegans.PRJNA13758.WS253.annotations.gff2")
    gff3_records = P("ws253/c_elegans.PRJNA13758.WS253.annotations.gff3")
    gene_ids = P("ws253/annotation/c_elegans.PRJNA13758.WS253.geneIDs.txt")
    ref = P("refs")

class MapperConfig():

    reads_dir = P("reads")
    alignments_dir = P("alignments")

class PreprocessConfig():

    tmp_dir = P("tmp")

    force_preprocess = False

    # Trim 5' barcode
    trim_barcode = False # Set to False when analysing ScriptMiner reads
    barcode_len = 4

    # For ScriptMiner Protocol
    template = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG"

    # Clip 3' adapter
    min_overlap = 15
    min_seq_len = 17
    max_seq_len = 27


class BowtieConfig():

    genome_ref = P("refs/bt/genome")
    coding_transcripts_ref = P("refs/bt/coding_transcripts")
    transcript_bed = P("refs/primary_transcripts.bed")
    transcript_parents = P("refs/transcript_parents")
    exon_coords = P("refs/exon_coords")


class SourceFinderConfig():

    exons = P("refs/exons.bed")


class PostprocessConfig():

    pirna_mirna_records = P("refs/pirna_mirna")


class GeneIntersectConfig():

    annotate_only = True
    genes = P("refs/genes")


class AnalysisConfig():

    genes = P("refs/genes")
    bin_size = 400
    exclude_multi_mapped = True
    exclude_pirna_mirna = True
    exclude_multi_gene = True
    alignments_dir = P("alignments")
    output_file = P("data_400.csv")
    output_bam = P("data_400.bam")

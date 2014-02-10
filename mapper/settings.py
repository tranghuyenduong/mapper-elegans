class PreprocessConfig():

    raw_input = "/home/sba/noelle/reads/Mut_Odor_1.fastq"
    force_preprocess = False

    # Clip 3' adapter
    barcode_seq = "ACTTGA"
    adapter_seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG" % barcode_seq
    min_overlap = 10
    min_seq_len = 15

    # Filter out low-count reads
    min_expr = 0 # Count per 1 million filtered reads

    # Filter out multi-mapped reads
    multi_map_limit = 10
    genome_ref = "/home/sba/noelle/refs/bt/genome"


class BowtieConfig():

    primary_transcripts_ref = "/home/sba/noelle/refs/bt/primary_transcripts"
    coding_transcripts_ref = "/home/sba/noelle/refs/bt/coding_transcripts"
    transcript_bed = "/home/sba/noelle/refs/primary_transcripts.bed"
    transcript_parents = "/home/sba/noelle/refs/transcript_parents"
    public_names = "/home/sba/noelle/refs/public_names"
    exon_coords = "/home/sba/noelle/refs/exon_coords"


class Output():

    alignments = "/home/sba/noelle/alignments/Mut_Odor_1_alignments.txt"

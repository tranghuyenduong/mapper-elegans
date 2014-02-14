class PreprocessConfig():

    raw_input = "/home/sba/noelle/reads/WT_Odor_3.fastq"
    force_preprocess = False

    # Clip 3' adapter
    barcode_seq = "CAGATC"
    adapter_seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC%sATCTCGTATGCCGTCTTCTGCTTG" % barcode_seq
    min_overlap = 15
    min_seq_len = 17
    max_seq_len = 27


class BowtieConfig():

    genome_ref = "/home/sba/noelle/refs/bt/genome"
    coding_transcripts_ref = "/home/sba/noelle/refs/bt/coding_transcripts"
    transcript_bed = "/home/sba/noelle/refs/primary_transcripts.bed"
    transcript_parents = "/home/sba/noelle/refs/transcript_parents"
    exon_coords = "/home/sba/noelle/refs/exon_coords"


class Output():

    alignments = "/home/sba/noelle/alignments/WT_Odor_3_alignments.txt"

import modules
import os
import settings


def wrapper():
    config = settings.RefsConfig

    print "Generating reference files for genome..."
    print "Please be patient. This may take awhile."

    # Create refs folder
    if not os.path.exists(config.ref):
        os.makedirs(config.ref)

    # Extract protein coding genes
    print "\nExtracting protein coding genes..."
    modules.protein_coding_genes(config.gff3_records, config.ref + "/genes")

    # Extract primary transcripts
    print "\nExtracting primary transcripts..."
    modules.extract_bed(config.gff2_records, config.ref + "/primary_transcripts.bed", "protein_coding_primary_transcript")
    modules.extract_sequences(config.genome, config.ref + "/primary_transcripts.bed", config.ref + "/primary_transcripts.fa")

    # Extract introns
    print "\nExtracting introns..."
    modules.extract_bed(config.gff2_records, config.ref + "/introns.bed", "intron")
    modules.order_records(config.ref + "/introns.bed", config.ref + "/introns_sorted.bed")
    modules.extract_sequences(config.genome, config.ref + "/introns_sorted.bed", config.ref + "/introns_sorted.fa")

    # Extract transcript parents
    print "\nExtracting transcript parents..."
    modules.transcript_parents(config.gff3_records, config.ref + "/transcript_parents")

    # Extract piRNA/miRNA records
    print "\nExtracting piRNA and miRNA records..."
    modules.pirna_mirna_records(config.gff3_records, config.ref + "/pirna_mirna")

    # Extract public names for transcripts
    print "\nExtracting gene public names..."
    modules.public_names(config.gene_ids, config.ref + "/public_names")

    # Part primary transcripts
    print "\nParting primary transcripts..."
    modules.part_primary_transcripts(config.ref + "/introns_sorted.fa", config.ref + "/primary_transcripts.fa", config.ref + "/parted_primary_transcripts")

    # Generate coding transcripts
    print "\nGenerating coding transcripts..."
    modules.coding_transcripts(config.ref + "/introns_sorted.fa", config.ref + "/parted_primary_transcripts", config.ref + "/coding_transcripts.fa")

    # Extract exon coordinates relative to primary transcript
    print "\nExtracting exon coordinates from primary transcript..."
    modules.exon_coords(config.ref + "/introns_sorted.fa", config.ref + "/parted_primary_transcripts", config.ref + "/exon_coords")

    # Generate bowtie indices for genome, primary transcripts, and coding transcripts
    print "\nBuilding the necessary bowtie indices..."

    # Create folder to hold bowtie indices
    bt_indices = config.ref + "/bt"
    if not os.path.exists(bt_indices):
        os.makedirs(bt_indices)
    modules.build_bowtie_index(config.genome, bt_indices + "/genome")
    modules.build_bowtie_index(config.ref + "/coding_transcripts.fa", bt_indices + "/coding_transcripts")
    print "\nReference files generated in %s." % config.ref


if __name__ == "__main__":
    wrapper()

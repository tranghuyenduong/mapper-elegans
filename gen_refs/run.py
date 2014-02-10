import modules
import os
import settings


def wrapper():
    print "Generating reference files for genome..."
    print "Please be patient. This may take awhile."

    # Create refs folder
    if not os.path.exists(settings.REF):
        os.makedirs(settings.REF)

    # Extract primary transcripts
    print "\nExtracting primary transcripts..."

    modules.extract_bed(settings.GFF2_RECORDS, settings.REF + "/primary_transcripts.bed", "protein_coding_primary_transcript")

    modules.extract_sequences(settings.GENOME, settings.REF + "/primary_transcripts.bed", settings.REF + "/primary_transcripts.fa")

    # Extract introns
    print "\nExtracting introns..."

    modules.extract_bed(settings.GFF2_RECORDS, settings.REF + "/introns.bed", "intron")

    modules.order_records(settings.REF + "/introns.bed", settings.REF + "/introns_sorted.bed")

    modules.extract_sequences(settings.GENOME, settings.REF + "/introns_sorted.bed", settings.REF + "/introns_sorted.fa")

    # Extract transcript parents
    print "\nExtracting transcript parents..."

    modules.transcript_parents(settings.GFF3_RECORDS, settings.REF + "/transcript_parents")

    # Extract public names for transcripts
    print "\nExtracting gene public names..."

    modules.public_names(settings.GENE_IDS, settings.REF + "/public_names")

    # Part primary transcripts
    print "\nParting primary transcripts..."

    modules.part_primary_transcripts(settings.REF + "/introns_sorted.fa", settings.REF + "/primary_transcripts.fa", settings.REF + "/parted_primary_transcripts")

    # Generate coding transcripts
    print "\nGenerating coding transcripts..."

    modules.coding_transcripts(settings.REF + "/introns_sorted.fa", settings.REF + "/parted_primary_transcripts", settings.REF + "/coding_transcripts.fa")

    # Extract exon coordinates relative to primary transcript
    print "\nExtracting exon coordinates from primary transcript..."

    modules.exon_coords(settings.REF + "/introns_sorted.fa", settings.REF + "/parted_primary_transcripts", settings.REF + "/exon_coords")

    # Generate bowtie indices for genome, primary transcripts, and coding transcripts
    print "\nBuilding the necessary bowtie indices..."

    # Create folder to hold bowtie indices
    bt_indices = settings.REF + "/bt"

    if not os.path.exists(bt_indices):
        os.makedirs(bt_indices)

    modules.build_bowtie_index(settings.GENOME, bt_indices + "/genome")
    modules.build_bowtie_index(settings.REF + "/primary_transcripts.fa", bt_indices + "/primary_transcripts")
    modules.build_bowtie_index(settings.REF + "/coding_transcripts.fa", bt_indices + "/coding_transcripts")

    print "\nReference files generated in %s." % settings.REF


if __name__ == "__main__":
    wrapper()

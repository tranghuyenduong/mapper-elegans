import modules
import os
import settings
import time

def wrapper():
    config = settings.RefsConfig

    print "Generating reference files for genome..."
    print "Please be patient. This may take awhile."

    # Create refs folder
    if not os.path.exists(config.ref):
        os.makedirs(config.ref)

    def P(path, base_dir = config.ref):
        return os.path.join(base_dir, path)

    # Extract gff2 data
    print "\nExtracting data from gff2 annotation file..."
    modules.extract_gff2_data(config.gff2_records)

    # Extract data from gff3 annotation file
    print "\nExtracting data from gff3 annotation file..."
    modules.extract_gff3_data(config.gff3_records)

    # Extract primary transcripts
    print "\nExtracting primary transcripts..."
    modules.extract_sequences(config.genome, P("primary_transcripts.bed"), P("primary_transcripts.fa"))

    # Extract introns
    print "\nExtracting introns..."
    modules.order_records(P("introns.bed"), P("introns_sorted.bed"))
    modules.extract_sequences(config.genome, P("introns_sorted.bed"), P("introns_sorted.fa"))

    # Extract public names for transcripts
    print "\nExtracting gene public names..."
    modules.public_names(config.gene_ids, P("public_names"))

    # Extract sequence ids for transcripts
    print "\nExtracting gene sequence ids..."
    modules.sequence_ids(config.gene_ids, P("sequence_ids"))

    # Part primary transcripts
    print "\nParting primary transcripts..."
    modules.part_primary_transcripts(P("introns_sorted.fa"), P("primary_transcripts.fa"), P("parted_primary_transcripts"))

    # Generate coding transcripts
    print "\nGenerating coding transcripts..."
    modules.coding_transcripts(P("introns_sorted.fa"), P("parted_primary_transcripts"), P("coding_transcripts.fa"))

    # Extract exon coordinates relative to primary transcript
    print "\nExtracting exon coordinates from primary transcript..."
    modules.exon_coords(P("introns_sorted.fa"), P("parted_primary_transcripts"), P("exon_coords"))

    # Generate bowtie indices for genome, primary transcripts, and coding transcripts
    print "\nBuilding the necessary bowtie indices..."

    # Create folder to hold bowtie indices
    bt_indices = config.ref + "/bt"
    if not os.path.exists(bt_indices):
        os.makedirs(bt_indices)
    modules.build_bowtie_index(config.genome, P("genome", bt_indices))
    modules.build_bowtie_index(P("coding_transcripts.fa"), P("coding_transcripts", bt_indices))
    print "\nReference files generated in %s." % config.ref


if __name__ == "__main__":
    start_time = time.time()

    wrapper()

    print "Script completed in %i seconds." % (time.time() - start_time)

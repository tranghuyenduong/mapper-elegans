import os
import re
import subprocess
import settings
import sys

from Bio import SeqIO
from collections import defaultdict
from formats import Bed, GFF2, GFF3, Wormbase

def P(path):
    return os.path.join(settings.RefsConfig.ref, path)

def extract_gff2_data(gff2_filename):

    # Inner method for writing records to disk that match the given source and type
    def handle_record(gff2_record, record_source, record_type, output_handle):
        if gff2_record.source != record_source or gff2_record.type != record_type:
            return

        ident = re.search(
            r'[Transcript|Gene] \"(?P<ident>.+?)\"',
            gff2_record.attr
        ).group('ident')

        output_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
            re.search("CHROMOSOME_(.*)", gff2_record.seqid).group(1),
            gff2_record.start-1,
            gff2_record.end,
            ident,
            gff2_record.strand
        ))

    # Extract transposons
    transposons_handle = open(P("transposons.bed"), "w")
    def extract_transposons(gff2_record):
        handle_record(gff2_record, "transposon_gene", "gene", transposons_handle)

    # Extract snRNA
    snRNA_handle = open(P("sn_rna.bed"), "w")
    def extract_snRNA(gff2_record):
        handle_record(gff2_record, "snRNA", "snRNA", snRNA_handle)

    # Extract snoRNA
    snoRNA_handle = open(P("sno_rna.bed"), "w")
    def extract_snoRNA(gff2_record):
        handle_record(gff2_record, "snoRNA", "snoRNA", snoRNA_handle)

    # Extract ncRNA
    ncRNA_handle = open(P("nc_rna.bed"), "w")
    def extract_ncRNA(gff2_record):
        handle_record(gff2_record, "ncRNA", "ncRNA", ncRNA_handle)

    # Extract tRNA transcripts
    tRNA_transcripts_handle = open(P("t_rna.bed"), "w")
    def extract_tRNA_transcripts(gff2_record):
        handle_record(gff2_record, "tRNA", "tRNA", tRNA_transcripts_handle)

    # Extract rRNA transcripts
    rRNA_transcripts_handle = open(P("r_rna.bed"), "w")
    def extract_rRNA_transcripts(gff2_record):
        handle_record(gff2_record, "rRNA", "rRNA", rRNA_transcripts_handle)

    # Extract primary transcripts
    primary_transcripts_handle = open(P("primary_transcripts.bed"), "w")
    def extract_primary_transcripts(gff2_record):
        handle_record(gff2_record, "Coding_transcript", "protein_coding_primary_transcript", primary_transcripts_handle)

    # Extract introns
    introns_handle = open(P("introns.bed"), "w")
    def extract_introns(gff2_record):
        handle_record(gff2_record, "Coding_transcript", "intron", introns_handle)

    # Extract coding exons
    coding_exons_handle = open(P("exons.bed"), "w")
    def extract_coding_exons(gff2_record):
        handle_record(gff2_record, "Coding_transcript", "exon", coding_exons_handle)

    # Parse the file, handle each record
    for gff2_record in GFF2.parse(open(gff2_filename, "rU")):
        extract_transposons(gff2_record)
        extract_snRNA(gff2_record)
        extract_snoRNA(gff2_record)
        extract_ncRNA(gff2_record)
        extract_tRNA_transcripts(gff2_record)
        extract_rRNA_transcripts(gff2_record)
        extract_primary_transcripts(gff2_record)
        extract_introns(gff2_record)
        extract_coding_exons(gff2_record)

def extract_gff3_data(gff3_filename):
    record_handlers = ()

    # Extract protein coding genes
    protein_coding_genes_handle = open(P("genes"), "w")
    def extract_protein_coding_gene(gff3_record):
        if gff3_record.source != "WormBase" or gff3_record.type != "gene" or gff3_record.attr["biotype"].pop() != "protein_coding":
            return

        name = gff3_record.attr["Name"].pop() if gff3_record.attr["Name"] else ""

        protein_coding_genes_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
            gff3_record.seqid,
            gff3_record.start-1,
            gff3_record.end,
            name,
            gff3_record.strand
        ))

    # Extract transcript parents
    transcript_parents_handle = open(P("transcript_parents"), "w")
    def extract_transcript_parents(gff3_record):
        if gff3_record.source != "WormBase" or gff3_record.type != "mRNA":
            return

        transcript = [val.split(":")[1] for val in gff3_record.attr["ID"] if val.startswith("Transcript:")].pop()
        geneid = [val.split(":")[1] for val in gff3_record.attr["Parent"] if val.startswith("Gene:")].pop()
        name = gff3_record.attr["Name"].pop() if gff3_record.attr["Name"] else ""

        transcript_parents_handle.write("%s,%s,%s\n" % (
            transcript,
            name,
            geneid
        ))

    # Extract piRNA and miRNA
    pirna_mirna_records_handle = open(P("pirna_mirna"), "w")
    def extract_pirna_mirna_records(gff3_record):
        if gff3_record.source != "WormBase":
            return

        biotype = gff3_record.attr["biotype"].pop() if gff3_record.attr["biotype"] else ""
        name = gff3_record.attr["Name"].pop() if gff3_record.attr["Name"] else ""

        if gff3_record.type == "gene" and biotype == "piRNA":
            pirna_mirna_records_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                gff3_record.seqid,
                gff3_record.start-1,
                gff3_record.end,
                name,
                gff3_record.strand
            ))

        elif gff3_record.type == "miRNA":
            pirna_mirna_records_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                gff3_record.seqid,
                gff3_record.start-5,
                gff3_record.end+4,
                name,
                gff3_record.strand
            ))

    # Parse the file, handle each record
    for gff3_record in GFF3.parse(open(gff3_filename, "rU")):
        extract_protein_coding_gene(gff3_record)
        extract_transcript_parents(gff3_record)
        extract_pirna_mirna_records(gff3_record)

def order_records(bed_file, output):
    records = defaultdict(list)
    for bed_record in Bed.parse(open(bed_file, "rU")):
        records[bed_record.strand].append(bed_record)

    records["+"].sort(key=lambda p: (p.name, p.chrom_start))
    records["-"].sort(key=lambda p: (p.name, -p.chrom_start))

    with open(output, "w") as output_handle:
        for _, bed_records in records.iteritems():
            output_handle.writelines(("%s\n" % str(bed_record)) for bed_record in bed_records)

def extract_sequences(genome, bed_file, output):
    extract = subprocess.Popen(
        [
            "bedtools",
            "getfasta",
            "-name",
            "-s",
            "-fi",
            genome,
            "-bed",
            bed_file,
            "-fo",
            output
        ]
    )
    extract.wait()

def public_names(wormbase_records, output):
    with open(output, "w") as output_handle:
        for wormbase_record in Wormbase.parse(open(wormbase_records, "rU")):
            if wormbase_record.taxid != 6239:
                continue

            output_handle.write("%s,%s\n" % (
                wormbase_record.geneid,
                wormbase_record.public_name,
            ))

def sequence_ids(wormbase_records, output):
    with open(output, "w") as output_handle:
        for wormbase_record in Wormbase.parse(open(wormbase_records, "rU")):
            if wormbase_record.taxid != 6239:
                continue

            output_handle.write("%s,%s\n" % (
                wormbase_record.geneid,
                wormbase_record.seqid,
            ))

def part_primary_transcripts(sorted_introns, primary_transcripts, output):
    introns = defaultdict(list)
    for intron_record in SeqIO.parse(open(sorted_introns, "rU"), "fasta"):
        introns[intron_record.name].append(str(intron_record.seq))

    with open(output, "w") as output_handle:
        for primary_transcript in SeqIO.parse(open(primary_transcripts, "rU"), "fasta"):
            if not introns[primary_transcript.name]:
                continue

            parts = [str(primary_transcript.seq)]

            for intron in introns[primary_transcript.name]:
                suffix = parts.pop(-1)

                for new_part in suffix.partition(intron):
                    parts.append(new_part)

            output_handle.write(">%s\n%s\n" % (
                primary_transcript.name,
                ",".join(parts)
            ))

def coding_transcripts(sorted_introns, parted_pri_transcripts, output):
    introns = defaultdict(list)
    for intron in SeqIO.parse(open(sorted_introns, "rU"), "fasta"):
        introns[intron.name].append(str(intron.seq))

    with open(output, "w") as output_handle:
        for pri_transcript in SeqIO.parse(open(parted_pri_transcripts, "rU"), "fasta"):
            pri_transcript_seq = str(pri_transcript.seq).strip().split(",")

            exons = []

            for block in pri_transcript_seq:
                if not introns[pri_transcript.name]:
                    exons.append(block)

                elif block not in introns[pri_transcript.name]:
                    exons.append(block)

            output_handle.write(">%s\n%s\n" % (
                pri_transcript.name,
                "".join(exons)
            ))

def exon_coords(sorted_introns, parted_pri_transcripts, output):
    introns = defaultdict(list)
    for intron in SeqIO.parse(open(sorted_introns, "rU"), "fasta"):
        introns[intron.name].append(str(intron.seq))

    with open(output, "w") as output_handle:
        for pri_transcript in SeqIO.parse(open(parted_pri_transcripts, "rU"), "fasta"):
            pri_transcript_seq = str(pri_transcript.seq).strip().split(",")

            ranges = []
            start_index = sys.maxsize
            pri_transcript_coord = 0

            # Always capture range for first coordinate
            ranges.append("0:1")

            # Capture ranges for all other sequences
            for block in pri_transcript_seq:
                for base in block:
                    pri_transcript_coord += 1

                    if not introns[pri_transcript.name] or block not in introns[pri_transcript.name]:
                        # Start a new range
                        if start_index == sys.maxsize:
                            start_index = pri_transcript_coord
                    else:
                        # Finish an existing range
                        if start_index != sys.maxsize:
                            length = pri_transcript_coord - start_index
                            ranges.append("%d:%d" % (start_index, length))
                            start_index = sys.maxsize

            # Capture range that includes the last coordinate.
            if start_index != sys.maxsize:
                # Bump up coord by 1 so the last pri_transcript_coord is included in the range
                length = pri_transcript_coord - start_index + 1
                ranges.append("%d:%d" % (start_index, length))

            output_handle.write("%s\t%s\n" % (
                pri_transcript.name,
                " ".join(map(str, ranges))
            ))

def build_bowtie_index(fasta, bt_index):
    build = subprocess.Popen(
        [
            "bowtie-build",
            fasta,
            bt_index
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    build.communicate()

def write_alignments(records, output):
    if not output:
        return

    all_alignments = sorted(records, key=lambda s: (s.chrom, s.chrom_start))

    write_count = 0
    with open(output, "w") as output_handle:
        for alignment in all_alignments:
            output_handle.write(alignment.summary)
            write_count += 1

    return write_count

def is_existing_file(file_name):
    try:
        with open(file_name):
            return True
    except IOError:
        return False

def new_or_existing_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

    return path

def base_file_name(path):
    return os.path.splitext(os.path.basename(path))[0]

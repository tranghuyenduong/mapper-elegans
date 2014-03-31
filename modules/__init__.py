import re
import subprocess

from Bio import SeqIO
from collections import defaultdict
from formats import Bed, GFF2, GFF3, Wormbase


def protein_coding_genes(gff3, output):
    with open(output, "w") as output_handle:
        for gff3_record in GFF3.parse(open(gff3, "rU")):
            if gff3_record.source != "WormBase" or gff3_record.type != \
                    "gene" or gff3_record.attr["biotype"].pop() != \
                    "protein_coding":
                continue

            output_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                gff3_record.seqid,
                gff3_record.start-1,
                gff3_record.end,
                gff3_record.attr["Name"].pop(),
                gff3_record.strand
            ))

def extract_bed(annotations, output, record_type):
    with open(output, "w") as output_handle:
        for record in GFF2.parse(open(annotations, "rU")):
            if record.source != "Coding_transcript":
                continue

            if record.type != record_type:
                continue

            transcript_id = re.search('(?<=Transcript ")[^"]*(?=")', record.attr).group(0)

            output_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                re.search("CHROMOSOME_(.*)", record.seqid).group(1),
                record.start-1,
                record.end,
                transcript_id,
                record.strand
            ))

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

def transcript_parents(gff3, output):
    with open(output, "w") as output_handle:
        for gff3_record in GFF3.parse(open(gff3, "rU")):
            if gff3_record.source != "WormBase" or gff3_record.type != "mRNA":
                continue

            transcript = [val.split(":")[1] for val in gff3_record.attr["ID"] if val.startswith("Transcript:")].pop()
            geneid = [val.split(":")[1] for val in gff3_record.attr["Parent"] if val.startswith("Gene:")].pop()
            name = gff3_record.attr["Name"].pop()

            output_handle.write("%s,%s,%s\n" % (
                transcript,
                name,
                geneid
            ))

def pirna_mirna_records(gff3, output):
    with open(output, "w") as output_handle:
        for gff3_record in GFF3.parse(open(gff3, "rU")):
            if gff3_record.source != "WormBase":
                continue

            if gff3_record.type == "gene" and gff3_record.attr["biotype"].pop() == "piRNA":
                output_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                    gff3_record.seqid,
                    gff3_record.start-1,
                    gff3_record.end,
                    gff3_record.attr["Name"].pop(),
                    gff3_record.strand
                ))

            elif gff3_record.type == "miRNA":
                output_handle.write("%s\t%i\t%i\t%s\t.\t%s\n" % (
                    gff3_record.seqid,
                    gff3_record.start-5,
                    gff3_record.end+4,
                    gff3_record.attr["Name"].pop(),
                    gff3_record.strand
                ))

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

            coords = []
            pri_transcript_coord = 0

            coords.append(pri_transcript_coord)

            for block in pri_transcript_seq:
                for base in block:
                    pri_transcript_coord += 1

                    if not introns[pri_transcript.name]:
                        coords.append(pri_transcript_coord)

                    elif block not in introns[pri_transcript.name]:
                        coords.append(pri_transcript_coord)

            output_handle.write("%s\t%s\n" % (
                pri_transcript.name,
                " ".join(map(str, coords))
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

    total_records = sorted(
        [x.strip().split("\t") for x in records],
        key=lambda s: (s[0], int(s[1]))
    )

    write_count = 0

    with open(output, "w") as output_handle:
        for record in total_records:
            output_handle.write("%s\n" % "\t".join(record))

            write_count += 1

    return write_count

def is_existing_file(fileName):
    try:
        with open(fileName):
            return True
    except IOError:
        return False
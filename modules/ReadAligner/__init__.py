import subprocess

from collections import defaultdict
from formats.Bowtie import BowtieRecord


class RecordFormatter():

    def __init__(self, transcript_bed, transcript_parents, exon_coords):
        self.complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

        self.transcript_bed = defaultdict(tuple)
        self._parse_transcript_bed(transcript_bed)

        self.transcript_parents = defaultdict(str)
        self._parse_transcript_parents(transcript_parents)

        self.exon_coords = defaultdict(list)
        self._parse_exon_coords(exon_coords)

    def _parse_transcript_bed(self, transcript_bed):
        print("Reading transcript BED from %s" % transcript_bed)
        # python 2
        # for line in open(transcript_bed, "rU").xreadlines():
        for line in open(transcript_bed, "rU"):
            chromosome, start, end, transcript_id, _, strand = line.strip().split("\t")

            if strand == "+":
                self.transcript_bed[transcript_id] = chromosome, int(start), 1
            else:
                self.transcript_bed[transcript_id] = chromosome, int(end), -1

    def _parse_transcript_parents(self, transcript_parents):
        print("Reading transcript parents from %s" % transcript_parents)
        # python 2
        # for line in open(transcript_parents, "rU").xreadlines():
        for line in open(transcript_parents, "rU"):
            transcript_id, name, gene_id = line.strip().split(",")

            self.transcript_parents[transcript_id] = gene_id
            self.transcript_parents[name] = gene_id

    def _parse_exon_coords(self, exon_coords):
        print("Reading exon coordinates from %s" % exon_coords)
        # python 2
        # for line in open(exon_coords, "rU").xreadlines():
        for line in open(exon_coords, "rU"):
            transcript_id, range_list_string = line.strip().split("\t")
            
            coords = []
            range_list = range_list_string.split(" ")
            for range_string in range_list:
                range_array = range_string.split(":")
                start_index = int(range_array[0])
                length = int(range_array[1])
                for pri_coord in range(start_index, start_index + length):
                    coords.append(pri_coord)
            self.exon_coords[transcript_id] = coords

    def format_record(self, alignment_record, coding_transcript, result=None):
        if not coding_transcript:
            return

        try:
            # Transforming from gene-based coordinates to chromosome-based coordinates
            chromosome, start, strand = self.transcript_bed[alignment_record.ref]
            alignment_record.start_coord = \
                start + self.exon_coords[alignment_record.ref][alignment_record.start_coord]*strand
            alignment_record.end_coord = \
                start + self.exon_coords[alignment_record.ref][alignment_record.end_coord]*strand

            # Replace the locus reference (e.g. ZK512.10) with chromosome
            alignment_record.ref = chromosome

            # Reverse the coordinates if strand is not +
            if strand == -1:
                alignment_record.strand = "-" if alignment_record.strand == "+" \
                    else "+"
                alignment_record.start_coord, alignment_record.end_coord = \
                    alignment_record.end_coord, alignment_record.start_coord
                alignment_record.seq = "".join(self.complement[base] for base in\
                        reversed(alignment_record.seq))

        except IndexError:
            print("Unable to format this line %s" % result)


class ReadAligner():

    def __init__(self, config, reads):
        self.config = config
        self.reads = reads

        self.formatter = RecordFormatter(
            self.config.transcript_bed,
            self.config.transcript_parents,
            self.config.exon_coords
        )

    def align_to_sam(self, ref, coding_transcript=False):

        bt_params = [
            # original version was bowtie, not bowtie2  
             "bowtie",
            #"bowtie2",
            "-f",
            "-v 0",
            "-S",
            # "--trimm5 4",
            # "--trim3 4",
            "--all",
            "--best",
            "--strata",
            ref,       # example: data/refs/bt/genome
            self.reads # example: data/tmp/WT_Cont_1_GGCTAC_clipped_collapsed.fq
        ]

        print("Running bowtie with the following parameters.")
        print("%% %s" % (' '.join(str(p) for p in bt_params)))

        align = subprocess.Popen(
            bt_params,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        std_out, std_err = align.communicate()

        filename = "bowtie_sam"
        if coding_transcript:
            filename = "bowtie_coding_sam"

        bowtie_intermediate = open(filename,"wb")
        for result in std_out.splitlines():
            bowtie_intermediate.write(result)
            bowtie_intermediate.write(str.encode("\n"))


    def align_to(self, ref, coding_transcript=False):

        bt_params = [
            # original version was bowtie, not bowtie2  
             "bowtie",
            #"bowtie2",
            "-f",
            "-v 0",
            # "--trimm5 4",
            # "--trim3 4",
            "--all",
            "--best",
            "--strata",
            ref,       # example: data/refs/bt/genome
            self.reads # example: data/tmp/WT_Cont_1_GGCTAC_clipped_collapsed.fq
        ]

        print("Running bowtie with the following parameters.")
        print("%% %s" % (' '.join(str(p) for p in bt_params)))

        align = subprocess.Popen(
            bt_params,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        std_out, std_err = align.communicate()

        #bowtie_intermediate = open("bowtie_intermediate","wb")
        for result in std_out.splitlines():
            #bowtie_intermediate.write(result)
            #bowtie_intermediate.write(str.encode("\n"))

            # result = result.strip()
            result = result.strip().decode("utf-8")
            if result[0] == "@" :
                continue

            try:
                # original
                # bt_record = BowtieRecord(*result.splits("\t"))

                # trang's attempt
                # bt_record = BowtieRecord(*result("\t"))

                # troy
                bt_record = BowtieRecord(*result.split("\t"))

                if bt_record.ref == "*":
                    continue

                self.formatter.format_record(bt_record, coding_transcript, result)
                yield bt_record.to_alignment()
            except TypeError:
                print("Failed to process this line: %s\n" % result)

        print(std_err)

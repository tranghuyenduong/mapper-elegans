# mapper-elegans
First create empty data folder that has 5 empty folders name: reads, alignments, ws253, refs, tmp
You will need to download ws253 genome in .fa format, gff2, gff3, gtf files from wormbase. The gtf file on wormbase might need some minor modification before but you would be able to identify the "error" once running thr "ref.py", then modify according to the instruction. 
I do have the "ready to go" gtf, ws253.fa, gff2, gff3 files if you want to be transfer it to you via external hardrive, since github doenst accept data transfering
Now you can run ref.py: the empty folder will be filled up with data once you run ref.py. 
One you finished running refs.py you should get a data set from GEO or any small RNA data in .fa format, deposite this data set into "reads" folder in "data" folder
identify the barcode of the sequencing data set and change it acordingly in setting.py: template = "AGATCGGAAGAGCACACGTCTGAACTCC"
Then run map.py
The alignment data should be populated in the mapper_elegans folder as 8 different folder: 4 for finnal alignment result which included: intron_inter, intron_exon, exon_exon, exon, and 4 for unfilter of miRNA, piRNA result which included: unfilter...

#! /bin/bash
WD="/user_data/men/sepseq/dummy114/"
cd $WD
mkdir $WD/2022_10_24_Barcode_count/
OutD=$WD/2022_10_24_Barcode_count/

raw_data=$WD/data/lib114dum/20221020_1339_X1_FAV42340_d58b518a/fastq_pass/

single_barcodes="/user_data/men/sepseq/dummy114/Single_barcodes_rev_for.fasta"

sample_name="dummy114"

mkdir $OutD/$sample_name
sample_dir=$OutD/$sample_name
cd $sample_dir
gunzip $raw_data/*.fastq.gz

cat $raw_data/*.fastq > $sample_dir/passed_reads.fastq

module load seqtk/1.3-foss-2018a
seqtk seq -a $sample_dir/passed_reads.fastq > $sample_dir/$sample_name.fasta
module purge

module load SeqKit/2.0.0
seqkit fx2tab -nl $sample_dir/passed_reads.fastq > $sample_dir/$sample_name\_read_lengths.tsv
module purge

mkdir $sample_dir/index
cd $sample_dir/index
module load Bowtie2/2.4.2-foss-2020b
bowtie2-build -f $sample_dir/$sample_name.fasta $sample_name --threads 50

bowtie2 -x $sample_name \
-U $single_barcodes \
-f \
--local \
-D 30 -R 4 -N 0 -L 5 -i S,1,0.25 \
-p 70 \
--norc \
-a \
--mp 2,2 \
--score-min C,31 \
--no-hd \
--met-file report_$sample_name.csv \
-S $sample_dir/$sample_name.sam 

module purge

#Evt. lav matrix i R med reads i rows og barcodes i columns med antal fundne barcodes i reads som værdier.

#Strategi for filtrering af resultater i R:
#	Der skal enten: 
#		1. Findes et set af to ens barcodes uden barcodes imellem eller 
#		2. Findes én barcode med meget høj alignment score. 


#Status 02/09: Nyt script for alle sekventerings runs kører på serveren.

#Include something in R about mean distance from end of a read.

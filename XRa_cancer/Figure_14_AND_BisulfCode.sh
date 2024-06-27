#Get the reference genome
#Download and get the human reference genome (the same as illumina EPIC has, to facilitate its combination).
mkdir reference_genome
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*' 
#I eliminated manually al the residual files that were installed, although it could just be done decompressing with the following line and then remove all the files except the concatenated reference genome.
for a in *.gz; do gunzip $a; done
#Concatenate all the files into one unique chromosome
cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrM.fa chrX.fa chrY.fa > hg19_ref_genome.fa

#Download the entry: GSE124027 

mkdir PersonalExample4
cd PersonalExample4/

mkdir reads_GSE124027

cd reads_GSE124027

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/005/SRR8346115/SRR8346115_1.fastq.gz -o SRR8346115_GSM3518982_C1_RRBS_Homo_sapiens_Bisulfite-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/005/SRR8346115/SRR8346115_2.fastq.gz -o SRR8346115_GSM3518982_C1_RRBS_Homo_sapiens_Bisulfite-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/005/SRR8346125/SRR8346125_1.fastq.gz -o SRR8346125_GSM3518987_A1_RRBS_Homo_sapiens_Bisulfite-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR834/005/SRR8346125/SRR8346125_2.fastq.gz -o SRR8346125_GSM3518987_A1_RRBS_Homo_sapiens_Bisulfite-Seq_2.fastq.gz

#Now untar everything

for file in *.gz; do
    # Untar the .gz file. The tool automatically eliminates the .gz files.
    gunzip "$file"
done

#The reads must be located each pair in a single folder. We need to create the new folders and move all the files to their respective folders.

	#First, the folders
folder_names=$(ls | cut -d "_" -f1 | uniq)
for new_folder in $folder_names 
do
	mkdir $new_folder
done

	#Now, we move the files
for i in *.fastq
do
	dest_folder=$(echo $i | cut -d "_" -f1)
	mv $i $dest_folder 
done

# Inspect the downloaded reads (FastQC). The fastqc codeline should be automatized to enter directly inside all the directories listed and process forward and reverse reads.

module load lang/Java/11.0.2
module load bio/FastQC/0.11.9-Java-11
module load vis/fontconfig/2.14.2-GCCcore-12.3.0

mkdir fastQC_OUT
fastqc -o fastQC_OUT/ -j /soft/modules/software/Java/11.0.2/bin/java *

#Load all the necessary tools

module load bio/SAMtools/1.15-GCC-11.2.0 
module load bio/Bowtie2/2.4.5-GCC-11.3.0 
module load lang/Java/15.0.1
alias bicycle="bicycle-1.8.2/cmd/bicycle" #Path to the algorithm
BOWTIE2=Bowtie2/2.4.5-GCC-11.3.0/bin #Path to the algorithm
SAMTOOLS=SAMtools/1.15-GCC-11.2.0/bin #Path to the algorithm

DATA_DIR=/PROJECTES/GENOMICS/aalegret/SequencingPipe/PersonalExample4
PROJECT_DIR=$DATA_DIR"/project_dir_GSE124027"
SAMPLES_DIR=$DATA_DIR"/reads_GSE124027"
REFERENCE_GENOME="/PROJECTES/GENOMICS/aalegret/SequencingPipe/PersonalExample/reference_genome/"

#Trimming procedure with trimgalore! Trim and filter the reads. Use TrimGalore! -o is the output directory where the results will be written (it writes a lot of different things, the remaining .fq files are the processed reads). I have to change the .fq for .fastq!

cd /PROJECTES/GENOMICS/aalegret/SequencingPipe/PersonalExample4/reads_GSE124027/ 
#Enter inside each folder, trim galore, enter in the folder, move the files and change the name. THE LOOP CODE DOES NOT WORK! THE TRIMGALORE FUNCTION YES.
for dir in SRR*
do 
	cd $dir
	trim_galore --paired --phred33 --fastqc -o trim_galore_out/ *_1.fastq *_2.fastq
	cd trim_galore_out/
	for file in *.fq
	do
		name_file=$(basename $file .fq)
		mv $file ../$name_file".fastq"
	done
	cd ../..
done


trim_galore --paired --phred33 --fastqc -o trim_galore_out/ SRR8346125_GSM3518987_A1_RRBS_Homo_sapiens_Bisulfite-Seq_1.fastq SRR8346125_GSM3518987_A1_RRBS_Homo_sapiens_Bisulfite-Seq_2.fastq

#Step 1.- Create the project. Most of the reads are paired end, and we need to specify the forward and the program will automatically detect both. Indicate eiher _1.fastq or _R1.fastq depending on if we have changed the name or not for previous preprocessing.

bicycle create-project -p $PROJECT_DIR -r $REFERENCE_GENOME -f $SAMPLES_DIR -b2 $BOWTIE2 -s $SAMTOOLS --paired-mate1-regexp _1.fastq

#Step 2 and 3 only if it is necessary (de novo reference genome).

#Step 2: Reference genome bisulfitation (Watson and Crick genomes). OUTPUT in $REFERENCE_GENOME.

bicycle reference-bisulfitation -p $PROJECT_DIR

#Step 3: Genomic index creation. OUTPUT in $REFERENCE_GENOME. -t is the number of threads.

bicycle reference-index -p $PROJECT_DIR -t 2

#Step 4.- Alignment
bicycle align -p $PROJECT_DIR -t 4 --bowtie2-quals phred33

# Post alignemnt quality control. First the sam files must be converted to bam files (the -O is the output, choose a name to it, and the sam file is the input). Three steps must be done: conversion SAM --> BAM format, sort BAM file and INDEX bam file. For all samtools files, the option -@ n can be used to modify the threads used (use -@ 4 in the cluster).

samtools view -Sb -o WATSON.bam bisulfited_CT_SRR7508945_against_hg19_NO_Y_ref_genome.fa_WATSON.sam
samtools sort -O bam -o sorted_WATSON.bam WATSON.bam
samtools index sorted_WATSON.bam

#Now we use qualymap. This tool is installed in "SeqPipe" conda environment. Download all the folder to see correctly the plots and the .html created.


#Step 5.- Methylation calling
bicycle analyze-methylation -p $PROJECT_DIR -t 4 -n 2 --remove-ambiguous --only-with-one-alignment

#Reorder the results within the output directory to use qualimap. bam files are generated after the methylation calling.
cd /PROJECTES/GENOMICS/aalegret/SequencingPipe/PersonalExample4/project_dir_GSE124027/output/

for file in *.sam.sorted.sam 
do
	name_file=$(basename "$file" .sam.sorted.sam)
	mkdir $name_file
	mv $name_file.* $name_file
done

#Now, we enter each folder and do the qualimap inside each. Rmeember to enter inside the conda environment SeqPipe!
for dir in bisulfited*
do
	echo $dir
	cd $dir
	qualimap bamqc -bam *.bam -outdir "qualimap_results_"$dir -c -gd HUMAN
	cd ..
done

#!/bin/bash
#Downloads fastq files from SRA, trims, maps and generates plotfiles

#ASSUMES PAIRED END ILLUMINA WITH TRUSEQ ADAPTERS

#REQUIREMENTS
#path to genome in fasta format
#SRX folder path
#SRR file
#output folder

SRA=$1
GENOME=$2

if [ ! -d trimmed ]; then
mkdir trimmed
fi

fastq-dump --split-3 ${SRA} &&

java -jar ~/inst/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 ${SRA}_1.fastq ${SRA}_2.fastq trimmed/${SRA}_1_paired.fastq trimmed/${SRA}_1_unpaired.fastq trimmed/${SRA}_2_paired.fastq trimmed/${SRA}_2_unpaired.fastq ILLUMINACLIP:/home/beth/inst/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


#Build index of genome if necessary

if [ ! -d index ]; then
bowtie2-build ${GENOME}.fasta ${GENOME} &&
mkdir index && mv *.bt2* index/
fi

#mapping

bowtie2 -p 20 -x index/${GENOME} -1 trimmed/${SRA}_1_paired.fastq -2 trimmed/${SRA}_2_paired.fastq -S ${SRA}.sam

samtools view -bS -@ 10 ${SRA}.sam > ${SRA}.bam
samtools sort -@ 10 ${SRA}.bam > ${SRA}.sorted.bam

#fastq2plot takes sorted BAM file and genome and generates plotfiles for artemis

fastq2plot3.sh ${SRA} 

rm *.bam
rm *.mpileup
rm *.bai

if [ ! -d fastq ]; then
mkdir fastq
fi

mv *.fastq fastq/

#To-do
#add filename/directory and install checks
#add opts for directory outputs
#make functions to call with opts
#write a small readme

#!/bin/sh

#hello, this will be my bash script
#python project.py
#star --runMode genomeGenerate --genomeDir /home/arosa/project/ref_genome/star_indexed --genomeFastaFiles /home/arosa/project/ref_genome/GRCh38_r90.all.fa --runThreadN 8 --sjdbGTFfile /home/arosa/project/ref_genome/annotations/Homo_sapiens.GRCh38.90.gtf
#STAR  --runThreadN 12 --genomeDir /proj/sllstore2017067/project6/STAR_indexed_genome --readFilesIn /proj/sllstore2017067/project6/Andre/project/datasets/sample2.fastaq
#bowtie2 -x bwtie2 -U /proj/sllstore2017067/project6/Andre/datasets/control/control.fastq -S control_bwtie2_aligned.sam
#samtools faidx /proj/sllstore2017067/project6/Andre/ref_genome/GRCh38_r90.all.fa
#samtools import /proj/sllstore2017067/project6/Andre/ref_genome/GRCh38_r90.all.fa.fai /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_aligned.sam /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_aligned.bam
#samtools sort /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_aligned.bam -o /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_sorted.bam
#bamToBed -i /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_sorted.bam > /proj/sllstore2017067/project6/Andre/datasets/control/control_bwtie2_sorted.bed
#module add MACS

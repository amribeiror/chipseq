#!/bin/sh

#hello, this will be my bash script

module add bioinfo-tools
module add FastQC

fastqc -o /proj/sllstore2017067/project6/andre/datasets/control/trimmed /proj/sllstore2017067/project6/andre/datasets/control/control_trimmed_new.fastq 

#star --runMode genomeGenerate --genomeDir /home/arosa/project/ref_genome/star_indexed --genomeFastaFiles /home/arosa/project/ref_genome/GRCh38_r90.all.fa --runThreadN 8 --sjdbGTFfile /home/arosa/project/ref_genome/annotations/Homo_sapiens.GRCh38.90.gtf
#STAR  --runThreadN 12 --genomeDir /proj/sllstore2017067/project6/STAR_indexed_genome --readFilesIn /proj/sllstore2017067/project6/Andre/project/datasets/sample2.fastaq
#module add cutadapt
#cutadapt -g GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACGATCTCGTAT -o /proj/sllstore2017067/project6/andre/datasets/control/control_trimmed_new.fastq /proj/sllstore2017067/project6/andre/datasets/control/control.fastq

#PASTE COMMAND TO DONWLOAD GENOME HERE!!

module add bowtie2

#ALIGNS READS TO INDEXED REFERENCE GENOME
#for control
bowtie2 -x bowtie2_index_refgen_GRCh38 -U /proj/sllstore2017067/project6/andre/datasets/control/control.fastq -S control_andre_bw2_aligned.sam

#for sample1
bowtie2 -x bowtie2_index_refgen_GRCh38 -U /proj/sllstore2017067/project6/andre/datasets/sample1/sample1.fastaq -S sample1_andre_bw2_aligned.sam

#for sample2
bowtie2 -x bowtie2_index_refgen_GRCh38 -U /proj/sllstore2017067/project6/andre/datasets/sample2/sample2.fastaq -S sample2_andre_bw2_aligned.sam

module add samtools

#REMOVES UNMAPPED READS AND KEEPS MAPPED READS
#for control
samtools view -F4 /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.sam -b -o /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped3.bam

#for sample1
samtools view -F4 /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_aligned.sam -b -o /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped3.bam

#for sample2
samtools view -F4 -b /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_aligned.sam -o /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped3.bam
samtools flagstat /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped3.bam
samtools flagstat /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped3.bam

#STATISTICS WITH SAMVIEW
#counts UNmapped reads on SAM file
samtools view -f4 /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.sam
#counts mapped reads on SAM file
samtools view -F4 /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.sam
#counts UNmapped reads on BAM file
samtools view -f 0x04 -b /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.bam 
#counts UNmapped reads on mapped BAM file, should be 0
samtools view -f4 /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped.bam
#counts mapped reads on BAM file, should be the same as in the SAM file
samtools view -F4 /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped.bam

samtools faidx /proj/sllstore2017067/project6/ref_genome/all_chr_ref_GRCh38.fa

samtools import /proj/sllstore2017067/project6/tools/samtools/all_chr_ref_GRCh38.fa.fai /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.sam /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.bam
samtools import /proj/sllstore2017067/project6/tools/samtools/all_chr_ref_GRCh38.fa.fai /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_aligned.sam /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_aligned.bam
samtools import /proj/sllstore2017067/project6/tools/samtools/all_chr_ref_GRCh38.fa.fai /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_aligned.sam /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_aligned.bam

samtools sort /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_aligned.bam -o /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_sorted.bam
samtools sort /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_aligned.bam -o /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_sorted.bam
samtools sort /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_aligned.bam -o /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_sorted.bam
bamToBed -i /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_sorted.bam > /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_sorted.bed
bamToBed -i /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_sorted.bam > /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_sorted.bed
bamToBed -i /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_sorted.bam > /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_sorted.bed

#ACCESSES REPRODUCIBILITY OF DATA
#multiBamSummary requires index BAM files. This indexes both sample1 and sample2 .BAM files.
#sample1
samtools sort -o /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped_sorted.bam /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped.bam
samtools index /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped_sorted.bam /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped.bai
#javi's sample1
samtools index /proj/sllstore2017067/project6/Javi/data/alignments/sample1_mapped_bowtie2.sorted.duplicates_filtered.bam /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_mapped_bowtie2.sorted.duplicates_filtered.bai

#sample2
samtools sort -o /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped_sorted.bam /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped.bam
samtools index /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped_sorted.bam  /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped.bai
samtools index /proj/sllstore2017067/project6/Javi/data/alignments/sample2_mapped_bowtie2.sorted.duplicates_filtered.bam /proj/sllstore2017067/project6/Javi/data/alignments/sample2_mapped_bowtie2.sorted.duplicates_filtered.bai

multiBamSummary bins -b /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped_sorted.bam /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped_sorted.bam -out /proj/sllstore2017067/project6/andre/tools/deeptools/correlation.npz
#javi's correlation
multiBamSummary bins -b /proj/sllstore2017067/project6/Javi/data/alignments/sample1_mapped_bowtie2.sorted.duplicates_filtered.bam /proj/sllstore2017067/project6/Javi/data/alignments/sample2_mapped_bowtie2.sorted.duplicates_filtered.bam -out /proj/sllstore2017067/project6/andre/tools/deeptools/correlation_javi.npz

module add MACS

#RUNS MACS2
#For sample1 against control
macs2 callpeak -t /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped.bam -c /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped.bam -f BAM -n sample1 --outdir /proj/sllstore2017067/project6/andre/datasets/callpeaks 

#For sample2 against control
macs2 callpeak -t /proj/sllstore2017067/project6/andre/datasets/sample2/sample2_andre_bw2_mapped.bam -c /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped.bam -f BAM -n sample2 --outdir /proj/sllstore2017067/project6/andre/datasets/callpeaks

#RUNS HOMER(NOW ON UPPMAX)
export PATH=$PATH:/proj/sllstore2017067/project6/Evgeny/scripts/bin
/proj/sllstore2017067/project6/Evgeny/scripts/bin/proj/sllstore2017067/project6/Evgeny/scripts/bin/findPeaks /proj/sllstore2017067/project6/andre/datasets/sample1/sample1_andre_bw2_mapped.bam -style factor -o auto -i /proj/sllstore2017067/project6/andre/datasets/control/control_andre_bw2_mapped.bam
perl /proj/sllstore2017067/project6/Evgeny/scripts/bin/findMotifsGenome.pl /proj/sllstore2017067/project6/andre/datasets/callpeaks/sample1_summits_chr_plus.bed hg38 /proj/sllstore2017067/project6/andre/datasets/motifs/sample1_Upp

perl /proj/sllstore2017067/project6/Evgeny/scripts/bin/findMotifsGenome.pl /proj/sllstore2017067/project6/andre/datasets/callpeaks/sample2_summits_chr_plus.bed hg38 /proj/sllstore2017067/project6/andre/datasets/motifs/sample2_Upp

/proj/sllstore2017067/project6/Evgeny/scripts/bin/annotatePeaks.pl /proj/sllstore2017067/project6/andre/datasets/callpeaks/sample1_summits_chr_plus.bed hg38  > /proj/sllstore2017067/project6/andre/datasets/annotatedpeaks/sample1

/proj/sllstore2017067/project6/Evgeny/scripts/bin/annotatePeaks.pl /proj/sllstore2017067/project6/andre/datasets/callpeaks/sample2_summits_chr_plus.bed hg38  > /proj/sllstore2017067/project6/andre/datasets/annotatedpeaks/sample2

bedtools genomecov -ibam bedtools genomecov -ibam /proj/sllstore2017067/projects6/andre/datasets/sample1/sample1_andre_bw2_mapped_sorted.bam -bg > /proj/sllstore2017067/projects6/andre/datasets/sample1/sample1_andre_bw2_mapped_sorted.bedgraph

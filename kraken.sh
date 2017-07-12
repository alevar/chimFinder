#!/usr/bin/env bash

inputDir=$1
file=$2
outputDir=$3
krakenDB=$4
hivDB=$5
humanDB=$6

mkdir -p ${outputDir}
mkdir -p ${outputDir}/krakenOut
mkdir -p ${outputDir}/krakenOut/selected
mkdir -p ${outputDir}/localAlignments
mkdir -p ${outputDir}/fullAlignments
mkdir -p ${outputDir}/consensusHIV
mkdir -p ${outputDir}/hisat

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz

echo RUNNING KRAKEN: $sample
kraken --threads 10 --fastq-input --gzip-compressed --output ${outputDir}/krakenOut/${sample}.kraken --paired ${inputDir}/${sampleR1Base} ${inputDir}/${sampleR2Base} --db ${krakenDB}
kraken-report --db ${krakenDB} ${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/${sample}.report

#extract chimeric read names
echo EXTRACTING READS
awk -F '\t' 'match($5,/9606*/) && match($5,/11676*/) {print $2}' ${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/selected/${sample}.chim

#parse names and extract reads into fastq
touch ${outputDir}/krakenOut/selected/${sample}_R1.fastq
touch ${outputDir}/krakenOut/selected/${sample}_R2.fastq
zcat ${inputDir}/${sampleR1Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R1.fastq
zcat ${inputDir}/${sampleR2Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R2.fastq

#align with bowtie2 to hg38
echo ALIGNING TO HG38
bowtie2 --no-unal --local --phred33 -p 10 -x ${humanDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.chim.hum.sam
#align with bowtie2 to hiv89.6
echo ALIGNING TO HIV89.6
bowtie2 --no-unal --local --phred33 -p 10 -x ${hivDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.chim.hiv.sam

echo ALIGNING ALL HUMAN READS
bowtie2 --very-sensitive --no-unal --local --phred33 -p 12 -x ${humanDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}.full.hum.sam
samtools view -S -@ 10 -b ${outputDir}/fullAlignments/${sample}.full.hum.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}.full.hum.bam -
samtools index -@ 10 ${outputDir}/fullAlignments/${sample}.full.hum.bam

echo ALIGNING ALL HIV READS
bowtie2 --very-sensitive --no-unal --local --phred33 -p 10 -x ${hivDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}.full.hiv.sam
samtools view -S -@ 10 -b ${outputDir}/fullAlignments/${sample}.full.hiv.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}.full.hiv.bam -
samtools index -@ 10 ${outputDir}/fullAlignments/${sample}.full.hiv.bam

#build consensus hiv sequence
echo BUILDING CONSENSUS SEQUENCE
samtools mpileup -uf ./refs/HIV-1_89.6_sequence.fa ${outputDir}/fullAlignments/${sample}.full.hiv.bam | bcftools call -c | vcfutils.pl vcf2fq > ${outputDir}/consensusHIV/${sample}.fq
echo CONVERT TO FASTA
# awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' ${outputDir}/consensusHIV/${sample}.fq 20 > ${outputDir}/consensusHIV/${sample}.fa
seqtk seq -a  ${outputDir}/consensusHIV/${sample}.fq > ${outputDir}/consensusHIV/${sample}.fa
echo BUILDING HISAT INDEX
hisat2-build -p 10 ${outputDir}/consensusHIV/${sample}.fa ${outputDir}/consensusHIV/${sample}
echo ALIGNING WITH hisat2
hisat2 -p 10 --no-unal --dta --phred33 -x ${outputDir}/consensusHIV/${sample} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} --novel-splicesite-outfile ${outputDir}/hisat/${sample}.junctions -S ${outputDir}/hisat/${sample}.sam
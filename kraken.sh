#!/usr/bin/env bash

inputDir=$1
file=$2
outputDir=$3
krakenDB=$4
hivDB=$5
humanDB=$6
annotation=$7
threads=$8

mkdir -p ${outputDir}
mkdir -p ${outputDir}/krakenOut
mkdir -p ${outputDir}/krakenOut/selected
mkdir -p ${outputDir}/localAlignments
mkdir -p ${outputDir}/fullAlignments
mkdir -p ${outputDir}/consensusHIV
mkdir -p ${outputDir}/hisat
mkdir -p ${outputDir}/splices
mkdir -p ${outputDir}/beds

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz

echo "$(date)"
SECONDS=0
echo RUNNING KRAKEN: $sample
kraken --preload --threads ${threads} --fastq-input --gzip-compressed --output ${outputDir}/krakenOut/${sample}.kraken --paired ${inputDir}/${sampleR1Base} ${inputDir}/${sampleR2Base} --db ${krakenDB}
kraken-report --db ${krakenDB} ${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/${sample}.report
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}

SECONDS=0
TOTAL_TIME=0
# extract chimeric read names
echo EXTRACTING READS
awk -F '\t' 'match($5,/9606*/) && match($5,/11676*/) {print $2}' ${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/selected/${sample}.chim

#parse names and extract reads into fastq
touch ${outputDir}/krakenOut/selected/${sample}_R1.fastq
touch ${outputDir}/krakenOut/selected/${sample}_R2.fastq
zcat ${inputDir}/${sampleR1Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R1.fastq
zcat ${inputDir}/${sampleR2Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R2.fastq

DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
# align with bowtie2 to hg38
echo ALIGNING TO HG38
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.chim.hum.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
#align with bowtie2 to hiv89.6
echo ALIGNING TO HIV89.6
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.chim.hiv.sam

DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
echo ALIGNING ALL HUMAN READS
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p 12 -x ${humanDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}.full.hum.sam

bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}.full.hiv.sam

samtools view -S -@ ${threads} -b ${outputDir}/fullAlignments/${sample}.full.hum.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}.full.hum.bam -
samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}.full.hum.bam
samtools view -S -@ ${threads} -b ${outputDir}/fullAlignments/${sample}.full.hiv.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}.full.hiv.bam -
samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}.full.hiv.bam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
echo MAPPING AGAINST ANNOTATION
./get_hum_hiv.pl ${outputDir}/fullAlignments/${sample}.full.hum.sam ${outputDir}/fullAlignments/${sample}.full.hiv.sam ${annotation} > ${outputDir}/splices/${sample}.full.txt
./get_hum_hiv.pl ${outputDir}/localAlignments/${sample}.chim.hum.sam ${outputDir}/localAlignments/${sample}.chim.hiv.sam ${annotation} > ${outputDir}/splices/${sample}.chim.txt
#
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
#build consensus hiv sequence
echo BUILDING CONSENSUS SEQUENCE
samtools mpileup -uf /ccb/salz7-data/genomes/HIV1/HIV1_genomes.fa ${outputDir}/fullAlignments/${sample}.full.hiv.bam | bcftools call -c | vcfutils.pl vcf2fq > ${outputDir}/consensusHIV/${sample}.fq
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo CONVERTING TO FASTA
# awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' ${outputDir}/consensusHIV/${sample}.fq 20 > ${outputDir}/consensusHIV/${sample}.fa
seqtk seq -a  ${outputDir}/consensusHIV/${sample}.fq > ${outputDir}/consensusHIV/${sample}.fa
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo BUILDING HISAT INDEX
hisat2-build -p ${threads} ${outputDir}/consensusHIV/${sample}.fa ${outputDir}/consensusHIV/${sample}
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo ALIGNING WITH hisat2
hisat2 -p ${threads} --no-unal --dta --phred33 -x ${outputDir}/consensusHIV/${sample} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} --novel-splicesite-outfile ${outputDir}/hisat/${sample}.junctions -S ${outputDir}/hisat/${sample}.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
SECONDS=0

T_DUR="$(($TOTAL_TIME / 60)) minutes and $(($TOTAL_TIME % 60)) seconds"
echo TOTAL TIME ELAPSED: ${T_DUR}

echo "$(date)"

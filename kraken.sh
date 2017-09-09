#!/usr/bin/env bash

inputDir=$1
file=$2
outputDir=$3
krakenDB=$4
hivDB=$5
humanDB=$6
annotation=$7
threads=$8
baseEnd=$9

mkdir -p ${outputDir}
mkdir -p ${outputDir}/krakenOut
mkdir -p ${outputDir}/krakenOut/selected
mkdir -p ${outputDir}/localAlignments
mkdir -p ${outputDir}/fullAlignments
mkdir -p ${outputDir}/consensusHIV
mkdir -p ${outputDir}/hisat
mkdir -p ${outputDir}/splices
mkdir -p ${outputDir}/beds
mkdir -p ${outputDir}/tempF
mkdir -p ${outputDir}/tmp

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz

echo "$(date)"
SECONDS=0
echo RUNNING KRAKEN: $sample
kraken --preload --threads ${threads} --fastq-input --gzip-compressed --output ${outputDir}/krakenOut/${sample}${baseEnd}.kraken --paired ${inputDir}/${sampleR1Base} ${inputDir}/${sampleR2Base} --db ${krakenDB}
kraken-report --db ${krakenDB} ${outputDir}/krakenOut/${sample}${baseEnd}.kraken > ${outputDir}/krakenOut/${sample}${baseEnd}.report
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}

SECONDS=0
TOTAL_TIME=0
# extract chimeric read names
echo EXTRACTING READS
awk -F '\t' 'match($5,/9606*/) && match($5,/11676*/) {print $2}' ${outputDir}/krakenOut/${sample}${baseEnd}.kraken > ${outputDir}/krakenOut/selected/${sample}${baseEnd}.chim

#parse names and extract reads into fastq
touch ${outputDir}/krakenOut/selected/${sample}_R1${baseEnd}.fastq
touch ${outputDir}/krakenOut/selected/${sample}_R2${baseEnd}.fastq
zcat ${inputDir}/${sampleR1Base} > ${outputDir}/tempF/${sample}_R1${baseEnd}.fastq
zcat ${inputDir}/${sampleR2Base} > ${outputDir}/tempF/${sample}_R2${baseEnd}.fastq
grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}${baseEnd}.chim ${outputDir}/tempF/${sample}_R1${baseEnd}.fastq > ${outputDir}/krakenOut/selected/${sample}_R1${baseEnd}.fastq
grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}${baseEnd}.chim ${outputDir}/tempF/${sample}_R2${baseEnd}.fastq > ${outputDir}/krakenOut/selected/${sample}_R2${baseEnd}.fastq

DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
# align with bowtie2 to hg38
echo ALIGNING TO HG38
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1${baseEnd}.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2${baseEnd}.fastq -S ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
#align with bowtie2 to hiv89.6
echo ALIGNING TO HIV89.6
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1${baseEnd}.fastq -2 ${outputDir}/krakenOut/selected/${sample}_R2${baseEnd}.fastq -S ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
echo ALIGNING ALL HUMAN READS
bowtie2 --very-sensitive-local --no-unal --local --phred33 -p 12 -x ${humanDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.sam

bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

SECONDS=0
echo REMOVING DUPLICATES

samtools view -S -@ ${threads} -b ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.bam -
samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.bam
# java -Xmx4g -jar /ccb/salz7-data/sw/packages/picard-tools-1.119/MarkDuplicates.jar INPUT=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.bam OUTPUT=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.no_dup.rem.bam METRICS_FILE=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=${outputDir}/tmp
# samtools view -S -@ ${threads} -h -o ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.sam ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.no_dup.rem.bam
# samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.bam

samtools view -S -@ ${threads} -b ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.sam | samtools sort -o ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.bam -
samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.bam
# java -Xmx4g -jar /ccb/salz7-data/sw/packages/picard-tools-1.119/MarkDuplicates.jar INPUT=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.bam OUTPUT=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.no_dup.rem.bam METRICS_FILE=${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=${outputDir}/tmp
# samtools view -S -@ ${threads} -h -o ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.sam ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.no_dup.rem.bam
# samtools index -@ ${threads} ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.bam

samtools view -S -@ ${threads} -b ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.sam | samtools sort -o ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.bam -
samtools index -@ ${threads} ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.bam
# java -Xmx4g -jar /ccb/salz7-data/sw/packages/picard-tools-1.119/MarkDuplicates.jar INPUT=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.bam OUTPUT=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.no_dup.rem.bam METRICS_FILE=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=${outputDir}/tmp
# samtools view -S -@ ${threads} -h -o ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.sam ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.no_dup.rem.bam
# samtools index -@ ${threads} ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.bam

samtools view -S -@ ${threads} -b ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.sam | samtools sort -o ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.bam -
samtools index -@ ${threads} ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.bam
# java -Xmx4g -jar /ccb/salz7-data/sw/packages/picard-tools-1.119/MarkDuplicates.jar INPUT=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.bam OUTPUT=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.no_dup.rem.bam METRICS_FILE=${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=${outputDir}/tmp
# samtools view -S -@ ${threads} -h -o ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.sam ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.no_dup.rem.bam
# samtools index -@ ${threads} ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.bam

DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME+${SECONDS}))

SECONDS=0
echo MAPPING AGAINST ANNOTATION
./get_hum_hiv.pl ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hum.sam ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.sam ${annotation} ${outputDir}/tmp > ${outputDir}/splices/${sample}${baseEnd}.full.txt
./get_hum_hiv.pl ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hum.sam ${outputDir}/localAlignments/${sample}${baseEnd}.chim.hiv.sam ${annotation} ${outputDir}/tmp > ${outputDir}/splices/${sample}${baseEnd}.chim.txt
#
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
#build consensus hiv sequence
echo BUILDING CONSENSUS SEQUENCE
samtools mpileup -uf ${hivDB}.fa ${outputDir}/fullAlignments/${sample}${baseEnd}.full.hiv.bam | bcftools call -c | vcfutils.pl vcf2fq > ${outputDir}/consensusHIV/${sample}${baseEnd}.fq
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo CONVERTING TO FASTA
# awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' ${outputDir}/consensusHIV/${sample}.fq 20 > ${outputDir}/consensusHIV/${sample}.fa
seqtk seq -a  ${outputDir}/consensusHIV/${sample}${baseEnd}.fq > ${outputDir}/consensusHIV/${sample}${baseEnd}.fa
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo BUILDING HISAT INDEX
hisat2-build -p ${threads} ${outputDir}/consensusHIV/${sample}${baseEnd}.fa ${outputDir}/consensusHIV/${sample}${baseEnd}
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
#
SECONDS=0
echo ALIGNING WITH hisat2
hisat2 -p ${threads} --no-unal --dta --phred33 -x ${outputDir}/consensusHIV/${sample}${baseEnd} -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} --novel-splicesite-outfile ${outputDir}/hisat/${sample}${baseEnd}.junctions -S ${outputDir}/hisat/${sample}${baseEnd}.sam
DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
echo DONE IN ${DUR}
TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
SECONDS=0

T_DUR="$(($TOTAL_TIME / 60)) minutes and $(($TOTAL_TIME % 60)) seconds"
echo TOTAL TIME ELAPSED: ${T_DUR}

echo "$(date)"

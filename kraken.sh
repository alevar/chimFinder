#!/usr/bin/env bash

inputDir=$1
outputDir=$2
krakenDB=$3
hivDB=$4
humanDB=$5

mkdir -p ${outputDir}/krakenOut
mkdir -p ${outputDir}/krakenOut/selected
mkdir -p ${outputDir}/localAlignments
mkdir -p ${outputDir}/fullAlignments

for file in $inputDir/*R1_001.fastq.gz ; do
    sampleR1Base=$(basename ${file})
    sampleR1="${sampleR1Base%.*.*}"
    sample="${sampleR1%_R1*}"
    sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz
    echo $sampleR2Base
    kraken --threads 10 --fastq-input --gzip-compressed --output ./${outputDir}/krakenOut/${sample}.kraken --paired ${inputDir}/${sampleR1Base} ${inputDir}/${sampleR2Base} --db ${krakenDB}
    kraken-report --db ${krakenDB} ${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/${sample}.report

    #extract chimeric read names
    awk -F '\t' 'match($5,/9606*/) && match($5,/11676*/) {print $2}' ./${outputDir}/krakenOut/${sample}.kraken > ${outputDir}/krakenOut/selected/${sample}.chim

    #parse names and extract reads into fastq
    touch ${outputDir}/krakenOut/selected/${sample}_R1.fastq
    touch ${outputDir}/krakenOut/selected/${sample}_R2.fastq
    zcat ${inputDir}/${sampleR1Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R1.fastq
    zcat ${inputDir}/${sampleR2Base} | grep -A 3 --no-group-separator -Ff ${outputDir}/krakenOut/selected/${sample}.chim - > ${outputDir}/krakenOut/selected/${sample}_R2.fastq

    #align with bowtie2 to hg38
    bowtie2 --no-unal --local --phred33 -k 20 -x ${humanDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.hum.chim.sam
    #align with bowtie2 to hiv89.6
    bowtie2 --no-unal --local --phred33 -k 20 -x ${hivDB} -1 ${outputDir}/krakenOut/selected/${sample}_R1.fastq ${outputDir}/krakenOut/selected/${sample}_R2.fastq -S ${outputDir}/localAlignments/${sample}.hiv.chim.sam

    #Make a full alignment of HIV89.6
    hisat2 -x ${hivDB} -p 8 --dta -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ./fullAlignments/${sample}.hiv.sam
    #Make a full alignment of human
    hisat2 -x ${humanDB}.2 -p 8 --dta -1 ${inputDir}/${sampleR1Base} -2 ${inputDir}/${sampleR2Base} -S ./fullAlignments/${sample}.hum.sam

    #or
    bowtie2 --no-unal --local --phred33 -x ../refs/bowtie2/hg38 -1 ./data/108PI_pos_S2_R1_001.fastq.gz -2 ./data/108PI_pos_S2_R2_001.fastq.gz -S ./fullAlignments/108PI_pos_S2.hum.all.sam
    bowtie2 --no-unal --local --phred33 -x ../refs/bowtie2/hg38 -1 ./data/108PI_pos_S2_R1_001.fastq.gz -2 ./data/108PI_pos_S2_R2_001.fastq.gz -S ./fullAlignments/108PI_pos_S2.hum.all.sam
done
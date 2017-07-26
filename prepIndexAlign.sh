#!/usr/bin/ebv bash

outDir=$1
sample=$2
inDir=$3

mkdir -p ${outDir}/splicing
mkdir -p ${outDir}/splicing/${sample}

sampleOut=${outDir}/splicing/${sample}
spliceRefs=${outDir}/spliceRefs/${sample}

# extract hiv read names
awk -F '\t' '{print $1}' ${outDir}/fullAlignments/${sample}.full.hiv.sam >> ${sampleOut}/reads.txt

# extract previously extracted read names into fastq record ofall hiv reads
zcat ${inDir}/${sample}_R1_001.fastq.gz | grep -A 3 --no-group-separator -Ff ${sampleOut}/reads.txt - > ${sampleOut}/R1.fq
zcat ${inDir}/${sample}_R1_001.fastq.gz | grep -A 3 --no-group-separator -Ff ${sampleOut}/reads.txt - > ${sampleOut}/R2.fq
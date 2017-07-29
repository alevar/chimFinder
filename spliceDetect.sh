#!/usr/bin/env bash

outDir=$1
sample=$2
inDir=$3
pos=$4
posS=$5

sampleOut=${outDir}/splicing/${sample}
spliceRefs=${outDir}/spliceRefs/${sample}

# create bowtie2 index of the left and right excerpts from te reference
bowtie2-build ${spliceRefs}/i-1_${pos}.fa ${spliceRefs}/i-1_${pos}
bowtie2-build ${spliceRefs}/i-2_${pos}.fa ${spliceRefs}/i-2_${pos}
bowtie2-build ${spliceRefs}/u-${pos}_${posS}.fa ${spliceRefs}/u-${pos}_${posS}

bowtie2 --no-unal --local --phred33 -p 12 -x ${spliceRefs}/i-1_${pos} -1 ${sampleOut}/R1.fq -2 ${sampleOut}/R2.fq -S ${sampleOut}/i-1_${pos}.sam
bowtie2 --no-unal --local --phred33 -p 12 -x ${spliceRefs}/i-2_${pos} -1 ${sampleOut}/R1.fq -2 ${sampleOut}/R2.fq -S ${sampleOut}/i-2_${pos}.sam
bowtie2 --no-unal --local --phred33 -p 12 -x ${spliceRefs}/u-${pos}_${posS} -1 ${sampleOut}/R1.fq -2 ${sampleOut}/R2.fq -S ${sampleOut}/u-${pos}_${posS}.sam
bowtie2 --no-unal --phred33 -p 12 -x ${spliceRefs}/u-${pos}_${posS} -1 ${sampleOut}/R1.fq -2 ${sampleOut}/R2.fq -S ${sampleOut}/u-${pos}_${posS}.f.sam

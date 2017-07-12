#!/usr/bin/env bash

inputDir=$1
file=$2
pathChimReads=$3
outPath=$4

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz

zcat ${inputDir}/${sampleR1Base} | grep -A 3 --no-group-separator -Ff ${pathChimReads} - > ${outPath}
zcat ${inputDir}/${sampleR2Base} | grep -A 3 --no-group-separator -Ff ${pathChimReads} - > ${outPath}
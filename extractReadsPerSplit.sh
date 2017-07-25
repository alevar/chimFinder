#!/usr/bin/env bash

inputDir=$1
file=$2
pathChimReads=$3
outPath1=$4
outPath2=$5

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq
sampleR1Base="${sampleR1%_R1*}"_R1"${sampleR1##*_R1}".fastq

head -n -1 ${pathChimReads} >${pathChimReads}.t
rm ${pathChimReads}
mv ${pathChimReads}.t ${pathChimReads}

grep -A 3 --no-group-separator -Ff ${pathChimReads} ${inputDir}/${sampleR1Base} > ${outPath1}
grep -A 3 --no-group-separator -Ff ${pathChimReads} ${inputDir}/${sampleR2Base} > ${outPath2}
#!/usr/bin/env bash

outputDir=$1
file=$2
string=$3
outPath1=$4
outPath2=$5

sampleR1Base=$(basename ${file})
sampleR1="${sampleR1Base%.*.*}"
sample="${sampleR1%_R1*}"
sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq
sampleR1Base="${sampleR1%_R1*}"_R1"${sampleR1##*_R1}".fastq

# head -n -1 ${pathChimReads} >${pathChimReads}.t
# rm ${pathChimReads}
# mv ${pathChimReads}.t ${pathChimReads}

egrep -A 3 "${string}" ${outputDir}/${sampleR1Base} | seqtk seq -a - > ${outPath1}
egrep -A 3 "${string}" ${outputDir}/${sampleR2Base} | seqtk seq -a - > ${outPath2}
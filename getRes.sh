#!/usr/bin/env bash

inputDir=${1}

for file in ./out*/ ; do
        # echo ${file}
        sampleBase=$(basename ${file})
        sampleExt=${sampleBase: -2}
        # echo ${sampleExt}
        if [[ $sampleExt != "_1" && $sampleExt != "_2" ]]
        then
            sampleName=${sampleBase}
            echo $file
        fi
        # sampleR1Base=$(basename ${file})
        # sampleR1="${sampleR1Base%.*}"
        # sample="${sampleR1%_R1*}"
        # sampleR2Base="${sampleR1%_R1*}"_R2.fastq
        # baseEnd="${sampleR1##*_R1}"
        # echo ${sample}_R1${baseEnd}
#
        # dirP=${file%/*}
        # outputDir=outSample${dirP#*Sample}
        # echo ${dirP}
done

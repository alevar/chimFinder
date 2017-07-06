#!/usr/bin/env bash

outDir=$1

touch ./results.csv
echo "sample,krakenHiv%,krakenHum%,bowtie2HIV,bowtie2HUM,chimericReadsExtracted,chimericReadsFiltered" > ./results.csv

for file in ${outDir}/krakenOut/*.report ; do
    sampleBase=$(basename ${file})
    sample="${sampleBase%.*}"

    echo "${sample}"

    HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
    HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
    
    nReads="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
    if [ ! -f ${outDir}/localAlignments/${sample}.chim.csv ]; then
        bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        if [ ! -f ${outDir}/fullAlignments/${sample}.hiv.full.sam ]; then
            echo "${sample},${HIVP},${HUMP},0,${bowtieReadsHUM},${nReads},0" >> ./results.csv
        else
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.hiv.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        fi
        echo "${sample},${HIVP},${HUMP},${bowtieReadsHIV},${bowtieReadsHUM},${nReads},0" >> ./results.csv
    else
        postFilt="$(tail -n "$(wc -l ${outDir}/localAlignments/${sample}.chim.csv | awk -F ' ' '{print $1-1}')" ${outDir}/localAlignments/${sample}.chim.csv | awk -F "," '{print $2}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        if [ ! -f ${outDir}/fullAlignments/${sample}.hiv.full.sam ]; then
            echo "${sample},${HIVP},${HUMP},0,${bowtieReadsHUM},${nReads},${postFilt}" >> ./results.csv
        else
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.hiv.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        fi
        echo "${sample},${HIVP},${HUMP},${bowtieReadsHIV},${bowtieReadsHUM},${nReads},${postFilt}" >> ./results.csv
    fi
done

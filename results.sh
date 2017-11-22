#!/usr/bin/env bash

touch ./results2.csv
echo sequencingRun,R,sample,readsTotal,readsPreprocessed,numHIV_AllGenomes,numHIV_HXB2,numHUM,numIntegrationEvents > results.csv

for file in ./out*R1*/tmp/alns/*hum.bowtie.sam ; do
    seqRun1="${file%/*/*/*}"
    seqRun2="${seqRun1#*outYHo}"
    seqRun="${seqRun2%_R1}"
    sampleBase=$(basename ${file})
    sample="${sampleBase%.*.*.*}"
    sampleTrim="${sample%_*}"
    numHIVR1="$(samtools view outYHo${seqRun}_R1/tmp/alns/${sample}.hiv.bowtie.sam | wc -l)"
    numHIVR2="$(samtools view outYHo${seqRun}_R2/tmp/alns/${sample}.hiv.bowtie.sam | wc -l)"
    numHIVR1HX="$(samtools view outYHo${seqRun}_R1/tmp/alns/${sample}.hiv.hxb2.bowtie.sam | wc -l)"
    numHIVR2HX="$(samtools view outYHo${seqRun}_R2/tmp/alns/${sample}.hiv.hxb2.bowtie.sam | wc -l)"
    numHUMR1="$(samtools view outYHo${seqRun}_R1/tmp/alns/${sample}.hum.bowtie.sam | wc -l)"
    numHUMR2="$(samtools view outYHo${seqRun}_R2/tmp/alns/${sample}.hum.bowtie.sam | wc -l)"
    numTotalR1="$(tail -n 3 outYHo${seqRun}_R1/tmp/trim/${sampleTrim}_R1_001.fastq.gz_trimming_report.txt | head -n 1 - | awk '{print $1}')"
    numTotalR2="$(tail -n 3 outYHo${seqRun}_R2/tmp/trim/${sampleTrim}_R2_001.fastq.gz_trimming_report.txt | head -n 1 - | awk '{print $1}')"
    numPreprocessedR1="$(tail -n 2 outYHo${seqRun}_R1/tmp/trim/${sampleTrim}_R1_001.fastq.gz_trimming_report.txt | head -n 1 - | awk '{print $14}')"
    numPreprocessedR2="$(tail -n 2 outYHo${seqRun}_R2/tmp/trim/${sampleTrim}_R2_001.fastq.gz_trimming_report.txt | head -n 1 - | awk '{print $14}')"
    numIntegrationsR1=0
    numIntegrationsR2=0
    if [ ! -f ./outYHo${seqRun}_R1/${sample}_Pos.clean.csv  ]; then
        numIntegrationsR1=0
        echo ${seqRun},R1,${sample},${numTotalR1},${numPreprocessedR1},${numHIVR1},${numHIVR1HX},${numHUMR1},$((numIntegrationsR1)) >> ./results.csv
    else
        numIntegrationsR1="$(wc -l ./outYHo${seqRun}_R1/${sample}_Pos.clean.csv | awk '{print $1}')"
        echo ${seqRun},R1,${sample},${numTotalR1},${numPreprocessedR1},${numHIVR1},${numHIVR1HX},${numHUMR1},$((numIntegrationsR1-1)) >> ./results.csv
    fi
    if [ ! -f ./outYHo${seqRun}_R2/${sample}_Pos.clean.csv  ]; then
        numIntegrationsR2=0
        echo ${seqRun},R2,${sample},${numTotalR2},${numPreprocessedR2},${numHIVR2},${numHIVR2HX},${numHUMR2},$((numIntegrationsR2)) >> ./results.csv
    else
        numIntegrationsR2="$(wc -l ./outYHo${seqRun}_R2/${sample}_Pos.clean.csv | awk '{print $1}')"
        echo ${seqRun},R2,${sample},${numTotalR2},${numPreprocessedR2},${numHIVR2},${numHIVR2HX},${numHUMR2},$((numIntegrationsR2-1)) >> ./results.csv
    fi
done

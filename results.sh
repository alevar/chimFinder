#!/usr/bin/env bash

# The following format shall be used for the results page:
# 1.  - sampleName              - name of the sample
# 2.  - krakenHIV%              - % reads classified as Human immunodeficiency virus type 1 by Kraken
# 3.  - krakenHUM%              - % reads classified as Homo sapience by Kraken
# 4.  - chimericKrakenExtracted - number of unique chimeric reads extracted based on the kraken output
# 5.  - chimericKrakenFiltered  - number of unique chimeric reads after aligning with bowtie2 and filtering sam output with a python script
# 6.  - bowtie2HIV              - number of unique reads aligned to the HIV89.6 reference using bowtie2
# 7.  - bowtie2HUM              - number of unique reads aligned to the hg38 using bowtie2
# 8.  - chimericBowtie          - number of unique chimeric reads detected using bowtie2 alignments of all reads against both hg38 and hiv89.6 references
# 9.  - numSplits               - number of unique split points calculated from chimericKraken and chimericBowtie2 reads
# 10. - numReads                - number of unique chimericKraken and chimericBowtie2 reads which support potential split locations
# 11. - numSpliceJunctionsHIV   - number of unique splice junctions from Hisat2 output
# 12. - totalNumberReads        - total number of raw reads in the sample

# what needs to be done
# if a file such as *Pos.csv is absent - record

outDir=$1
inDir=$2

touch ./results.csv
echo "sample,krakenHiv%,krakenHum%,chimericReadsExtracted,chimericReadsFiltered,bowtie2HIV,bowtie2HUM,chimericBowtie,numSplits,numReads,numSpliceJunctionsHIV,totalNumberReads" > ./results.csv

for file in ${outDir}/krakenOut/*.report ; do
    sample=""
    HIVP=""
    HUMP=""
    nReadsKraken="0"
    postFiltKraken="0"
    bowtieReadsHIV="0"
    bowtieReadsHUM="0"
    chimericBowtie="0"
    numSplits="0"
    numReads="0"
    numSpliceJunctionsHIV="0"

    sampleBase=$(basename ${file})
    sample="${sampleBase%.*}"

    echo "totalNumberReads"
    # totalNumberReads="$(zcat ${inDir}/${sample}_R1_001.fastq.gz | echo $((`wc -l`/4)))"
    totalNumberReads=""

    HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
    HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
    echo "bowtieReadsHUM"
    bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.full.hum.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"

    if [ -f ${outDir}/${sample}.chim.csv ] && [ "$(wc -l ${outDir}/${sample}.chim.csv | awk -F ' ' '{print $1}')" > 0 ]; then
        echo "1.nReadsKraken"
        nReadsKraken="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
        echo "1.postFiltKraken"
        postFiltKraken="$(tail -n "$(wc -l ${outDir}/${sample}.chim.csv | awk -F ' ' '{print $1-1}')" ${outDir}/${sample}.chim.csv | awk -F "," '{print $2}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
        echo "1.numSplits"
        numSplits="$(wc -l ${outDir}/${sample}_Pos.csv | awk -F ' ' '{print $1}')"
        echo "1.numReads"
        numReads="$(awk -F ',' '{sum += $5} END {print sum}' ${outDir}/${sample}_Pos.csv)"
        echo "${numReads}"
        if [ -f ${outDir}/fullAlignments/${sample}.full.hiv.csv ] && [ "$(wc -l ${outDir}/fullAlignments/${sample}.full.hiv.csv | awk -F ' ' '{print $1}')" > 0 ]; then
            echo "1.1.bowtieHIV"
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.full.hiv.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            if [ -f ${outDir}/${sample}.full.csv ] && [ "$(wc -l ${outDir}/${sample}.full.csv | awk -F ' ' '{print $1}')" > 0 ]; then
                echo "1.1.1.chimericBowtie"
                chimericBowtie="$(awk -F ',' '$4==0.0' ${outDir}/${sample}.full.csv | wc -l | awk -F ' ' '{print $1}')"
                echo "1.1.1.numSpliceJunctionsHIV"
                if [ -f ${outDir}/hisat/${sample}.junctions ]; then 
                    numSpliceJunctionsHIV="$(wc -l ${outDir}/hisat/${sample}.junctions | awk -F ' ' '{print $1}')"
                else
                    numSpliceJunctionsHIV="0"
                fi
            else
                echo "1.1.2chimericBowtie=0"
                chimericBowtie="0"
                echo "1.1.2numSpliceJunctionsHIV=0"
                numSpliceJunctionsHIV="0"
            fi
        else
            echo "1.2bowtieReadsHIV=0"
            bowtieReadsHIV="0"
            echo "1.2chimericBowtie=0"
            chimericBowtie="0"
            echo "1.2numSpliceJunctionsHIV=0"
            numSpliceJunctionsHIV="0"
        fi
    else
        echo "2.nReadsKraken=0"
        nReadsKraken="0"
        echo "2.postFiltKraken=0"
        postFiltKraken="0"
        if [ -f ${outDir}/fullAlignments/${sample}.full.hiv.csv ] && [ "$(wc -l ${outDir}/fullAlignments/${sample}.full.hiv.csv | awk -F ' ' '{print $1}')" > 0 ]; then
            echo "2.1.bowtieReadsHIV"
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.full.hiv.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            if [ -f ${outDir}/${sample}.full.csv ] && [ "$(wc -l {outDir}/${sample}.full.csv | awk -F ' ' '{print $1}')" > 0 ]; then
                echo "2.1.1.chimericBowtie"
                chimericBowtie="$(awk -F ',' '$4==0.0' ${outDir}/${sample}.full.csv | wc -l | awk -F ' ' '{print $1}')"
                echo "2.1.1.numSpliceJunctionsHIV"
                if [ -f ${outDir}/hisat/${sample}.junctions ]; then 
                    numSpliceJunctionsHIV="$(wc -l ${outDir}/hisat/${sample}.junctions | awk -F ' ' '{print $1}')"
                else
                    numSpliceJunctionsHIV="0"
                fi
                echo "2.1.1.numSplits"
                numSplits="$(wc -l ${outDir}/${sample}_Pos.csv | awk -F ' ' '{print $1}')"
                echo "2.1.1.numreads"
                numReads="$(awk -F ',' '{sum += $5} END {print sum}' ${outDir}/${sample}_Pos.csv)"
            else
                echo "2.1.2.chimericBowtie=0"
                chimericBowtie="0"
                echo "2.1.2.numSpliceJunctionsHIV=0"
                numSpliceJunctionsHIV="0"
                echo "2.1.2.numSplits=0"
                numSplits="0"
                echo "2.1.2.numreads=0"
                numReads="0"
            fi
        else
            echo "2.2.bowtieReadsHIV=0"
            bowtieReadsHIV="0"
            echo "2.2.chimericBowtie=0"
            chimericBowtie="0"
            echo "2.2.numSpliceJunctionsHIV=0"
            numSpliceJunctionsHIV="0"
            echo "2.2.numSplits=0"
            numSplits="0"
            echo "2.2.numreads=0"
            numReads="0"
        fi
    fi
    echo "${sample},${HIVP},${HUMP},${nReadsKraken},${postFiltKraken},${bowtieReadsHIV},${bowtieReadsHUM},${chimericBowtie},${numSplits},${numReads},${numSpliceJunctionsHIV},${totalNumberReads}" >> ./results.csv
done
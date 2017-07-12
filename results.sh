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
# 11. - numSpliceJunctionsHIV      - number of unique splice junctions from Hisat2 output

outDir=$1

touch ./results.csv
echo "sample,krakenHiv%,krakenHum%,chimericReadsExtracted,chimericReadsFiltered,bowtie2HIV,bowtie2HUM,chimericBowtie,numSplits,numReads,numSpliceJunctionsHIV" > ./results.csv

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

    if [ ! -f ${outDir}/localAlignments/${sample}.chim.csv ]; then # if chimeric kraken reads do not exist
        if [ ! -f ${outDir}/fullAlignments/${sample}.hiv.full.sam ]; then # if full hiv alignment is empty
            HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
            HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
            nReadsKraken="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
            bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            numSplits="$()"
            numReads="$()"
            numSpliceJunctionsHIV="$()"
        else # if full hiv alignment is not empty
            HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
            HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
            nReadsKraken="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.hiv.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            chimericBowtie="$(awk -F ',' '$4==0.0' Y430_pos_12_S51_Pos.csv | wc -l)"
            numSplits="$(wc -l Y430_pos_12_S51_Pos.csv | awk -F ' ' '{print $1}')"
            numReads="$(awk -F ',' '{sum += $3} END {print sum}' Y430_pos_12_S51_Pos.csv)"
            numSpliceJunctionsHIV="$(wc -l hisat/...junctions)"
        fi
    else # if local chimeric kraken reads do exist
        if [ ! -f ${outDir}/fullAlignments/${sample}.hiv.full.sam ]; then # if full hiv alignment is empty
            HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
            HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
            nReadsKraken="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
            postFiltKraken="$(tail -n "$(wc -l ${outDir}/localAlignments/${sample}.chim.csv | awk -F ' ' '{print $1-1}')" ${outDir}/localAlignments/${sample}.chim.csv | awk -F "," '{print $2}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            numSplits="$(wc -l Y430_pos_12_S51_Pos.csv | awk -F ' ' '{print $1}')"
            numReads="$(awk -F ',' '{sum += $3} END {print sum}' Y430_pos_12_S51_Pos.csv)"
            numSpliceJunctionsHIV="$(wc -l hisat/...junctions)"
        else # if full hiv alignment is not empty
            HIVP="$(grep 'Human immunodeficiency' ${file} | awk -F '\t' '{print $1}')"
            HUMP="$(grep 'Homo sapiens' ${file} | awk -F '\t' '{print $1}')"
            nReadsKraken="$(wc -l ${outDir}/krakenOut/selected/${sample}.chim | awk -F ' ' '{print $1}')"
            postFiltKraken="$(tail -n "$(wc -l ${outDir}/localAlignments/${sample}.chim.csv | awk -F ' ' '{print $1-1}')" ${outDir}/localAlignments/${sample}.chim.csv | awk -F "," '{print $2}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            bowtieReadsHUM="$(samtools view ${outDir}/fullAlignments/${sample}.hum.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            numSplits="$(wc -l Y430_pos_12_S51_Pos.csv | awk -F ' ' '{print $1}')"
            numReads="$(awk -F ',' '{sum += $3} END {print sum}' Y430_pos_12_S51_Pos.csv)"
            numSpliceJunctionsHIV="$(wc -l hisat/...junctions)"
            bowtieReadsHIV="$(samtools view ${outDir}/fullAlignments/${sample}.hiv.full.sam | awk -F '\t' '{print $1}' | sort -u | wc -l | awk -F ' ' '{print $1}')"
            chimericBowtie="$(awk -F ',' '$4==0.0' Y430_pos_12_S51_Pos.csv | wc -l)"
        fi
    fi
    echo "${sample},${HIVP},${HUMP},${nReadsKraken},${postFiltKraken},${bowtieReadsHIV},${bowtieReadsHUM},${chimericBowtie},${numSplits},${numReads},${numSpliceJunctionsHIV}" >> ./results.csv
done

#!/usr/bin/env bash

inputDir=${1}
annotation=${2}

for file in ${inputDir}/*/*R1*fastq.gz ; do
	sampleR1Base=$(basename ${file})
	sampleR1="${sampleR1Base%.*}"
	sample="${sampleR1%_R1*}"
	sampleR2Base="${sampleR1%_R1*}"_R2.fastq
    baseEnd="${sampleR1##*_R1}"
    # echo ${sample}_R1${baseEnd}

    dirP=${file%/*}
    outputDir=outSample${dirP#*Sample}
    # outputDir=outAdditional
    echo ${dirP} 

    echo ANALYZING ${outputDir}_1/tmp/alns/${sample}.hiv.hxb2.bowtie.sam
    # ./chimFinder.py --input1r1 ${outputDir}_1/tmp/alns/${sample}.hiv.bowtie.sam \
                    # --input2r1 ${outputDir}_1/tmp/alns/${sample}.hum.bowtie.sam \
                    # --input1r2 ${outputDir}_2/tmp/alns/${sample}.hiv.bowtie.sam \
                    # --input2r2 ${outputDir}_2/tmp/alns/${sample}.hum.bowtie.sam \
                    # --splicedR1 ${outputDir}_1/tmp/alns/${sample}.hum.hisat.sam \
                    # --splicedR2 ${outputDir}_2/tmp/alns/${sample}.hum.hisat.sam \
                    # -o ${outputDir}/${sample}_2550_mc1 \
                    # -t 12 \
                    # --minLen 30 \
                    # --maxLenUnmapped 5 \
                    # --minCount 1 \
                    # --maxCountPenalty 1 \
                    # -a ${annotation} \
                    # --overlap 5 \
                    # --gap 5 \
                    # --minEntropy 0.84 \
                    # --close 30 \
                    # --score 0.8 \
                    # -q
    # samtools mpileup -uf /ccb/salz3/avaraby/genomes/HIV1/HXB2/HXB2.fa ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq | seqtk seq -a - > ${outputDir}/consensus/${sample}.hiv.hxb2.SEP.fa
    samtools view -h ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam > ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sam
    sam2splice.py -i ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sam -o ${outputDir}/${sample}.sj.csv --fasta
    sam2splice.py -i ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sam -o ${outputDir}/${sample}.sj.f.csv -f
done

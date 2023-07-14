#!/usr/bin/env bash

inputDir=${1}
outputDir=${2}
genome1DB=${3}
humanDB=${4}
genome2_annotation=${5}
threads=${6}
humanSplicing=${7}
trim_galore=${8}

mkdir -p ${outputDir}_R1
mkdir -p ${outputDir}_R1/tmp
mkdir -p ${outputDir}_R1/tmp/trim
mkdir -p ${outputDir}_R1/tmp/alns
mkdir -p ${outputDir}_R1/tmp/fq

mkdir -p ${outputDir}_R2
mkdir -p ${outputDir}_R2/tmp
mkdir -p ${outputDir}_R2/tmp/trim
mkdir -p ${outputDir}_R2/tmp/alns
mkdir -p ${outputDir}_R2/tmp/fq

mkdir -p ${outputDir}
mkdir -p ${outputDir}/tmp
mkdir -p ${outputDir}/alns
mkdir -p ${outputDir}/consensus

for file in ${inputDir}/*R1*fastq.gz ; do
	sampleR1Base=$(basename ${file})
	sampleR1="${sampleR1Base%.*.*}"
	sample="${sampleR1%_R1*}"
	sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz
	baseEnd="${sampleR1##*_R1}"

	echo ANALYZING ${sample}
	echo BEGIN DATE: ${date}
	TOTAL_TIME=0

	echo ADAPTOR AND QUALITY TRIMMING R1
	${trim_galore} --trim-n \
					-q 5 \
					--phred33 \
					--length 45 \
					-o ./${outputDir}_R1/tmp/trim/ \
					--dont_gzip ${inputDir}/${sampleR1Base}
	trimmedFQ=./${outputDir}_R1/tmp/trim/"${sampleR1%_R1*}"_R1"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	echo ADAPTOR AND QUALITY TRIMMING R2
	${trim_galore} --trim-n \
					-q 5 -\
					-phred33 \
					--length 45 \
					-o ./${outputDir}_R2/tmp/trim/ \
					--dont_gzip ${inputDir}/${sampleR2Base}
	trimmedFQ=./${outputDir}_R2/tmp/trim/"${sampleR1%_R1*}"_R2"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R1 TO first genome index
	bowtie2 --very-sensitive-local \
			--no-unal \
			--local \
			--phred33 \
			-p ${threads} \
			-x ${genome1DB} \
			-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq \
			-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i1.bowtie.sam \
			--al ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R2 TO first genome index
	bowtie2 --very-sensitive-local \
			--no-unal \
			--local \
			--phred33 \
			-p ${threads} \
			-x ${genome1DB} \
			-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq \
			-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i1.bowtie.sam \
			--al ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
	
	# extract read names into separate file for r1
	sed -n '1~4p' ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq > ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.reads
	# extract read names into separate file for r2
	sed -n '1~4p' ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq > ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.reads

	# get list of uniq names
	cat ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.reads ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.reads | awk -F ' ' '{print $1}' | sort | uniq > ${outputDir}/${sample}${baseEnd}.reads

	# now use this list to extract data for the human alignments
	grep -A 3 --no-group-separator -w -Ff ${outputDir}/${sample}${baseEnd}.reads ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq > ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq
	grep -A 3 --no-group-separator -w -Ff ${outputDir}/${sample}${baseEnd}.reads ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq > ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq


	SECONDS=0
	echo BOWTIE2 ALIGNING R1 TO the second genome index
	bowtie2 --very-sensitive-local \
			--no-unal \
			--local \
			--phred33 \
			-p ${threads} \
			-x ${genome2DB} \
			-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq \
			-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i2.bowtie.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R1 TO the second genome index
	hisat2 -p ${threads} \
			--very-sensitive \
			--end-to-end \
			--known-splicesite-infile ${humanSplicing} \
			--no-unal \
			--phred33 \
			-x ${genome2DB} \
			-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq \
			-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i2.hisat.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R2 TO the second genome index
	bowtie2 --very-sensitive-local \
			--no-unal \
			--local \
			--phred33 \
			-p ${threads} \
			-x ${genome2DB} \
			-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq \
			-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i2.bowtie.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R2 TO the second genome index
	hisat2 -p ${threads} \
			--very-sensitive \
			--end-to-end \
			--known-splicesite-infile ${humanSplicing} \
			--no-unal \
			--phred33 \
			-x ${genome2DB} \
			-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq \
			-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i2.hisat.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo ANALYZING
	./chimFinder.py --input1r1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i1.bowtie.sam \
					--input2r1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i2.bowtie.sam \
					--input1r2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i1.bowtie.sam \
					--input2r2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i2.bowtie.sam \
					--splicedR1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.i2.hisat.sam \
					--splicedR2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.i2.hisat.sam \
					-o ${outputDir}/${sample}${baseEnd} \
					-t 12 \
					--minLen 20 \
					--maxLenUnmapped 30 \
					-a ${genome2_annotation} \
					--overlap 5 \
					--gap 5 \
					--minEntropy 0.84 \
					--close 5 \
					--score 0.75 \
					-q
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
	SECONDS=0
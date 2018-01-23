#!/usr/bin/env bash

inputDir=${1}
outputDir=${2}
hivDB=${3}
humanDB=${4}
annotation=${5}
threads=${6}
humanSplicing=${7}
hivSplicing=${8}
HXB2=${9}
HXB2fa=${10}
trim_galore=${11}

echo ${HXB2fa}

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

	# R1=============================================R1
	echo ADAPTOR AND QUALITY TRIMMING R1
	${trim_galore} --trim-n -q 5 --phred33 --length 45 -o ./${outputDir}_R1/tmp/trim/ --dont_gzip ${inputDir}/${sampleR1Base}
	trimmedFQ=./${outputDir}_R1/tmp/trim/"${sampleR1%_R1*}"_R1"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R1 TO HIV
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam --al ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R1 TO HG38
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R1 TO HG38
	hisat2 -p ${threads} --very-sensitive --end-to-end --known-splicesite-infile ${humanSplicing} --no-unal --phred33 -x ${humanDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R1 TO HXB2
	hisat2 -p ${threads} --very-sensitive --end-to-end --novel-splicesite-outfile ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.novelSites --known-splicesite-infile ${hivSplicing} --no-unal --phred33 -x ${HXB2} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam --un ${outputDir}_R1/tmp/fq/${sample}${baseEnd}.un.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R1 TO HXB2
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${HXB2} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}.un.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
	rm ${outputDir}_R1/tmp/fq/${sample}${baseEnd}.un.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo INDEXING R1
	samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.bam -
	samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam -
	samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo ANALYZING R1
	./chimFinder.py -i1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam -i2 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam --spliced ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam -o ${outputDir}_R1/${sample}${baseEnd} -t 12 --minLen 20 -a ${annotation} --overlap 5 --gap 5 --minEntropy 0.84 --close 5 --score 0.75 -q #> ${outputDir}_R1/${sample}${baseEnd}.report
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
	SECONDS=0

	#R2=============================================R2
	echo ADAPTOR AND QUALITY TRIMMING R2
	${trim_galore} --trim-n -q 5 --phred33 --length 45 -o ./${outputDir}_R2/tmp/trim/ --dont_gzip ${inputDir}/${sampleR2Base}
	trimmedFQ=./${outputDir}_R2/tmp/trim/"${sampleR1%_R1*}"_R2"${sampleR1##*_R1}"_trimmed.fq
	mv ${trimmedFQ} ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R2 TO HIV
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam --al ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R2 TO HG38
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R2 TO HG38
	hisat2 -p ${threads} --very-sensitive --end-to-end --known-splicesite-infile ${humanSplicing} --no-unal --phred33 -x ${humanDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo HISAT ALIGNING R2 TO HXB2
	hisat2 -p ${threads} --very-sensitive --end-to-end --novel-splicesite-outfile ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.novelSites --known-splicesite-infile ${hivSplicing} --no-unal --phred33 -x ${HXB2} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam --un ${outputDir}_R2/tmp/fq/${sample}${baseEnd}.un.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BOWTIE2 ALIGNING R2 TO HXB2
	bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${HXB2} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}.un.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
	rm ${outputDir}_R2/tmp/fq/${sample}${baseEnd}.un.fastq
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo INDEXING R2
	samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.bam -
	samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam -
	samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam

	samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam -
	samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo ANALYZING R2
	./chimFinder.py -i1 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam -i2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam --spliced ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam -o ${outputDir}_R2/${sample}${baseEnd} -t 12 --minLen 20 -a ${annotation} --overlap 5 --gap 5 --minEntropy 0.84 --close 5 --score 0.75 -q #> ${outputDir}_R2/${sample}${baseEnd}.report
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
	SECONDS=0

	#LANDSCAPES=============================================LANDSCAPES
	SECONDS=0
	echo MERGING SEP
	samtools merge -f ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam
	samtools sort -o ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sorted.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam
	mv ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sorted.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam
	sleep 10s
	samtools index ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam.bai
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	SECONDS=0
	echo BUILDING CONSENSUS SEQUENCE FROM SEP
	samtools mpileup -uf ${HXB2fa} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq | seqtk seq -a - > ${outputDir}/consensus/${sample}${baseEnd}.hiv.hxb2.SEP.fa
	DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
	echo DONE IN ${DUR}
	TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

	T_DUR="$(($TOTAL_TIME / 60)) minutes and $(($TOTAL_TIME % 60)) seconds"
	echo TOTAL TIME ELAPSED: ${T_DUR}
	echo "$(date)"
done
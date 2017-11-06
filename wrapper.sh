#!/usr/bin/env bash

inputDir=$1
outputDir=$2
hivDB=$3
humanDB=$4
annotation=$5
threads=$6
humanSplicing=$7
hivSplicing=$8
HXB2=$9
HXB2fa=${10}

echo ${HXB2fa}

mkdir -p ${outputDir}_R1
mkdir -p ${outputDir}_R1/tmp
mkdir -p ${outputDir}_R1/tmp/alns
mkdir -p ${outputDir}_R1/tmp/fq

mkdir -p ${outputDir}_R2
mkdir -p ${outputDir}_R2/tmp
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

	if [[ $sample == *"pos"* || $sample == *"Pos"* ]]; then
		echo ANALYZING ${sample}
		echo BEGIN DATE: ${date}
		TOTAL_TIME=0

		#R1=============================================R1
		
		# echo DECOMPRESSING R1
		# SECONDS=0
		# zcat ${inputDir}/${sampleR1Base} > ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R1 TO HG38
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R1 TO HIV
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R1 TO HXB2
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${HXB2} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo HISAT ALIGNING R1 TO HG38
		# hisat2 -p ${threads} --very-sensitive --end-to-end --known-splicesite-infile ${humanSplicing} --no-unal --phred33 -x ${humanDB} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo HISAT ALIGNING R1 TO HXB2
		# hisat2 -p ${threads} --very-sensitive --end-to-end --novel-splicesite-outfile ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.novelSites --known-splicesite-infile ${hivSplicing} --no-unal --phred33 -x ${HXB2} -U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo FIXING HXB2 REFERENCE NAMES IN R1
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo INDEXING R1
		# samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.bam -
		# samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam | samtools sort -o ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam -
		# samtools index -@ ${threads} ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo ANALYZING R1
		./hiv.py -i1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam -i2 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam --spliced ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hum.hisat.sam -o ${outputDir}_R1/${sample}${baseEnd} -t 12 --minLen 30 -a ${annotation} --overlap 5 --gap 5
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
		SECONDS=0


		#R2=============================================R2


		# echo DECOMPRESSING R2
		# SECONDS=0
		# zcat ${inputDir}/${sampleR2Base} > ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R2 TO HG38
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${humanDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R2 TO HIV
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${hivDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BOWTIE2 ALIGNING R2 TO HXB2
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${HXB2} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo HISAT ALIGNING R2 TO HG38
		# hisat2 -p ${threads} --very-sensitive --end-to-end --known-splicesite-infile ${humanSplicing} --no-unal --phred33 -x ${humanDB} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo HISAT ALIGNING R2 TO HXB2
		# hisat2 -p ${threads} --very-sensitive --end-to-end --novel-splicesite-outfile ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.novelSites --known-splicesite-infile ${hivSplicing} --no-unal --phred33 -x ${HXB2} -U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo FIXING HXB2 REFERENCE NAMES IN R2
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo INDEXING R2
		# samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.bam -
		# samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.sam | samtools sort -o ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam -
		# samtools index -@ ${threads} ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo ANALYZING R2
		./hiv.py -i1 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.sam -i2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.bowtie.sam --spliced ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hum.hisat.sam -o ${outputDir}_R2/${sample}${baseEnd} -t 12 --minLen 30 -a ${annotation} --overlap 5 --gap 5
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
		SECONDS=0


		#PAIRED=============================================PAIRED

		# SECONDS=0
		# echo BOWTIE2 ALIGNING PAIRED TO HXB2
		# bowtie2 --very-sensitive-local --no-unal --local --phred33 -p ${threads} -x ${HXB2} -1 ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -2 ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo HISAT ALIGNING PAIRED TO HXB2
		# hisat2 -p ${threads} --very-sensitive --novel-splicesite-outfile ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.novelSites --known-splicesite-infile ${hivSplicing} --no-unal --phred33 -x ${HXB2} -1 ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq -2 ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq -S ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo FIXING HXB2 REFERENCE NAMES PAIRED
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.sam
		# sed -i 's/K03455|HIVHXB2CG/HIVHXB2CG/g' ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.sam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo INDEXING PAIRED
		# samtools view -S -@ ${threads} -b ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.sam | samtools sort -o ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.bam -
		# samtools index -@ ${threads} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.bam

		# samtools view -S -@ ${threads} -b ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.sam | samtools sort -o ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.bam -
		# samtools index -@ ${threads} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.bam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo MERGING PAIRED
		# samtools merge ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.bowtie.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.paired.hisat.bam
		# samtools sort -o ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam.sorted ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam
		# mv ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam.sorted ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam
		# samtools index -@ ${threads} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo FILTERING SEP ALIGNMENTS
		# samtools view ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | awk -F '\t' '$6 ~ /N/ {print $1}' > ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.splicedReads.txt
		# samtools view ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | awk -F '\t' '$6 ~ /N/ {print $1}' > ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.splicedReads.txt
		
		# samtools view -h ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | awk -F '\t' '$6 ~ /N/ {print}' > ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.sam
		# samtools view -h ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | awk -F '\t' '$6 ~ /N/ {print}' > ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.sam

		# samtools view -H ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | cat - ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.sam | samtools view -h -bS -o ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.bam -
		# samtools view -H ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.hisat.bam | cat - ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.sam | samtools view -h -bS -o ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.bam -
		
		# samtools view -h ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam | grep -vf ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.splicedReads.txt | samtools view -bS -o ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.bowtie.bam -
		# samtools view -h ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.hiv.hxb2.bowtie.bam | grep -vf ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.splicedReads.txt | samtools view -bS -o ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.bowtie.bam -
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo MERGING SEP
		# samtools merge ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.bowtie.bam ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.bam ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.bowtie.bam ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.bam
		# samtools sort -o ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sorted.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam
		# mv ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.sorted.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam
		# sleep 10s
		# samtools index ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam.bai
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# rm ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.splicedReads.txt
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.splicedReads.txt
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.bam
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.bam
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.bowtie.bam
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.bowtie.bam
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R1.hiv.hxb2.hisat.sam
		# rm ${outputDir}/tmp/${sample}${baseEnd}_R2.hiv.hxb2.hisat.sam

		# SECONDS=0
		# echo BUILDING CONSENSUS SEQUENCE FROM PAIRED
		# samtools mpileup -uf ${HXB2fa} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.PAIRED.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > ${outputDir}/consensus/${sample}${baseEnd}.hiv.hxb2.PAIRED.fq
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# SECONDS=0
		# echo BUILDING CONSENSUS SEQUENCE FROM SEP
		# samtools mpileup -uf ${HXB2fa} ${outputDir}/alns/${sample}${baseEnd}.hiv.hxb2.SEP.bam | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq | seqtk seq -a - > ${outputDir}/consensus/${sample}${baseEnd}.hiv.hxb2.SEP.fa
		# DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		# echo DONE IN ${DUR}
		# TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		# T_DUR="$(($TOTAL_TIME / 60)) minutes and $(($TOTAL_TIME % 60)) seconds"
		# echo TOTAL TIME ELAPSED: ${T_DUR}
		# echo "$(date)"
	fi
done
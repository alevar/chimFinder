#!/usr/bin/env bash

HIV_ref=$1 #file
HUM_ref=$2 #pth to wherer chromosome references are
outDir=$3

#need to get basename

HIV_ref_base=$(basename ${HIV_ref})
HUM_ref_base=$(basename ${HUM_ref})
HIV_ref_base="${HIV_ref_base%.*}"
HUM_ref_base="${HUM_ref_base%.*}"

mkdir -p ${outDir}
mkdir -p ${outDir}/refs

#dust the references and build indices

dustmasker -in ${HIV_ref} -infmt fasta -out ${outDir}/refs/${HIV_ref_base}.dusted.fa -outfmt fasta
sed -e '/^>/! s/[[:lower:]]/N/g' ${outDir}/refs/${HIV_ref_base}.dusted.fa > ${outDir}/refs/${HIV_ref_base}.fa
bowtie2-build -f ${outDir}/refs/${HIV_ref_base}.fa ${HIV_ref_base}
hisat2-build -p 12 ${outDir}/refs/${HIV_ref_base}.fa ${outDir}/refs/${HIV_ref_base}

kraken-build --add-to-library ${outDir}/refs/${HIV_ref_base}.fa --db ${outDir}/refs/customKrakenDB #add dusted hiv reference to krakenDB

#for human refs need to iterate over the chromosomes

for file in HUM_ref/*fa;do
	base=$(basename ${file})
	base="${base%.*}"

	dustmasker -in ${file} -infmt fasta -out ${outDir}/refs/human/${base}.dusted.fa -outfmt fasta
	sed -e '/^>/! s/[[:lower:]]/N/g' ${outDir}/refs/human/${base}.dusted.fa > ${outDir}/refs/human/${base}.fa
	rm ${outDir}/refs/human/${base}.dusted.fa
	# add dusted human reference to the krakenDB
	kraken-build --add-to-library ${outDir}/refs/${base}.fa --db ${outDir}/refs/customKrakenDB
done

humList= #comma-separated list of human reference chromosomes to be used when building a reference for bowtie2 and hisat2

bowtie2-build -f ${outDir}/refs/${HUM_ref_base}.fa ${HUM_ref_base}
hisat2-build -p 12 ${outDir}/refs/${HUM_ref_base}.fa ${outDir}/refs/${HUM_ref_base}
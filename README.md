## Motivation

chimFinder is designed to identify reads from RNA and DNA sequencing experiments which align to two distinct genomes. This tool was developed as part of the work on identification of fusion transcripts resulting from HIV type 1 integrations in human CD4+ T cells.

## Prerequisites
Along with chimFinder this repository contains several scripts designed to automate all steps of the pipeline described in the manuscript:
1. database cleanup (contamination.py)
2. read trimming and preprocessing
3. alignment to the genomes
4. full spliced landscape generation of the integrated genome
5. identification integration sites
6. final report

Several prerequisites are assumed to be installed and properly configured to be able to run the full pipeline
1. unix environment 
2. bowtie 
3. hisat
4. samtools
5. bowtie2 databases
6. hisat2 databases
7. trim galore (including fastq and cutadapt)
8. ncbi-toolkit

chimFinder uses as input two sam-formatted files. Users may choose to perform upstream analysis using other tools and pipelines which generate alignments in the appropriate sam format.

## Installation
chimFinder is implemented in Python and assumes minimum requirements.
1. git clone
2. pip install -r requirements.txt
(Optional) - wrapper.sh may be used to run the entire pipeline in an automated fashion. Please ensure all software is installed and properly configured on the machine
3. ./chimFinder.py -i1 <hiv.bowtie.sam> -i2 <hum.bowtie.sam> --spliced <hum.hisat.sam> -o <> -a <> <other parameters>

## Parameters
-h, --help, help page

-i1r1 INPUT1R1, --input1r1 INPUT1R1
                    first alignment r1

-i1r2 INPUT1R2, --input1r2 INPUT1R2
                    first alignment r2

-i2r1 INPUT2R1, --input2r1 INPUT2R1
                    second alignment r1

-i2r2 INPUT2R2, --input2r2 INPUT2R2
                    second alignment r2

-o OUT, --out OUT     
					output file

-t THREADS, --threads THREADS
                    the number of threads to use in the computation

--minCount MINCOUNT   the minimum number of reads supporting an integration
                    site.

--maxCountPenalty MAXCOUNTPENALTY
                    the maximum penalty to give for minCount when hit.

--weightCount WEIGHTCOUNT
                    Weight of the number of reads score in the cumulative
                    score equation

--minEntropy MINENTROPY
                    minimum alignment entropy score.

--weightEntropy WEIGHTENTROPY
                    Weight of the entropy score in the cumulative score
                    equation

--minQual MINQUAL     minimum mean mapping quality of sequenced read.
                    Everything below this threshold will be reported as
                    multialignment for which no assumption can be made
                    from the annotation

--weightQual WEIGHTQUAL
                    Weight of the mapping score in the cumulative score
                    equation

--minLen MINLEN       the minimum number of nucleotides in alignment to keep
                    a read.

--maxLenUnmapped MAXLENUNMAPPED
                    minimum percent of the read to be aligned correctly in
                    spanning read search

--weightLen WEIGHTLEN
                    Weight of the alignment length score in the cumulative
                    score equation

--steepSlopeAL STEEPSLOPEAL
                    Weight of the alignment length score in the cumulative
                    score equation

--maxAlLenPenalty MAXALLENPENALTY
                    the maximum penalty to give for minLen when hit.

-s SCORE, --score SCORE
                    the minimum overall score to keep.

--overlap OVERLAP     overlap threshold

--gap GAP             gap threshold

--close CLOSE         distance between two integration sites to group
                    together

-a ANNOTATION, --annotation ANNOTATION
                    annotation for the human genome

-w, --writeReads      write reads to fasta files

-p, --plot            plot snapshots

-q, --quiet           do not print to std out. report will still be saved

--splicedR1 SPLICEDR1
                    spliced end-to-end alignment r1. Reads will be
                    subtracted from the main alignments and will not be
                    reported in the final report of integrations sites

--splicedR2 SPLICEDR2
                    spliced end-to-end alignment r2. Reads will be
                    subtracted from the main alignments and will not be
                    reported in the final report of integrations sites

## Output
chimFinder outputs the following three files per run:
1. .csv - standard ouput of integration sites grouped by position
2. .full.csv - same as .csv but contains additional information about the reads in each group, such as read name
3. .report - read statistics for each of the applied filters
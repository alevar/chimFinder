## Current Plan:
1. Run Kraken with a database for Human and HIV genomes
2. Extract all reads which might contain both HIV and Human segments in accordance with the Kraken output
3. Build fastq files containing all candidate paired-end reads
4. Align extracted paired-end reads using bowtie2 to both the human and hiv databases and output alignments in SAM format
5. Identify potential candidates for chimeric reads
6. For each candidate read find overlapping reads (for which the alignment spans portion of the alignment)
7. For each candidate read check that neighboring reads belong to the same alignment
8. for each candidate read find reads which extend the HIV section. Check how far the read ends

### Possible alterations:
1. do not run Kraken. Instead run bowtie2 against both human and hiv databases using a full set of reads.

## Forseeable Future:
1. Based on the qualities computed above score the confidence of the integration site
2. Build an automated pipeline

## Distant Future:
1. build an efficient standalone computational tool to replace this pipeline. Given two different references the tool will identify integration sites.

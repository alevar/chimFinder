##Motivation for complexity of the application:
Initially we tried Kraken and extracting chimeric reads from the kraken output
Then we tried adding full bowtie2 alignments to the pipeline and observed an increase in true positives
With an increase in true positives there was  an enormous increase in false positives observed (~40K groupings)
It would be next to impossible to find true positives manually among all this data hence we tried to do some additional preprocessing of the reads:
Removing duplicated reads using picard tools did drastically decrease the number of false positives but along introduced false negatives
Dusting low-complexity regions in the reference had the same effect of introducing false negatives
In order to prevent introduction of false negatives to the results and reduce the false positives, making it feasible to sort through the results manually we implemented the algorithm described in greater detail below. The competitive advantage of the algorithm comes from the fact that several parameters are computed from the alignment information for all possible chimeric reads before any filtering is implemented. Such approach allows to push all the top hits to the top based on the qualities typically used in deciding whether a read is chimeric (overlap length of unique alignment, complexity, mapping quality, sequencing quality etc.)

##Detailed description of the software

Load alignment data for hiv and human
Extract sam flag parameters and check if paired or unpaired. Downstream analysis will depend partially on pairedness.
If hisat2 spliced alignment is provided through the command line arguments:
Load
Extract sam flag parameters and discard any reads which are not spliced (do not contain “N” in the CIGAR string)
Create a set of all resulting reads
Subtract this set from the bowtie alignments
Extract reference and template start and end from the CIGAR string of each alignment
Perform initial filtering
Remove secondary alignments
Remove PCR duplicates as marked in alignment
Remove suppAl == 0
Remove noPassFilter == 0
Remove all end-to-end reads
If paired, R1 and R2 are allowed to remain if end-to-ed in distinct genomes
Combine records from two distinct alignments preserving all information
Mark hiv-human alignment orientation on each template read
Calculate the overlap between hiv and human alignments on the template read
Create records unique to each read
Each record contains information specific to its alignment conformation
For each such record the following information is extracted
Entropy of hiv and human alignments outside the overlapping region
Length of hiv and human alignments outside the overlapping region
Unique ID is assigned to each record from the HIV-human integration position. The structure is as follows:
End position of the first genome
:
Start position of the last genome
@
ID of the first genome
:
ID of the last genome
Collapse records based on the following information:
comb
HUM_ID
R
Orient
For each of the resulting groupings compute the following:
Count of all reads
List of all read names that went into the grouping
Minimum value for the start of the alignment on the human reference
Minimum value for the end of the alignment on the human reference
Minimum value for the start of the alignment on the HIV reference
Minimum value for the end of the alignment on the HIV reference
Sum of HIV alignment length
Sum of human alignment length
Sum of read lengths
Sum of entropy scores of the template region outside the alignment overlap region which aligned against human reference
Sum of entropy score of the template region outside the alignment overlap region which aligned against hiv reference
Sum of mapping qualities of HIV alignments for reads in the grouping as reported by the aligner
Sum of mapping qualities of human alignments for reads in the grouping as reported by the aligner
Annotate records by the closest human splice site to the human alignment of the read
Sort records hierarchically based on the human nearest splice site and the start of the human alignment
For each human splice site sorting compute distance between human alignments
Separate the resulting alignments into two groups based on the MAPQ property conserved from SAM output
lowMAPQ contains records below MAPQ threshold
highMAPQ contains records above or equal to the MAPQ threshold
Collapse records further based on identical nearest human splice site if the distance between alignments does not exceed the closeness parameter (30nt default)
Compute mean of each summed parameter based on how many reads went into group
Compute raw score in accordance with the equation
Apply the sigmoid power function based on the count score to produce the final scoring of records
Sort records by the score
Remove all records which do not pass thresholds specified in the parameters.

Notes and future additional improvements:
Splicing
As can be seen from several examples, UCSC genome browser indicates human splicing to be present in samples characterized as chimeric.
It would be interesting to run alignments with hisat instead to see if a spliced aligner helps resolve the issue
Sigmoid parameter transformation
So far only the count score was applied as a sigmoid function
However as can be seen from the outYHo29 KANSL3 record (indicated by Dr. Ho as a very likely true chimera) it is pushed down the list because the alignment length is very close to the threshold.
By applying the alignment length as a sigmoid [0,1] function we would push the record up significantly and would be able to apply a more strict filtering to the final scoring
The sigmoid function helps app;ly soft thresholding to separate reads below and above the preset threshold value. However does not have a great effect on extremely high values, thus allowing other parameters to make greater effect on the function outcome. Although the effect of the function decreases logarithmically once above the threshold it still has minor effect, which might be useful to make fine sorting.
Improve producing fasta files corresponding to each suggested integration site
Needs to be implemented without relying on the raw fastq files but instead only on the information presented in the sam records
Another idea for the sigmoid function is to allow user to make fine adjustment to the function. The fine adjustments could include the following:
Setting the denominator of the threshold
Setting the minimum penalty as the number of reads (x) approaches 0.

Other observations:
Does a good job at clustering R1 and R2 reads which belong to the same insertion

Main Sources of false positives:
Splicing (currently running spliced aligners)
Sequence similarity between human and hiv (related to splicing)
Defective references

One of the bigger challenges, as can be seen from the comparison of clean1 (multi-fasta reference) and clean2 (lab isolate with known sequence, crispr/cas9 insertion) there are unwanted detections with a multi-fasta reference used. It would help immensely if we could narrow down the reference. Need to speak with Ya-Chi - the question was already brought up a few times, and we will likely narrow down the reference for the final run of the analysis. Yet this still prompts a question of how to remove such false positives.

For a while there was yet another challenge with the analysis pipeline. Namely, the high rate of false positives due to splicing.
However, HISAT2 can not perform local alignments and Bowtie2 can not perform spliced alignments
The solution was to run hisat2 --very-sensitive, end-to-end, given a splicing infile generated from the UCSC gtf annotation.
Set of all hisat2 aligned reads with splicing detected (N present in CIGAR string) were subtracted from the set of aligned reads generated by bowtie2
The results was a great decrease in the rate of false positives

Another possible improvement which has the potential to reduce false positives (but has a potential to introduce false negatives):
Remove alignments where the end of the alignment opposite of the chimeric position has “significant” soft-clipping
This option would be similar to end-to-end
Meaning that the chimera should be aligned end-to-end, which is the case for many
Ex. HUM:9S32M110S HIV:41S110M would be discarded if the proposed parameter is set to less than 9
Via this option the argument would simply filter out anything that has 
Here is an example read:
CGCCAGACTGGAGTGCAATGGTATGATGTCGGCTCACTGCACTCCAGCCTGGGCGACAGTGCCAGACTCCATTTCAAAAAATAAATAAATAAATAAAAGAATAAGTTGCACAGAAAGAGCTTTTTTAAATTAGGAAGAATAGAAACCAAAT
For this read the HIV alignment is 40-50nts short of the end-to-end.
Could be splicing
It is confusing, since the alignment is reported by the UCSC but far behind the spliced alignment
Perhaps that is due to the second pair
Perhaps running single-end alignment would be better in this instance
But interestingly HISAT2 does not align as spliced either, even provided the spliced site

Another option would be to introduce percent identity parameter, which would allow to specify the minimum percent identity required for reporting the chimeric read.
This parameter could also be discussed with Ya-Chi to be set in accordance with her manual verification workflow.
The necessary information could be computed from the optional MD:Z field in the sam output, if present

I have sploken with a student who is participating in the research conducted within the ALIVE study at the JHU.
They too have RNA-seq data, and might be interested in applying this tool.
Once this workflow is finalized, I am very curious to run it on the HIV_HCV coinfection dataset provided by Ashwin, to see if we can detect any HIV chimeras.
However reads are very short (100nts)
Previous analysis did not yield promising results

Bowtie still sometimes produces results which are very different from the UCSC example:

mRNA:FANCI
chr15
R2
164
7.0
22.000000
151.0
NB501749:135:HK2VTBGX3:2:12106:23891:6642
1
0.896241
0.830482
1.0
1769:89314828@gi|3163929|emb|AJ006287.1|:chr15
31.000000
0.000000
NB501749:135:HK2VTBGX3:3:12403:10023:3674;NB50...
51.000000
0.998308
0.750000
0.577938
0.4362

When the human alignment is extracted and aligned with the UCSC BLAT - a single hit to a different chromosome
When full read is aligned - multiple results not a single one matches bowtie2

Yet another idea is to pass inserted genome size as a parameter:
If two unique suggested integration sites are found on the same chromosome with a relative distance in compliance with the size of the inserted genome, we would be more confident in the suggestion and could report such cases separately as very likely
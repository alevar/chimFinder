#!/usr/bin/env python

import sys
import argparse
import Step1
import Step2
import Step3

def hiv(argv):

    parser = argparse.ArgumentParser(description='''Help Page''')
    subparsers = parser.add_subparsers(help='sub-command help')

#==========================================
#==================STEP1===================
#==== Take two alignments and output ======
#=========== suggested chimeras ===========
#==========================================
#./hiv.py step1 -i1 hiv.sam -i2 hum.sam -o outDir-t 12 --minLen 30 -a annotation -w --spliced hum.hisat.sam

    parser_step1 = subparsers.add_parser('step1',
                                help='step1 help')
    parser_step1.add_argument('-i1',
                                '--input1',
                                required=True,
                                type=str,
                                help="first alignment")
    parser_step1.add_argument('-i2',
                                '--input2',
                                required=True,
                                type=str,
                                help="second alignment")
    parser_step1.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./out",
                                help="output file")
    parser_step1.add_argument('-t',
                              '--threads',
                              required=False,
                              default=1,
                              type=int,
                              help="the number of threads to use in the computation")
    parser_step1.add_argument('--minCount',
                              required=False,
                              default=5,
                              type=float,
                              help="the minimum number of reads supporting an integration site.")
    parser_step1.add_argument('--maxCountPenalty',
                              required=False,
                              default=0.85,
                              type=float,
                              help="the maximum penalty to give for minCount when hit.")
    parser_step1.add_argument('--weightCount',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the number of reads score in the cumulative score equation")
    parser_step1.add_argument('--minEntropy',
                              required=False,
                              default=0.75,
                              type=float,
                              help="minimum alignment entropy score.")
    parser_step1.add_argument('--weightEntropy',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the entropy score in the cumulative score equation")
    parser_step1.add_argument('--minQual',
                              required=False,
                              default=10,
                              type=float,
                              help="minimum mean mapping quality of sequenced read. Everything below this threshold will be reported as multialignment for which no assumption can be made from the annotation")
    parser_step1.add_argument('--weightQual',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the mapping score in the cumulative score equation")
    parser_step1.add_argument('--minLen',
                              required=False,
                              default=30,
                              type=float,
                              help="the minimum number of nucleotides in alignment to keep a read.")
    parser_step1.add_argument('--weightLen',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the alignment length score in the cumulative score equation")
    parser_step1.add_argument('--steepSlopeAL',
                              required=False,
                              default=0.1,
                              type=float,
                              help="Weight of the alignment length score in the cumulative score equation")
    parser_step1.add_argument('--maxAlLenPenalty',
                              required=False,
                              default=0.0,
                              type=float,
                              help="the maximum penalty to give for minLen when hit.")
    parser_step1.add_argument('-s',
                              '--score',
                              required=False,
                              default=0.6,
                              type=float,
                              help="the minimum overall score to keep.")
    parser_step1.add_argument('--close',
                              required=False,
                              default=30,
                              type=int,
                              help="distance between two integration sites to group together")
    parser_step1.add_argument('-a',
                              '--annotation',
                              required=True,
                              type=str,
                              help="annotation for the human genome")
    parser_step1.add_argument('-e',
                              '--end',
                              default="",
                              type=str,
                              help="suffix to append to the end of the ouput name")
    parser_step1.add_argument('-w',
                              '--writeReads',
                              action="store_true",
                              help="write reads to fasta files")
    parser_step1.add_argument('-p',
                              '--plot',
                              action="store_true",
                              help="plot snapshots")
    parser_step1.add_argument('--spliced',
                              required=False,
                              type=str,
                              help="spliced end-to-end alignment. Reads will be subtracted from the main alignments and will not be reported in the final report of integrations sites")
    parser_step1.set_defaults(func=Step1.main)

#========================================
#==================STEP2=================
#========================================

    parser_step2 = subparsers.add_parser('step2',
                                help='step2 help')
    parser_step2.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                help="path to the fastq.gz files")
    parser_step2.add_argument('-o',
                                '--out',
                                required=False,
                                default="./out",
                                type=str,
                                help="path to the output directory")
    parser_step2.add_argument('-m',
                              '--minLen',
                              type=int,
                              required=False,
                              default=50,
                              help="minimum length of the n segment")
    parser_step2.add_argument('-l',
                              '--len',
                              type=int,
                              required=False,
                              default=140,
                              help="maximum of contigs to consider")
    parser_step2.add_argument('-c',
                              '--cLen',
                              type=int,
                              required=False,
                              default=50,
                              help="minimum length of contig to consider")
    parser_step2.set_defaults(func=Step2.main)

#========================================
#==================STEP3=================
#========================================

    parser_step3 = subparsers.add_parser('step3',
                                help='step3 help')
    parser_step3.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                help="path to the fastq.gz files")
    parser_step3.add_argument('-o',
                                '--out',
                                required=False,
                                default="./out",
                                type=str,
                                help="path to the output directory")
    parser_step3.set_defaults(func=Step3.main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    hiv(sys.argv[1:])

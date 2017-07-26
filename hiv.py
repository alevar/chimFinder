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
#==Run Krakem, Hisat2, Bowtie2 and other===
#==tools and output unfiltered data =======
#==========================================
#./hiv.py step1 -i ./data/YHo060517 -k ./customDB -v ./refs/hiv89.6/hiv89.6 -g ./refs/bowtie2/hg38 -a /ccb/salz7-data/genomes/hg38/annotation/hg38_p8.refseq.gff3 -t 12 -o ./out 

    parser_step1 = subparsers.add_parser('step1',
                                help='step1 help')
    parser_step1.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                help="directory which contains fastq.gz files")
    parser_step1.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./out",
                                help="output directory")
    parser_step1.add_argument('-k',
                                '--krakenDB',
                                required=True,
                                type=str,
                                help="path to kraken database")
    parser_step1.add_argument('-v',
                                '--hivDB',
                                required=True,
                                type=str,
                                help="path to the hiv reference")
    parser_step1.add_argument('-g',
                                '--humDB',
                                required=True,
                                type=str,
                                help="path to the hg38 reference")
    parser_step1.add_argument('-t',
                              '--threads',
                              required=False,
                              default=1,
                              type=int,
                              help="the number of threads to use in the computation")
    parser_step1.add_argument('-x',
                              '--minLen1',
                              required=False,
                              default="30:30",
                              type=str,
                              help="the minimum number of nucleotides in alignment to keep a read. REF1:REF2")
    parser_step1.add_argument('-y',
                              '--minLen2',
                              required=False,
                              default="0:0",
                              type=str,
                              help="the minimum number of nucleotides in alignment to keep a read. REF1:REF2")
    parser_step1.add_argument('-a',
                              '--annotation',
                              required=True,
                              type=str,
                              help="annotation for the human genome")
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

#=========================================
#==================STEP3==================
#=========================================

    parser_step3 = subparsers.add_parser('step3',
                                help='step3 help')
    parser_step3.add_argument('-i',
                                '--input',
                                required=False,
                                type=str,
                                help="path to the stats.log generated by or in the format of eva assemble")
    parser_step3.set_defaults(func=Step3.main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    hiv(sys.argv[1:])

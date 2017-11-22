#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import shutil
import sys
import glob
import time
import itertools
import argparse
import scipy
import math

def main():
	# UCSC Blat approach
	columns=["matches",
	"misMatches",
	"repMatches",
	"nCount",
	"qNumInsert",
	"qBaseInsert",
	"tNumInsert",
	"tBaseInsert",
	"strand",
	"qName",
	"qSize",
	"qStart",
	"qEnd",
	"tName",
	"tSize",
	"tStart",
	"tEnd",
	"blockCount",
	"blockSizes",
	"qStarts",
	"tStarts"]
	data=pd.read_csv("./out.psl",sep="\t",comment='@',names=columns)
	toRemove=list(set(data[data["matches"]>100]["qName"]))

	fout=open("./seqFiltered.fasta","w+")
	with open("./sequence.fasta") as fin:
	    curSeq=""
	    seqLines=fin.readlines()
	    for ln in seqLines:
	        if ln[0] == '>':  # new sequence
	            name = ln[1:].split()[0]
	            curSeq=name
	            if curSeq not in toRemove:
	                fout.write(ln)
	        else:
	            if curSeq not in toRemove:
	                fout.write(ln)
	                
	fout.close()

if __name__=="__main__":
	main()
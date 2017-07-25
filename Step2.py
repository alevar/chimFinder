import pandas as pd
import numpy as np
import numba
import os
import signal
import multiprocessing
import shutil
import sys
import glob
import time
import itertools
import argparse
import scipy
import warnings
warnings.filterwarnings('ignore')

# get positions where n's occur
def getContigPos(consensus):
    consensus=consensus.translate({ord(c): None for c in '\n'})
    if "n" in consensus:
        enumeration=[i for i,x in enumerate(consensus) if ((x=="n") and not i==len(consensus)-1 and not consensus[i+1]=="n") or ((x=="n" and not consensus[i-1]=="n"))]
        if not consensus[len(consensus)-1]=="n":
            enumeration.append(len(consensus)-1)
        pairs=list(zip(enumeration[0::2],enumeration[1::2]))
        return pairs
    else:
        return []

# filter out short distances.
# for instance if the threshold is set to 70 - do not count any contigPairs which such limited length
def filterDistances(contigPairs,threshold):
    for contigPair in contigPairs:
        if contigPair[1][0]-contigPair[0][1]<threshold:
            contigPairs.remove(contigPair)
    return contigPairs

# filter contigs smaller in length than set by a parameter
def filterLength(x,y):
	return x

# build all possible forward pairs of positions which could be potential links between splices
def buildPairs(contigs):
    return list(itertools.combinations(contigs, 2))

# trim positions to specified length on each side
def trimContigPos(contigsPaired,length,totalLen):
    newContigList=[]
    for contigPair in contigsPaired:
        if contigPair[0][1]>length and totalLen-contigPair[1][0]>length:
            newcontigPair=((contigPair[0][1]-length,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]<=length and totalLen-contigPair[1][0]>length:
            newcontigPair=((contigPair[0][0]+1,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]>length and totalLen-contigPair[1][0]<=length:
            newcontigPair=((contigPair[0][1]-length,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][1]))
            newContigList.append(newcontigPair)

        if contigPair[0][1]<=length and totalLen-contigPair[1][0]<=length:
            newcontigPair=((contigPair[0][0]+1,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][1]))
            newContigList.append(newcontigPair)

    return newContigList

# extract string for each contig
def getStrings(coordinates,sequence):
    return [(sequence[cord[0][0]:cord[0][1]],sequence[cord[1][0]:cord[1][1]],cord[0][1],cord[1][0]) for cord in coordinates]

# build a possible spliced read
def getUnifiedRef(contigs):
    return [("".join([contigPair[0],contigPair[1]]),contigPair[2],contigPair[3]) for contigPair in contigs]

# write splicesitedata to fasta files to be used as reference in alignment
def writeRefs(finalOut,finalOutUni,header,outDir):
	# create a directory for the sample
    # in the directory save each possible splice fasta as endPosPre:startPosPost
    if not os.path.exists(os.path.abspath(outDir)):
        os.mkdir(os.path.abspath(outDir))

    # write unified splice sites first
    print(finalOutUni)
    for spliceSite in finalOutUni:
        spliceSiteUniPath=os.path.abspath(outDir)+"/u-"+str(spliceSite[1])+":"+str(spliceSite[2])
        readsFile=open(spliceSiteUniPath,'w+')
        readsFile.write(header)
        readsFile.write(spliceSite[0])
        readsFile.close()

    print(finalOut)
    for spliceSite in finalOut:
        spliceSitePath1=os.path.abspath(outDir)+"/i-1_"+str(spliceSite[2])+":"+str(spliceSite[3])
        spliceSitePath2=os.path.abspath(outDir)+"/i-2_"+str(spliceSite[2])+":"+str(spliceSite[3])
        readsFile1=open(spliceSitePath1,'w+')
        readsFile1.write(header)
        readsFile1.write(spliceSite[0])
        readsFile1.close()
        readsFile2=open(spliceSitePath2,'w+')
        readsFile2.write(header)
        readsFile2.write(spliceSite[1])
        readsFile2.close()

# run the shell script to build indeces and perform alignments
def indexAlign(outDir,finalOut):
    for pos in finalOut:
        os.system("./spliceDetect.sh"+" "++""++" "++" "++" ")

# read .sam files and compare results
def confirm(finalOut,finalOutUni):
    pass

def confirmSplice(filePath,args):
    header=""
    finalOut=[]
    finalOutUni=[]
    count=0
    for line in open(filePath,"r"):
        if count==0:
            header=line
            # write first line to the new file
            count=count+1
        else:
            contigs=getContigPos(line)
            if len(contigs)>0:
                contigsPaired=buildPairs(contigs)
                print(contigsPaired)
                contigsPaired=filterDistances(contigsPaired,args.minLen)
                contigsPaired=filterLength(contigsPaired,args.cLen)
                print(contigsPaired)
                contigsPairedTrimmed=trimContigPos(contigsPaired,args.len,len(line))
                contigsStrs=getStrings(contigsPairedTrimmed,line)
                contigStrsUni=getUnifiedRef(contigsStrs)
                finalOut=contigsStrs
                finalOutUni=contigStrsUni

    writeRefs(finalOut,finalOutUni,header,args.out)
    indexAlign(args.out,finalOut)
    return [(x[1],x[2]) for x in finalOutUni]

def main(args):
    poss=confirmSplice("./t.fa",args)
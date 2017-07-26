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
import re
import warnings
warnings.filterwarnings('ignore')

# get positions where n's occur
def getContigPos(consensus):
    consensus=consensus.replace("\n","")
    if "n" in consensus:
        #enumeration=[i for i,x in enumerate(consensus) if ((x=="n") and not i==len(consensus)-1 and not consensus[i+1]=="n") or (x=="n" and not consensus[i-1]=="n") and not i==0]
        #if not consensus[len(consensus)-1]=="n":
        #    enumeration.append(len(consensus)-1)
        #pairs=list(zip(enumeration[0::2],enumeration[1::2]))
        #return pairs
        return [(m.start(0)-1, m.end(0)) for m in re.finditer(re.compile(r'([^n]+)'), consensus)]
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
def filterLength(contigs,length):
	return [contig for contig in contigs if contig[1]-contig[0]>length]

# build all possible forward pairs of positions which could be potential links between splices
def buildPairs(contigs):
    return list(itertools.combinations(contigs, 2))

# trim positions to specified length on each side
def trimContigPos(contigsPaired,length,totalLen):
    newContigList=[]
    for contigPair in contigsPaired:
        ex1=contigPair[0][1]-contigPair[0][0]
        ex2=contigPair[1][1]-contigPair[1][0]
        if contigPair[0][1]>length and ex1>length and totalLen-contigPair[1][0]>length and ex2>length:
            newcontigPair=((contigPair[0][1]-length,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]<=length and totalLen-contigPair[1][0]>length and ex2>length:
            newcontigPair=((contigPair[0][0]+1,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]>length and ex1>length and totalLen-contigPair[1][0]<=length:
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
def writeRefs(finalOut,finalOutUni,header,outDir,baseName):
    outDir=os.path.abspath(outDir)+"/spliceRefs/"
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    outDir=outDir+baseName
	# create a directory for the sample
    # in the directory save each possible splice fasta as endPosPre:startPosPost
    if not os.path.exists(os.path.abspath(outDir)):
        os.mkdir(os.path.abspath(outDir))

    # write unified splice sites first
    for spliceSite in finalOutUni:
        spliceSiteUniPath=os.path.abspath(outDir)+"/u-"+str(spliceSite[1])+":"+str(spliceSite[2])+".fa"
        readsFile=open(spliceSiteUniPath,'w+')
        readsFile.write(header)
        readsFile.write(spliceSite[0])
        readsFile.close()

    print(finalOut)
    for spliceSite in finalOut:
        spliceSitePath1=os.path.abspath(outDir)+"/i-1_"+str(spliceSite[2])+":"+str(spliceSite[3])+".fa"
        spliceSitePath2=os.path.abspath(outDir)+"/i-2_"+str(spliceSite[2])+":"+str(spliceSite[3])+".fa"
        readsFile1=open(spliceSitePath1,'w+')
        readsFile1.write(header)
        readsFile1.write(spliceSite[0])
        readsFile1.close()
        readsFile2=open(spliceSitePath2,'w+')
        readsFile2.write(header)
        readsFile2.write(spliceSite[1])
        readsFile2.close()

# run the shell script to build indeces and perform alignments
def indexAlign(args,finalOut,baseName,posName):
    os.system("./spliceDetect.sh "+os.path.abspath(args.out)+" "+baseName+" "+args.input+" "+posName)

# read .sam files and compare results
def confirm(finalOut,finalOutUni):
    pass

def confirmSplice(baseName,filePath,args):
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
            contigs=filterLength(contigs,args.cLen)
            if len(contigs)>0:
                contigsPaired=buildPairs(contigs)
                contigsPaired=filterDistances(contigsPaired,args.minLen)
                contigsPairedTrimmed=trimContigPos(contigsPaired,args.len,len(line))
                contigsStrs=getStrings(contigsPairedTrimmed,line)
                contigStrsUni=getUnifiedRef(contigsStrs)
                finalOut=contigsStrs
                finalOutUni=contigStrsUni

    writeRefs(finalOut,finalOutUni,header,args.out,baseName)

    os.system("./prepIndexAlign.sh "+os.path.abspath(args.out)+" "+baseName+" "+args.input)
    for pos in finalOutUni:
        indexAlign(args,finalOut,baseName,str(pos[1])+":"+str(pos[2]))
    return [(x[1],x[2]) for x in finalOutUni]

def main(args):
    for fileN in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
        fullPath=os.path.abspath(fileN)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])
        baseName="_R1".join(fileName.split("_R1")[:-1])
        print("Working on: "+baseName)
        consensusPath=os.path.abspath(args.out)+"/consensusHIV/"+baseName+".fq.fa"
        if os.path.exists(consensusPath):
            poss=confirmSplice(baseName,consensusPath,args)

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

#This step is dedicated to producing final stats

#the function below produces a summarry of the samples
# The following format shall be used for the results page:
# 1.  - sampleName              - name of the sample
# 2.  - krakenHIV%              - % reads classified as Human immunodeficiency virus type 1 by Kraken
# 3.  - krakenHUM%              - % reads classified as Homo sapience by Kraken
# 4.  - chimericKrakenExtracted - number of unique chimeric reads extracted based on the kraken output
# 5.  - chimericKrakenFiltered  - number of unique chimeric reads after aligning with bowtie2 and filtering sam output with a python script
# 6.  - bowtie2HIV              - number of unique reads aligned to the HIV89.6 reference using bowtie2
# 7.  - bowtie2HUM              - number of unique reads aligned to the hg38 using bowtie2
# 8.  - chimericBowtie          - number of unique chimeric reads detected using bowtie2 alignments of all reads against both hg38 and hiv89.6 references
# 9.  - numSplits               - number of unique split points calculated from chimericKraken and chimericBowtie2 reads
# 10. - numReads                - number of unique chimericKraken and chimericBowtie2 reads which support potential split locations
# 11. - numSpliceJunctionsHIV   - number of unique splice junctions from Hisat2 output
# 12. - totalNumberReads        - total number of raw reads in the sample

# what needs to be done
# if a file such as *Pos.csv is absent - record

def getKrakenP(baseName,outDir):
    krakenFilePath=os.path.abspath(outDir)+"/krakenOut/"+baseName+".report"
    hiv=0
    hum=0
    if os.path.exists(krakenFilePath):
        for line in open(krakenFilePath,'r'):
            if 'Human immunodeficiency' in line:
                hiv=line.split('\t')[0]
            if 'Homo sapiens' in line:
                hum=line.split('\t')[0]
    return [hiv,hum]

def getNumChimericExtracted(baseName,outDir):
    selectedFilePath=os.path.abspath(outDir)+'/krakenOut/selected/'+baseName+'.chim'
    nReadsKraken=0
    if os.path.exists(selectedFilePath):
        i=0
        with open(selectedFilePath) as f:
            for i, l in enumerate(f):
                pass
        nReadsKraken=i+1
    return nReadsKraken

def getPostFiltKraken(baseName,outDir):
    postFiltKrakenPath=os.path.abspath(outDir)+'/'+baseName+'.chim.csv'
    postFiltKraken=0
    if os.path.exists(postFiltKrakenPath):
        i=0
        with open(postFiltKrakenPath) as f:
            for i, l in enumerate(f):
                pass
        postFiltKraken=i
    return postFiltKraken

def getFromPos(baseName,outDir):
    splitsPosPath=os.path.abspath(outDir)+'/'+baseName+'_Pos.csv'
    numSplits=0
    numReads=0
    if os.path.exists(splitsPosPath):
        dataPos=pd.read_csv(splitsPosPath)
        numSplits=len(dataPos)
        numReads=dataPos['totalCount'].sum()
    return [numSplits,numReads]

def getBowtieReadsHivFull(baseName,outDir):
    bowtieHivFullPath=os.path.abspath(outDir)+'/fullAlignments/'+baseName+'.full,hiv.sam'
    numBowtieHivFull=0
    if os.path.exists(bowtieHivFullPath):
        dataPos=pd.read_csv(bowtieHivFullPath,sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    return numBowtieHivFull

def getBowtieReadsHivLocal(baseName,outDir):
    bowtieHivLocalPath=os.path.abspath(outDir)+'/localAlignments/'+baseName+'.chim,hiv.sam'
    numBowtieHivLocal=0
    if os.path.exists(bowtieHivLocalPath):
        dataPos=pd.read_csv(bowtieHivLocalPath,sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    return numBowtieHivLocal

def getBowtieReadsHumFull(baseName,outDir):
    bowtieHumFullPath=os.path.abspath(outDir)+'/fullAlignments/'+baseName+'.full.hum.sam'
    numBowtieHumFull=0
    if os.path.exists(bowtieHumFullPath):
        dataPos=pd.read_csv(bowtieHumFullPath,sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    return numBowtieHumFull

def getBowtieReadsHumLocal(baseName,outDir):
    bowtieHumLocalPath=os.path.abspath(outDir)+'/localAlignments/'+baseName+'.chim.hum.sam'
    numBowtieHumLocal=0
    if os.path.exists(bowtieHumLocalPath):
        dataPos=pd.read_csv(bowtieHumLocalPath,sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    return numBowtieHumLocal

def getChimericBowtie(baseName,outDir):
    chimericBowtiePath=os.path.abspath(outDir)+'/'+baseName+'.full.csv'
    chimericBowtie=0
    if os.path.exists(chimericBowtiePath):
        i=0
        with open(chimericBowtiePath) as f:
            for i, l in enumerate(f):
                pass
        chimericBowtie=i
    return chimericBowtie

def main(args):
    data=pd.DataFrame([],columns=["sample",#
                                  "krakenHiv%",#
                                  "krakenHum%",#
                                  "chimericReadsExtracted",#
                                  "chimericReadsFiltered",#
                                  "numBowtieHivFull",#
                                  "numBowtieHivLocal",#
                                  "numBowtieHumFull",#
                                  "numBowtieHumLocal",#
                                  "chimericBowtie",
                                  "numSplits",#
                                  "numReads",#
                                  "numSpliceJunctionsHIV"])
    for fileN in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
        fullPath=os.path.abspath(fileN)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])
        baseName="_R1".join(fileName.split("_R1")[:-1])
        print("Working on: "+baseName)
        
        row=pd.DataFrame([],columns=["sample",#
                                  "krakenHiv%",#
                                  "krakenHum%",#
                                  "chimericReadsExtracted",#
                                  "chimericReadsFiltered",#
                                  "numBowtieHivFull",#
                                  "numBowtieHivLocal",#
                                  "numBowtieHumFull",#
                                  "numBowtieHumLocal",#
                                  "chimericBowtie",
                                  "numSplits",#
                                  "numReads",#
                                  "numSpliceJunctionsHIV"])
        row['sample']=[baseName]
        row[['krakenHiv%','krakenHum%']]=[getKrakenP(baseName,args.out)]
        row['chimericReadsExtracted']=[getNumChimericExtracted(baseName,args.out)]
        row['chimericReadsFiltered']=[getPostFiltKraken(baseName,args.out)]
        row['numBowtieHivFull']=[getBowtieReadsHivFull(baseName,args.out)]
        row['numBowtieHivLocal']=[getBowtieReadsHivLocal(baseName,args.out)]
        row['numBowtieHumFull']=[getBowtieReadsHumFull(baseName,args.out)]
        row['numBowtieHumLocal']=[getBowtieReadsHumLocal(baseName,args.out)]
        row['chimericBowtie']=[getChimericBowtie(baseName,args.out)]
        row[["numSplits","numReads"]]=[getFromPos(baseName,args.out)]
        row['numSpliceJunctionsHIV']=0 #none detected since no consensus sequence is built properly from the multiple reference files
        data=pd.concat([data,row])

    data.to_csv(os.path.abspath(args.out)+'/sampleSummary.csv')
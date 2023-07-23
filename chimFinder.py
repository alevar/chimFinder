#!/usr/bin/env python

#./chimFinder.py --pathogenR1 ./r1/154_pos_21_S7_001.i1.bowtie.sam --hostR1 ./r1/154_pos_21_S7_001.host.bowtie.sam --pathogenR2 ./r2/154_pos_21_S7_001.i1.bowtie.sam --hostR2 ./r2/154_pos_21_S7_001.host.bowtie.sam --splicedR1 ./r1/154_pos_21_S7_001.host.hisat.sam --splicedR2 ./r2/154_pos_21_S7_001.host.hisat.sam -o ./res -t 12 --minLen 20 --maxLenUnmapped 30 -a ./data/hg38_p8.refseq.gff3 --overlap 5 --gap 5 --minEntropy 0.84 --close 5 --score 0.75 -q

import pandas as pd
import numpy as np
import math
import os
import subprocess
import multiprocessing
import signal
import sys
import glob
import itertools
import argparse
from pybedtools import BedTool


gff3cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
sam_colnames = ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL']

reportDF=pd.DataFrame([[]])

def printReport():
    global reportDF
    if reportDF.empty:
        return
    for i in list(reportDF):
        print(i+str(": ")+str(reportDF[i].iloc[0]))

# Instead just write to a report file on separate lines
# those files can then be used to parse by the results.sh script

# Find chimeric reads (even if < 31nt alignment length to pathogen or host) from full alignments
def getFromSplits(datapathogen,dataHost):
    cleanHost=dataHost[(dataHost["paired"]==1)&(dataHost["secondaryAlignment"]==0)&(dataHost["aligned2Mates"]==0)&(dataHost["unmappedCurr"]==0)&(dataHost["unmappedMate"]==8)]
    cleanHost_First=cleanHost[(cleanHost["firstRead"]==64)&(cleanHost["lastRead"]==0)]
    cleanHost_Last=cleanHost[(cleanHost["firstRead"]==0)&(cleanHost["lastRead"]==128)]

    cleanPathogen=datapathogen[(datapathogen["paired"]==1)&(datapathogen["secondaryAlignment"]==0)&(datapathogen["aligned2Mates"]==0)&(datapathogen["unmappedCurr"]==0)&(datapathogen["unmappedMate"]==8)]
    cleanPathogen_First=cleanPathogen[(cleanPathogen["firstRead"]==64)&(cleanPathogen["lastRead"]==0)]
    cleanPathogen_Last=cleanPathogen[(cleanPathogen["firstRead"]==0)&(cleanPathogen["lastRead"]==128)]
    #find mates within paired-end reads which aligned to both pathogen and host
    setPathogen_First=set(cleanPathogen_First["QNAME"])
    setHost_First=set(cleanHost_First["QNAME"])
    splitMates=setPathogen_First.intersection(setHost_First)

# the function below will calculate and return start:end of the alignment within template read
def calcAlignmentStartEnd(row,reference,start):
    # first, break cigar into parseable entities
    cigarList=["".join(x) for _, x in itertools.groupby(row["CIGAR"], key=str.isdigit)]
    if 'M' in cigarList:
        idxfM=cigarList.index("M") #index of the first M
        idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
        if reference==False:
            start=start+sum(list(map(int,cigarList[:idxfM-1][::2])))
            end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
        else:
            end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
        if (row['reversedCurr']==16) and (reference==False):
            return str(len(row["SEQ"])-1-end)+":"+str(len(row["SEQ"])-1-start)
        else:
            return str(start)+":"+str(end)
    else:
        pass

# mark reads that have PATHOGEN on the right side
def leftRight(data,minLen):
    data["PATHOGEN"]=""
    data["PATHOGEN_R1"]=""
    data["PATHOGEN_R2"]=""
    data["sepR1"]=""
    data["sepR2"]=""

    data.loc[~(data["R1_HOST_ID"]==0)&~(data["R1_PATHOGEN_ID"]==0)&(data["R1_PATHOGEN_TS"]-data["R1_HOST_TS"]>minLen)&(data["R1_PATHOGEN_TE"]-data["R1_HOST_TE"]>minLen),"PATHOGEN_R1"]="R1:r"
    data.loc[~(data["R1_HOST_ID"]==0)&~(data["R1_PATHOGEN_ID"]==0)&(data["R1_HOST_TS"]-data["R1_PATHOGEN_TS"]>minLen)&(data["R1_HOST_TE"]-data["R1_PATHOGEN_TE"]>minLen),"PATHOGEN_R1"]="R1:l"

    data.loc[~(data["R2_HOST_ID"]==0)&~(data["R2_PATHOGEN_ID"]==0)&(data["R2_PATHOGEN_TS"]-data["R2_HOST_TS"]>minLen)&(data["R2_PATHOGEN_TE"]-data["R2_HOST_TE"]>minLen),"PATHOGEN_R2"]="R2:r"
    data.loc[~(data["R2_HOST_ID"]==0)&~(data["R2_PATHOGEN_ID"]==0)&(data["R2_HOST_TS"]-data["R2_PATHOGEN_TS"]>minLen)&(data["R2_HOST_TE"]-data["R2_PATHOGEN_TE"]>minLen),"PATHOGEN_R2"]="R2:l"

    data.loc[(~(data["R1_HOST_ID"].astype(str)=="0")&~(data["R2_PATHOGEN_ID"].astype(str)=="0")&(data["R2_HOST_ID"].astype(str)=="0")), "sepR1"]="sepR1:2"
    data.loc[(~(data["R2_HOST_ID"].astype(str)=="0")&~(data["R1_PATHOGEN_ID"].astype(str)=="0")&(data["R1_HOST_ID"].astype(str)=="0")), "sepR2"]="sepR2:1"

    data["PATHOGEN"]=data["PATHOGEN_R1"]+data["PATHOGEN_R2"]+data["sepR1"]+data["sepR2"]
    data=data[data['PATHOGEN'].str.len()>0]
    return data.drop(["PATHOGEN_R1","PATHOGEN_R2","sepR1","sepR2"],axis=1)

def leftRightUnpaired(data,minLen):
    data["PATHOGEN"]=""

    data.loc[~(data["HOST_ID"]==0)&~(data["PATHOGEN_ID"]==0)&(data["PATHOGEN_TS"]-data["HOST_TS"]>minLen)&(data["PATHOGEN_TE"]-data["HOST_TE"]>minLen),"PATHOGEN"]="r"
    data.loc[~(data["HOST_ID"]==0)&~(data["PATHOGEN_ID"]==0)&(data["HOST_TS"]-data["PATHOGEN_TS"]>minLen)&(data["HOST_TE"]-data["PATHOGEN_TE"]>minLen),"PATHOGEN"]="l"

    data=data[data['PATHOGEN'].str.len()>0]
    return data

def overlap(data):
    data.reset_index(inplace=True,drop=True)
    data["overlapR1"]=data[['R1_HOST_TE','R1_PATHOGEN_TE']].min(axis=1)-data[['R1_HOST_TS','R1_PATHOGEN_TS']].max(axis=1) #min end - max start
    data["gapR1"]=0
    data.loc[data['overlapR1']<0,['gapR1']]=data.loc[data['overlapR1']<0,['overlapR1']]["overlapR1"].abs()
    data.loc[data['overlapR1']<0,['overlapR1']]=0
    data["overlapR2"]=data[['R2_HOST_TE','R2_PATHOGEN_TE']].min(axis=1)-data[['R2_HOST_TS','R2_PATHOGEN_TS']].max(axis=1)
    data["gapR2"]=0
    data.loc[data['overlapR2']<0,['gapR2']]=data.loc[data['overlapR2']<0,['overlapR2']]["overlapR2"].abs()
    data.loc[data['overlapR2']<0,['overlapR2']]=0
    return data

def overlapUnpaired(data):
    data.reset_index(inplace=True,drop=True)
    data["overlap"]=data[['HOST_TE','PATHOGEN_TE']].min(axis=1)-data[['HOST_TS','PATHOGEN_TS']].max(axis=1) #min end - max start
    data["gap"]=0
    data.loc[data['overlap']<0,['gap']]=data.loc[data['overlap']<0,['overlap']]["overlap"].abs()
    data.loc[data['overlap']<0,['overlap']]=0
    return data

# extract flag information
def extractFlagBits(data):
    data["paired"]=data["FLAG"]               &1 #template having multiple segments in sequencing
    data["aligned2Mates"]=data["FLAG"]        &2 #each segment properly aligned according to the aligner
    data["unmappedCurr"]=data["FLAG"]         &4 #segment unmapped
    data["unmappedMate"]=data["FLAG"]         &8 #next segment in the template unmapped
    data["reversedCurr"]=data["FLAG"]         &16 #SEQ being reverse complemented
    data["reversedMate"]=data["FLAG"]         &32 #SEQ of the next segment in the template being reverse complemented
    data["firstRead"]=data["FLAG"]            &64 #the first segment in the template
    data["lastRead"]=data["FLAG"]             &128 #the last segment in the template
    data["secondaryAlignment"]=data["FLAG"]   &256 #secondary alignment
    data["noPassFilter"]=data["FLAG"]         &512 #not passing filters, such as platform/vendor quality controls
    data["PCRdup"]=data["FLAG"]               &1024 #PCR or optical duplicate
    data["suppAl"]=data["FLAG"]               &2048 #supplementary alignment

def extractStartEnd(data):
    data["CIGAR"].replace("*",np.nan,inplace=True)
    data.dropna(axis=0,inplace=True)
    data.reset_index(drop=True,inplace=True)

    data["READ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$",expand=False).replace(np.nan,0).astype(int)
    data["END"]=data.READ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[S]",expand=False).replace(np.nan,0).astype(int)

    data16=data[data["reversedCurr"]==16].reset_index(drop=True)
    data0=data[data["reversedCurr"]==0].reset_index(drop=True)
    data16["Template_start"]=data16.READ_LEN-data16.END
    data16["Template_end"]=data16.READ_LEN-data16.CIGAR_PRE-1
    data0["Template_start"]=data0.CIGAR_PRE
    data0["Template_end"]=data0.END

    data16["Reference_start"]=data16.READ_LEN-data16.END+data16.POS-data16.Template_start
    data16["Reference_end"]=data16.READ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE

    data=pd.concat([data16,data0]).reset_index(drop=True)
    data.drop(["CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
    return data

# filtering the reads based on the flags:
def filterReads(dataHost,dataPathogen):
    #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
    dataHost=dataHost[(dataHost["secondaryAlignment"]==0)&(dataHost["PCRdup"]==0)&(dataHost["suppAl"]==0)&(dataHost["noPassFilter"]==0)]
    dataPathogen=dataPathogen[(dataPathogen["secondaryAlignment"]==0)&(dataPathogen["PCRdup"]==0)&(dataPathogen["suppAl"]==0)&(dataPathogen["noPassFilter"]==0)]
    return dataHost, dataPathogen

def createData(data,dataHost,dataPathogen):
    dataHost=dataHost[dataHost['QNAME'].isin(set(data['QNAME']))]
    dataHostR1=pd.DataFrame([])
    dataHostR1[["QNAME",
                "R1_HOST_TS",
                "R1_HOST_TE",
                "R1_HOST_ID",
                "R1_HOST_RS",
                "R1_HOST_RE",
                "READ_LEN_R1",
                "R1_HOST_SEQ",
                "QUAL_R1",
                "R1_HOST_MAPQ",
                "R1_HOST_reversedCurr"]]=dataHost[dataHost['firstRead']==64][["QNAME",
                                                                            "Template_start",
                                                                            "Template_end",
                                                                            "RNAME",
                                                                            "Reference_start",
                                                                            "Reference_end",
                                                                            "READ_LEN",
                                                                            "SEQ",
                                                                            "QUAL",
                                                                            "MAPQ",
                                                                            "reversedCurr"]]
    dataHostR2=pd.DataFrame([])
    dataHostR2[["QNAME",
                "R2_HOST_TS",
                "R2_HOST_TE",
                "R2_HOST_ID",
                "R2_HOST_RS",
                "R2_HOST_RE",
                "READ_LEN_R2",
                "R2_HOST_SEQ",
                "QUAL_R2",
                "R2_HOST_MAPQ",
                "R2_HOST_reversedCurr"]]=dataHost[dataHost['lastRead']==128][["QNAME",
                                                                            "Template_start",
                                                                            "Template_end",
                                                                            "RNAME",
                                                                            "Reference_start",
                                                                            "Reference_end",
                                                                            "READ_LEN",
                                                                            "SEQ",
                                                                            "QUAL",
                                                                            "MAPQ",
                                                                            "reversedCurr"]]
    dataPathogenR1=pd.DataFrame([])
    dataPathogenR1[["QNAME",
                "R1_PATHOGEN_TS",
                "R1_PATHOGEN_TE",
                "R1_PATHOGEN_ID",
                "R1_PATHOGEN_RS",
                "R1_PATHOGEN_RE",
                "R1_PATHOGEN_SEQ",
                "R1_PATHOGEN_MAPQ",
                "R1_PATHOGEN_reversedCurr"]]=dataPathogen[dataPathogen['firstRead']==64][["QNAME",
                                                                            "Template_start",
                                                                            "Template_end",
                                                                            "RNAME",
                                                                            "Reference_start",
                                                                            "Reference_end",
                                                                            "SEQ",
                                                                            "MAPQ",
                                                                            "reversedCurr"]]
    dataPathogenR2=pd.DataFrame([])
    dataPathogenR2[["QNAME",
                "R2_PATHOGEN_TS",
                "R2_PATHOGEN_TE",
                "R2_PATHOGEN_ID",
                "R2_PATHOGEN_RS",
                "R2_PATHOGEN_RE",
                "R2_PATHOGEN_SEQ",
                "R2_PATHOGEN_MAPQ",
                "R2_PATHOGEN_reversedCurr"]]=dataPathogen[dataPathogen['lastRead']==128][["QNAME",
                                                                        "Template_start",
                                                                        "Template_end",
                                                                        "RNAME",
                                                                        "Reference_start",
                                                                        "Reference_end",
                                                                        "SEQ",
                                                                        "MAPQ",
                                                                        "reversedCurr"]]
    data=data.merge(dataHostR1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataHostR2,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataPathogenR1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataPathogenR2,left_on='QNAME', right_on='QNAME', how='left')
    data[["R1_HOST_ID",
        "R2_HOST_ID",
        "R1_PATHOGEN_ID",
        "R2_PATHOGEN_ID",
        "R1_HOST_MAPQ",
        "R2_HOST_MAPQ",
        "R1_PATHOGEN_MAPQ",
        "R2_PATHOGEN_MAPQ",
        "R1_HOST_reversedCurr",
        "R2_HOST_reversedCurr",
        "R1_PATHOGEN_reversedCurr",
        "R2_PATHOGEN_reversedCurr"]]=data[["R1_HOST_ID",
                                    "R2_HOST_ID",
                                    "R1_PATHOGEN_ID",
                                    "R2_PATHOGEN_ID",
                                    "R1_HOST_MAPQ",
                                    "R2_HOST_MAPQ",
                                    "R1_PATHOGEN_MAPQ",
                                    "R2_PATHOGEN_MAPQ",
                                    "R1_HOST_reversedCurr",
                                    "R2_HOST_reversedCurr",
                                    "R1_PATHOGEN_reversedCurr",
                                    "R2_PATHOGEN_reversedCurr"]].fillna('',inplace=True)
    data.fillna(0,inplace=True)
    return data

def createDataUnpaired(data,dataHost,dataPathogen):
    dataHost=dataHost[dataHost['QNAME'].isin(set(data['QNAME']))].reset_index(drop=True)
    dfHost=pd.DataFrame([])
    dfPathogen=pd.DataFrame([])
    dfHost[["QNAME",
            "HOST_TS",
            "HOST_TE",
            "HOST_ID",
            "HOST_RS",
            "HOST_RE",
            "READ_LEN",
            "HOST_SEQ",
            "QUAL",
            "HOST_MAPQ",
            "HOST_reversedCurr"]]=dataHost[["QNAME",
                                            "Template_start",
                                            "Template_end",
                                            "RNAME",
                                            "Reference_start",
                                            "Reference_end",
                                            "READ_LEN",
                                            "SEQ",
                                            "QUAL",
                                            "MAPQ",
                                            "reversedCurr"]]
    dfPathogen[["QNAME",
            "PATHOGEN_TS",
            "PATHOGEN_TE",
            "PATHOGEN_ID",
            "PATHOGEN_RS",
            "PATHOGEN_RE",
            "PATHOGEN_SEQ",
            "PATHOGEN_MAPQ",
            "PATHOGEN_reversedCurr"]]=dataPathogen[["QNAME",
                                            "Template_start",
                                            "Template_end",
                                            "RNAME",
                                            "Reference_start",
                                            "Reference_end",
                                            "SEQ",
                                            "MAPQ",
                                            "reversedCurr"]]
    data=data.merge(dfHost,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dfPathogen,left_on='QNAME', right_on='QNAME', how='left')
    data[["HOST_ID",
        "PATHOGEN_ID",
        "HOST_MAPQ",
        "PATHOGEN_MAPQ",
        "HOST_reversedCurr",
        "PATHOGEN_reversedCurr"]]=data[["HOST_ID",
                                    "PATHOGEN_ID",
                                    "HOST_MAPQ",
                                    "PATHOGEN_MAPQ",
                                    "HOST_reversedCurr",
                                    "PATHOGEN_reversedCurr"]].fillna(value='')
    data.fillna(0,inplace=True)
    return data

memo={}
# the function below should compute normalized entropy as described in https://academic.oup.com/bioinformatics/article/27/8/1061/227307/Topological-entropy-of-DNA-sequences
def topologicalNormalizedEntropy(s):
    def findCF(l):
        for i in range(1,l):
            lt=((4**i)+i)-1 #lesser term of the CF expression
            if lt>l:
                return i-1

    def countSubString(s,l):
        substrings=[]
        for i in range(len(s)-l):
            substrings.append(s[i:i+l])
        return len(set(substrings))

    # calculate the expected entropy
    def expectedEntropy(n):
        fourN=float(4**n)
        logTerm=(fourN-(fourN*((1.0-(1.0/fourN))**fourN)))
        return math.log(logTerm,4)/float(n)

    global memo
    if len(s) in memo:
        maxSubstringLen=memo[len(s)]
    else:
        maxSubstringLen=findCF(len(s))
        memo[len(s)]=maxSubstringLen
    n=countSubString(s,maxSubstringLen)
    res=math.log(n,4)/maxSubstringLen

    return res

# How to compute the mimimum entropy for a string of a given length and characters

#consider using technique described in https://academic.oup.com/bioinformatics/article/27/8/1061/227307/Topological-entropy-of-DNA-sequences
#another paper: Look at figure 2 https://www.nature.com/articles/srep19788#f2
#more https://www.xycoon.com/normalized_entropy.htm

def processAligns(seqHost,seqPathogen,qual,i1_g2,i2_g2,i1_g1,i2_g1,readLen,rHost,rPathogen):

    entropyScore_pathogen=0
    entropyScore_host=0
    meanQual_pathogen=0
    meanQual_host=0

    s_pathogen=""
    if rPathogen==True: # if reverse complemented take into account
        s_pathogen=seqPathogen[readLen-i2_g2:readLen-i1_g2]
    else:
        s_pathogen=seqPathogen[i1_g2:i2_g2]
    if not len(s_pathogen)==0:
        entropyScore_pathogen=topologicalNormalizedEntropy(s_pathogen)
        q_g2=qual[i1_g2:i2_g2]
        if len(q_g2)>0:
            meanQual_pathogen=sum([ord(x)-33 for x in q_g2])/len(q_g2)
        else:
            meanQual_pathogen=0

    s_host=""
    if rHost==True: # if reverse complemented take into account
        s_host=seqHost[readLen-i2_g1:readLen-i1_g1]
    else:
        s_host=seqHost[i1_g1:i2_g1]
    if not len(s_host)==0:
        entropyScore_host=topologicalNormalizedEntropy(s_host)
        q_g1=qual[i1_g1:i2_g1]
        if len(q_g1)>0:
            meanQual_host=sum([ord(x)-33 for x in q_g1])/len(q_g1)
        else:
            meanQual_host=0

    return pd.Series({"entropyScore_pathogen":entropyScore_pathogen,
                        "meanQual_pathogen":meanQual_pathogen,
                        "entropyScore_host":entropyScore_host,
                        "meanQual_host":meanQual_host})


# this function filters by overlap and flanking alignment length
# further it removes unnecessary data and combines into a single full dataframe with unified naming avoiding R1/R2 conventions
# output can be saved as .full.csv and then grouped by split position all at once

# this function allows identifying best reads
# However, since greater minLen values for PATHOGEN and HOST alignments will yield fewer reads but at higher confidence
# support reads could ignore the min len requirement as long as the split position is identical
# or perhaps the minLen requirement for the support reads should be lower
def filterOverlapCombine(data,args):
    dropList=["R1_HOST_TS",
              "R1_HOST_TE",
              "R1_HOST_ID",
              "R1_HOST_RS",
              "R1_HOST_RE",
              "R2_HOST_TS",
              "R2_HOST_TE",
              "R2_HOST_ID",
              "R2_HOST_RS",
              "R2_HOST_RE",
              "R1_PATHOGEN_TS",
              "R1_PATHOGEN_TE",
              "R1_PATHOGEN_ID",
              "R1_PATHOGEN_RS",
              "R1_PATHOGEN_RE",
              "R2_PATHOGEN_TS",
              "R2_PATHOGEN_TE",
              "R2_PATHOGEN_ID",
              "R2_PATHOGEN_RS",
              "R2_PATHOGEN_RE",
              "overlapR1",
              "overlapR2",
              "gapR1",
              "gapR2",
              "R1_HOST_SEQ",
              "R2_HOST_SEQ",
              "R1_PATHOGEN_SEQ",
              "R2_PATHOGEN_SEQ",
              "QUAL_R1",
              "QUAL_R2",
              "R1_HOST_MAPQ",
              "R1_PATHOGEN_MAPQ",
              "R2_HOST_MAPQ",
              "R2_PATHOGEN_MAPQ",
              "READ_LEN_R1",
              "READ_LEN_R2",
              "R1_HOST_reversedCurr",
              "R2_HOST_reversedCurr",
              "R1_PATHOGEN_reversedCurr",
              "R2_PATHOGEN_reversedCurr"]

    # R1 right
    data['entropyScore_pathogen']=0
    data['meanQual_pathogen']=0
    data['mapQual_pathogen']=0
    data['entropyScore_host']=0
    data['meanQual_host']=0
    data['mapQual_host']=0
    dataR1Right=data[data["PATHOGEN"].str.contains("R1:r")]
    dataR1Right=dataR1Right[~((dataR1Right["R1_HOST_TS"]<dataR1Right["R1_PATHOGEN_TS"])&(dataR1Right["R1_HOST_TE"]>dataR1Right["R1_PATHOGEN_TE"]))]
    dataR1Right=dataR1Right[~((dataR1Right["R1_PATHOGEN_TS"]<dataR1Right["R1_HOST_TS"])&(dataR1Right["R1_PATHOGEN_TE"]>dataR1Right["R1_HOST_TE"]))]
    dataR1Right["ins"]=dataR1Right["R1_PATHOGEN_TS"]-dataR1Right["R1_HOST_TE"]
    dataR1Right["split"]=dataR1Right['R1_HOST_RE'].astype(int).astype(str)+":"+dataR1Right['R1_PATHOGEN_RS'].astype(int).astype(str)
    dataR1Right['HOST']=dataR1Right['R1_HOST_RE'].astype(int)
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1_HOST_ID+":"+dataR1Right.R1_PATHOGEN_ID
    dataR1Right["orient"]="R1-host:pathogen"
    dataR1Right["overlap"]=dataR1Right["overlapR1"]
    dataR1Right["gap"]=dataR1Right["gapR1"]
    dataR1Right["HOST_TS"]=dataR1Right["R1_HOST_TS"].astype(int)
    dataR1Right["HOST_TE"]=dataR1Right["R1_HOST_TE"].astype(int)
    dataR1Right["HOST_RS"]=dataR1Right["R1_HOST_RS"].astype(int)
    dataR1Right["HOST_RE"]=dataR1Right["R1_HOST_RE"].astype(int)
    dataR1Right["PATHOGEN_TS"]=dataR1Right["R1_PATHOGEN_TS"].astype(int)
    dataR1Right["PATHOGEN_TE"]=dataR1Right["R1_PATHOGEN_TE"].astype(int)
    dataR1Right["PATHOGEN_RS"]=dataR1Right["R1_PATHOGEN_RS"].astype(int)
    dataR1Right["PATHOGEN_RE"]=dataR1Right["R1_PATHOGEN_RE"].astype(int)
    dataR1Right["HOST_ID"]=dataR1Right["R1_HOST_ID"]
    dataR1Right["PATHOGEN_ID"]=dataR1Right["R1_PATHOGEN_ID"]
    dataR1Right["SEQ"]=dataR1Right["R1_HOST_SEQ"]
#     dataR1Right["SEQ"]=dataR1Right[dataR1Right['R1_HOST_reversedCurr']==0]["R1_HOST_SEQ"]
#     dataR1Right["R_SEQ"]=dataR1Right[dataR1Right['R1_HOST_reversedCurr']==16]["R1_HOST_SEQ"]
    dataR1Right["PATHOGEN_MAPQ"]=dataR1Right["R1_PATHOGEN_MAPQ"].astype(int)
    dataR1Right["HOST_MAPQ"]=dataR1Right["R1_HOST_MAPQ"].astype(int)
    dataR1Right["PATHOGEN_AL"]=dataR1Right["PATHOGEN_TE"]-dataR1Right["PATHOGEN_TS"]-dataR1Right["overlap"]
    dataR1Right["HOST_AL"]=dataR1Right["HOST_TE"]-dataR1Right["HOST_TS"]-dataR1Right["overlap"]
    dataR1Right=dataR1Right[(dataR1Right['PATHOGEN_AL']>args.minLen)&(dataR1Right['HOST_AL']>args.minLen)]
    if len(dataR1Right>0):
        dataR1Right[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataR1Right.merge(dataR1Right.apply(lambda row: processAligns(row['R1_HOST_SEQ'],
                                                                                                    row['R1_PATHOGEN_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['PATHOGEN_TS']+row['overlap']),
                                                                                                    int(row['PATHOGEN_TE']),
                                                                                                    int(row['HOST_TS']),
                                                                                                    int(row['HOST_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1_HOST_reversedCurr"]==16,
                                                                                                    row["R1_PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                                                   'meanQual_pathogen_y',
                                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                                   'meanQual_host_y']]
    dataR1Right["READ_LEN"]=dataR1Right["READ_LEN_R1"]
    dataR1Right["R"]="R1"
    dataR1Right.drop(dropList,axis=1,inplace=True)

    # R1 left
    dataR1Left=data[data["PATHOGEN"].str.contains("R1:l")]
    dataR1Left=dataR1Left[~((dataR1Left["R1_HOST_TS"]<dataR1Left["R1_PATHOGEN_TS"])&(dataR1Left["R1_HOST_TE"]>dataR1Left["R1_PATHOGEN_TE"]))]
    dataR1Left=dataR1Left[~((dataR1Left["R1_PATHOGEN_TS"]<dataR1Left["R1_HOST_TS"])&(dataR1Left["R1_PATHOGEN_TE"]>dataR1Left["R1_HOST_TE"]))]
    dataR1Left["ins"]=dataR1Left["R1_HOST_TS"]-dataR1Left["R1_PATHOGEN_TE"]
    dataR1Left["split"]=dataR1Left['R1_PATHOGEN_RE'].astype(int).astype(str)+":"+dataR1Left['R1_HOST_RS'].astype(int).astype(str)
    dataR1Left['HOST']=dataR1Left['R1_HOST_RS'].astype(int)
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1_PATHOGEN_ID+":"+dataR1Left.R1_HOST_ID
    dataR1Left["orient"]="R1-pathogen:host"
    dataR1Left["overlap"]=dataR1Left["overlapR1"]
    dataR1Left["gap"]=dataR1Left["gapR1"]
    dataR1Left["HOST_TS"]=dataR1Left["R1_HOST_TS"].astype(int)
    dataR1Left["HOST_TE"]=dataR1Left["R1_HOST_TE"].astype(int)
    dataR1Left["HOST_RS"]=dataR1Left["R1_HOST_RS"].astype(int)
    dataR1Left["HOST_RE"]=dataR1Left["R1_HOST_RE"].astype(int)
    dataR1Left["PATHOGEN_TS"]=dataR1Left["R1_PATHOGEN_TS"].astype(int)
    dataR1Left["PATHOGEN_TE"]=dataR1Left["R1_PATHOGEN_TE"].astype(int)
    dataR1Left["PATHOGEN_RS"]=dataR1Left["R1_PATHOGEN_RS"].astype(int)
    dataR1Left["PATHOGEN_RE"]=dataR1Left["R1_PATHOGEN_RE"].astype(int)
    dataR1Left["HOST_ID"]=dataR1Left["R1_HOST_ID"]
    dataR1Left["PATHOGEN_ID"]=dataR1Left["R1_PATHOGEN_ID"]
    dataR1Left["SEQ"]=dataR1Left["R1_HOST_SEQ"]
    dataR1Left["PATHOGEN_MAPQ"]=dataR1Left["R1_PATHOGEN_MAPQ"].astype(int)
    dataR1Left["HOST_MAPQ"]=dataR1Left["R1_HOST_MAPQ"].astype(int)
    dataR1Left["PATHOGEN_AL"]=dataR1Left["PATHOGEN_TE"]-dataR1Left["PATHOGEN_TS"]-dataR1Left["overlap"]
    dataR1Left["HOST_AL"]=dataR1Left["HOST_TE"]-dataR1Left["HOST_TS"]-dataR1Left["overlap"]
    dataR1Left=dataR1Left[(dataR1Left["PATHOGEN_AL"]>args.minLen)&(dataR1Left["HOST_AL"]>args.minLen)]
    if len(dataR1Left)>0:
        dataR1Left[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataR1Left.merge(dataR1Left.apply(lambda row: processAligns(row['R1_HOST_SEQ'],
                                                                                                    row['R1_PATHOGEN_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['PATHOGEN_TS']),
                                                                                                    int(row['PATHOGEN_TE']-row['overlap']),
                                                                                                    int(row['HOST_TS']+row['overlap']),
                                                                                                    int(row['HOST_TE']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1_HOST_reversedCurr"]==16,
                                                                                                    row["R1_PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                                   'meanQual_pathogen_y',
                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                   'meanQual_host_y']]

    dataR1Left["READ_LEN"]=dataR1Left["READ_LEN_R1"]
    dataR1Left["R"]="R1"
    dataR1Left.drop(dropList,axis=1,inplace=True)

    # R2 right
    dataR2Right=data[data["PATHOGEN"].str.contains("R2:r")]
    dataR2Right=dataR2Right[~((dataR2Right["R2_HOST_TS"]<dataR2Right["R2_PATHOGEN_TS"])&(dataR2Right["R2_HOST_TE"]>dataR2Right["R2_PATHOGEN_TE"]))]
    dataR2Right=dataR2Right[~((dataR2Right["R2_PATHOGEN_TS"]<dataR2Right["R2_HOST_TS"])&(dataR2Right["R2_PATHOGEN_TE"]>dataR2Right["R2_HOST_TE"]))]
    dataR2Right["ins"]=dataR2Right["R2_PATHOGEN_TS"]-dataR2Right["R2_HOST_TE"]
    dataR2Right["split"]=dataR2Right['R2_HOST_RE'].astype(int).astype(str)+":"+dataR2Right['R2_PATHOGEN_RS'].astype(int).astype(str)
    dataR2Right['HOST']=dataR2Right['R2_HOST_RE'].astype(int)
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2_HOST_ID+":"+dataR2Right.R2_PATHOGEN_ID
    dataR2Right["orient"]="R2-host:pathogen"
    dataR2Right["overlap"]=dataR2Right["overlapR2"]
    dataR2Right["gap"]=dataR2Right["gapR2"]
    dataR2Right["HOST_TS"]=dataR2Right["R2_HOST_TS"].astype(int)
    dataR2Right["HOST_TE"]=dataR2Right["R2_HOST_TE"].astype(int)
    dataR2Right["HOST_RS"]=dataR2Right["R2_HOST_RS"].astype(int)
    dataR2Right["HOST_RE"]=dataR2Right["R2_HOST_RE"].astype(int)
    dataR2Right["PATHOGEN_TS"]=dataR2Right["R2_PATHOGEN_TS"].astype(int)
    dataR2Right["PATHOGEN_TE"]=dataR2Right["R2_PATHOGEN_TE"].astype(int)
    dataR2Right["PATHOGEN_RS"]=dataR2Right["R2_PATHOGEN_RS"].astype(int)
    dataR2Right["PATHOGEN_RE"]=dataR2Right["R2_PATHOGEN_RE"].astype(int)
    dataR2Right["HOST_ID"]=dataR2Right["R2_HOST_ID"]
    dataR2Right["PATHOGEN_ID"]=dataR2Right["R2_PATHOGEN_ID"]
    dataR2Right["SEQ"]=dataR2Right["R2_HOST_SEQ"]
    dataR2Right["PATHOGEN_MAPQ"]=dataR2Right["R2_PATHOGEN_MAPQ"].astype(int)
    dataR2Right["HOST_MAPQ"]=dataR2Right["R2_HOST_MAPQ"].astype(int)
    dataR2Right["PATHOGEN_AL"]=dataR2Right["PATHOGEN_TE"]-dataR2Right["PATHOGEN_TS"]-dataR2Right["overlap"]
    dataR2Right["HOST_AL"]=dataR2Right["HOST_TE"]-dataR2Right["HOST_TS"]-dataR2Right["overlap"]
    dataR2Right=dataR2Right[(dataR2Right["PATHOGEN_AL"]>args.minLen)&(dataR2Right["HOST_AL"]>args.minLen)]
    if len(dataR2Right)>0:
        dataR2Right[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataR2Right.merge(dataR2Right.apply(lambda row: processAligns(row['R2_HOST_SEQ'],
                                                                                                    row['R2_PATHOGEN_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['PATHOGEN_TS']+row['overlap']),
                                                                                                    int(row['PATHOGEN_TE']),
                                                                                                    int(row['HOST_TS']),
                                                                                                    int(row['HOST_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2_HOST_reversedCurr"]==16,
                                                                                                    row["R2_PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                                                   'meanQual_pathogen_y',
                                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                                   'meanQual_host_y']]
    dataR2Right["READ_LEN"]=dataR2Right["READ_LEN_R2"]
    dataR2Right["R"]="R2"
    dataR2Right.drop(dropList,axis=1,inplace=True)

    # R2 left
    dataR2Left=data[data["PATHOGEN"].str.contains("R2:l")]
    dataR2Left=dataR2Left[~((dataR2Left["R2_HOST_TS"]<dataR2Left["R2_PATHOGEN_TS"])&(dataR2Left["R2_HOST_TE"]>dataR2Left["R2_PATHOGEN_TE"]))]
    dataR2Left=dataR2Left[~((dataR2Left["R2_PATHOGEN_TS"]<dataR2Left["R2_HOST_TS"])&(dataR2Left["R2_PATHOGEN_TE"]>dataR2Left["R2_HOST_TE"]))]
    dataR2Left["ins"]=dataR2Left["R2_HOST_TS"]-dataR2Left["R2_PATHOGEN_TE"]
    dataR2Left["split"]=dataR2Left['R2_PATHOGEN_RE'].astype(int).astype(str)+":"+dataR2Left['R2_HOST_RS'].astype(int).astype(str)
    dataR2Left['HOST']=dataR2Left['R2_HOST_RS'].astype(int)
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2_PATHOGEN_ID+":"+dataR2Left.R2_HOST_ID
    dataR2Left["orient"]="R2-pathogen:host"
    dataR2Left["overlap"]=dataR2Left["overlapR2"]
    dataR2Left["gap"]=dataR2Left["gapR2"]
    dataR2Left["HOST_TS"]=dataR2Left["R2_HOST_TS"].astype(int)
    dataR2Left["HOST_TE"]=dataR2Left["R2_HOST_TE"].astype(int)
    dataR2Left["HOST_RS"]=dataR2Left["R2_HOST_RS"].astype(int)
    dataR2Left["HOST_RE"]=dataR2Left["R2_HOST_RE"].astype(int)
    dataR2Left["PATHOGEN_TS"]=dataR2Left["R2_PATHOGEN_TS"].astype(int)
    dataR2Left["PATHOGEN_TE"]=dataR2Left["R2_PATHOGEN_TE"].astype(int)
    dataR2Left["PATHOGEN_RS"]=dataR2Left["R2_PATHOGEN_RS"].astype(int)
    dataR2Left["PATHOGEN_RE"]=dataR2Left["R2_PATHOGEN_RE"].astype(int)
    dataR2Left["HOST_ID"]=dataR2Left["R2_HOST_ID"]
    dataR2Left["PATHOGEN_ID"]=dataR2Left["R2_PATHOGEN_ID"]
    dataR2Left["SEQ"]=dataR2Left["R2_HOST_SEQ"]
    dataR2Left["PATHOGEN_MAPQ"]=dataR2Left["R2_PATHOGEN_MAPQ"].astype(int)
    dataR2Left["HOST_MAPQ"]=dataR2Left["R2_HOST_MAPQ"].astype(int)
    dataR2Left["PATHOGEN_AL"]=dataR2Left["PATHOGEN_TE"]-dataR2Left["PATHOGEN_TS"]-dataR2Left["overlap"]
    dataR2Left["HOST_AL"]=dataR2Left["HOST_TE"]-dataR2Left["HOST_TS"]-dataR2Left["overlap"]
    dataR2Left=dataR2Left[(dataR2Left["PATHOGEN_AL"]>args.minLen)&(dataR2Left["HOST_AL"]>args.minLen)]
    if len(dataR2Left)>0:
        dataR2Left[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataR2Left.merge(dataR2Left.apply(lambda row: processAligns(row['R2_HOST_SEQ'],
                                                                                                    row['R2_PATHOGEN_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['PATHOGEN_TS']),
                                                                                                    int(row['PATHOGEN_TE']-row['overlap']),
                                                                                                    int(row['HOST_TS']+row['overlap']),
                                                                                                    int(row['HOST_TE']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2_HOST_reversedCurr"]==16,
                                                                                                    row["R2_PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_host_y',
                                                                                                                                                                   'meanQual_host_y',
                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                   'meanQual_host_y']]

    
    dataR2Left["READ_LEN"]=dataR2Left["READ_LEN_R2"]
    dataR2Left["R"]="R2"
    dataR2Left.drop(dropList,axis=1,inplace=True)

    frames=[dataR1Right,dataR1Left,dataR2Right,dataR2Left]
    df=pd.concat(frames).reset_index(drop=True)
    global reportDF
    dataLen=df["SEQ"].str.len()
    reportDF["number of reads after minimum alignment length cutoff"]=len(df)
    reportDF["read len mean after minimum alignment length cutoff"]=dataLen.mean()
    reportDF["read len std after minimum alignment length cutoff"]=dataLen.std()
    reportDF["read len min after minimum alignment length cutoff"]=dataLen.min()
    reportDF["read len max after minimum alignment length cutoff"]=dataLen.max()
    if not args.overlap<0 and not args.gap<0:
        df=df[(df["overlap"]<=args.overlap)&(df["gap"]<=args.gap)]
    elif not args.overlap<0:
        df=df[df["overlap"]<=args.overlap]
    elif not args.gap<0:
        df=df[df["gap"]<=args.gap]
    else:
        df=df

    # global reportDF
    dataLen=df["SEQ"].str.len()
    reportDF["number of reads after overlap/gap cutoff"]=len(df)
    reportDF["read len mean after overlap/gap cutoff"]=dataLen.mean()
    reportDF["read len std after overlap/gap cutoff"]=dataLen.std()
    reportDF["read len min after overlap/gap cutoff"]=dataLen.min()
    reportDF["read len max after overlap/gap cutoff"]=dataLen.max()

    k=float(args.minLen)
    ssAl=args.steepSlopeAL

    df["PATHOGEN_AL_score"]=((df['PATHOGEN_AL']-k)/((ssAl+(df['PATHOGEN_AL']-k)**2.0)**0.5)+1.0)/2.0
    df["HOST_AL_score"]=((df['HOST_AL']-k)/((ssAl+(df['HOST_AL']-k)**2.0)**0.5)+1.0)/2.0

    df['jointEntropy']=((df['entropyScore_pathogen'] \
                    +df['entropyScore_host']) \
                    /(2))
    df['jointAlLen']=((df['HOST_AL_score'] \
                    +df['PATHOGEN_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.HOST,df.PATHOGEN,df.overlap,df.gap))))

    # global reportDF
    dataLen=df["SEQ"].str.len()
    reportDF["number of reads after intermediate score cutoff"]=len(df)
    reportDF["read len mean after intermediate score cutoff"]=dataLen.mean()
    reportDF["read len std after intermediate score cutoff"]=dataLen.std()
    reportDF["read len min after intermediate score cutoff"]=dataLen.min()
    reportDF["read len max after intermediate score cutoff"]=dataLen.max()
    
    return df

def filterOverlapCombineUnpaired(data,args):
    #right
    data['entropyScore_pathogen']=0
    data['meanQual_pathogen']=0
    data['mapQual_pathogen']=0
    data['entropyScore_host']=0
    data['meanQual_host']=0
    data['mapQual_host']=0
    dataRight=data[data["PATHOGEN"].str.contains("r")]
    dataRight=dataRight[~((dataRight["HOST_TS"]<dataRight["PATHOGEN_TS"])&(dataRight["HOST_TE"]>dataRight["PATHOGEN_TE"]))]
    dataRight=dataRight[~((dataRight["PATHOGEN_TS"]<dataRight["HOST_TS"])&(dataRight["PATHOGEN_TE"]>dataRight["HOST_TE"]))]
    dataRight["ins"]=dataRight["PATHOGEN_TS"]-dataRight["HOST_TE"]
    dataRight["split"]=dataRight['HOST_RE'].astype(int).astype(str)+":"+dataRight['PATHOGEN_RS'].astype(int).astype(str)
    dataRight['HOST']=dataRight['HOST_RE'].astype(int)
    dataRight["comb"]=dataRight.split+"@"+dataRight.HOST_ID+":"+dataRight.PATHOGEN_ID
    dataRight["orient"]="host:pathogen"
    dataRight["overlap"]=dataRight["overlap"]
    dataRight["gap"]=dataRight["gap"]
    dataRight["HOST_TS"]=dataRight["HOST_TS"].astype(int)
    dataRight["HOST_TE"]=dataRight["HOST_TE"].astype(int)
    dataRight["HOST_RS"]=dataRight["HOST_RS"].astype(int)
    dataRight["HOST_RE"]=dataRight["HOST_RE"].astype(int)
    dataRight["PATHOGEN_TS"]=dataRight["PATHOGEN_TS"].astype(int)
    dataRight["PATHOGEN_TE"]=dataRight["PATHOGEN_TE"].astype(int)
    dataRight["PATHOGEN_RS"]=dataRight["PATHOGEN_RS"].astype(int)
    dataRight["PATHOGEN_RE"]=dataRight["PATHOGEN_RE"].astype(int)
    dataRight["PATHOGEN_MAPQ"]=dataRight["PATHOGEN_MAPQ"].astype(int)
    dataRight["HOST_MAPQ"]=dataRight["HOST_MAPQ"].astype(int)
    dataRight["SEQ"]=dataRight["HOST_SEQ"]
    dataRight["PATHOGEN_AL"]=dataRight["PATHOGEN_TE"]-dataRight["PATHOGEN_TS"]-dataRight["overlap"]
    dataRight["HOST_AL"]=dataRight["HOST_TE"]-dataRight["HOST_TS"]-dataRight["overlap"]
    dataRight=dataRight[(dataRight['PATHOGEN_AL']>args.minLen)&(dataRight['HOST_AL']>args.minLen)]
    if len(dataRight)>0:
        dataRight[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataRight.merge(dataRight.apply(lambda row: processAligns(row['HOST_SEQ'],
                                                                                                    row['PATHOGEN_SEQ'],
                                                                                                    row['QUAL'],
                                                                                                    int(row['PATHOGEN_TS']+row['overlap']),
                                                                                                    int(row['PATHOGEN_TE']),
                                                                                                    int(row['HOST_TS']),
                                                                                                    int(row['HOST_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN']),
                                                                                                    row["HOST_reversedCurr"]==16,
                                                                                                    row["PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                                                   'meanQual_pathogen_y',
                                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                                   'meanQual_host_y']]
    dataRight["READ_LEN"]=dataRight["READ_LEN"]
    dataRight["R"]=0

    # left
    dataLeft=data[data["PATHOGEN"].str.contains("l")]
    dataLeft=dataLeft[~((dataLeft["HOST_TS"]<dataLeft["PATHOGEN_TS"])&(dataLeft["HOST_TE"]>dataLeft["PATHOGEN_TE"]))]
    dataLeft=dataLeft[~((dataLeft["PATHOGEN_TS"]<dataLeft["HOST_TS"])&(dataLeft["PATHOGEN_TE"]>dataLeft["HOST_TE"]))]
    dataLeft["ins"]=dataLeft["HOST_TS"]-dataLeft["PATHOGEN_TE"]
    dataLeft["split"]=dataLeft['PATHOGEN_RE'].astype(int).astype(str)+":"+dataLeft['HOST_RS'].astype(int).astype(str)
    dataLeft['HOST']=dataLeft['HOST_RS'].astype(int)
    dataLeft["comb"]=dataLeft.split+"@"+dataLeft.PATHOGEN_ID+":"+dataLeft.HOST_ID
    dataLeft["orient"]="pathogen:host"
    dataLeft["overlap"]=dataLeft["overlap"]
    dataLeft["gap"]=dataLeft["gap"]
    dataLeft["HOST_TS"]=dataLeft["HOST_TS"].astype(int)
    dataLeft["HOST_TE"]=dataLeft["HOST_TE"].astype(int)
    dataLeft["HOST_RS"]=dataLeft["HOST_RS"].astype(int)
    dataLeft["HOST_RE"]=dataLeft["HOST_RE"].astype(int)
    dataLeft["PATHOGEN_TS"]=dataLeft["PATHOGEN_TS"].astype(int)
    dataLeft["PATHOGEN_TE"]=dataLeft["PATHOGEN_TE"].astype(int)
    dataLeft["PATHOGEN_RS"]=dataLeft["PATHOGEN_RS"].astype(int)
    dataLeft["PATHOGEN_RE"]=dataLeft["PATHOGEN_RE"].astype(int)
    dataLeft["PATHOGEN_MAPQ"]=dataLeft["PATHOGEN_MAPQ"].astype(int)
    dataLeft["HOST_MAPQ"]=dataLeft["HOST_MAPQ"].astype(int)
    dataLeft["SEQ"]=dataLeft["HOST_SEQ"]
    dataLeft["PATHOGEN_AL"]=dataLeft["PATHOGEN_TE"]-dataLeft["PATHOGEN_TS"]-dataLeft["overlap"]
    dataLeft["HOST_AL"]=dataLeft["HOST_TE"]-dataLeft["HOST_TS"]-dataLeft["overlap"]
    dataLeft=dataLeft[(dataLeft["PATHOGEN_AL"]>args.minLen)&(dataLeft["HOST_AL"]>args.minLen)]
    if len(dataLeft)>0:
        dataLeft[['entropyScore_pathogen',
                     'meanQual_pathogen',
                     'entropyScore_host',
                     'meanQual_host']]=dataLeft.merge(dataLeft.apply(lambda row: processAligns(row['HOST_SEQ'],
                                                                                            row['PATHOGEN_SEQ'],
                                                                                            row['QUAL'],
                                                                                            int(row['PATHOGEN_TS']),
                                                                                            int(row['PATHOGEN_TE']-row['overlap']),
                                                                                            int(row['HOST_TS']+row['overlap']),
                                                                                            int(row['HOST_TE']),
                                                                                            int(row['READ_LEN']),
                                                                                            row["HOST_reversedCurr"]==16,
                                                                                            row["PATHOGEN_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                                   'meanQual_pathogen_y',
                                                                                                                                                                   'entropyScore_host_y',
                                                                                                                                                                   'meanQual_host_y']]

    dataLeft["READ_LEN"]=dataLeft["READ_LEN"]
    dataLeft["R"]=0

    frames=[dataRight,dataLeft]
    df=pd.concat(frames).reset_index(drop=True)
    global reportDF
    dataLen=df["SEQ"].str.len()
    reportDF["number of reads after minimum alignment length cutoff"]=len(df)
    reportDF["read len mean after minimum alignment length cutoff"]=dataLen.mean()
    reportDF["read len std after minimum alignment length cutoff"]=dataLen.std()
    reportDF["read len min after minimum alignment length cutoff"]=dataLen.min()
    reportDF["read len max after minimum alignment length cutoff"]=dataLen.max()
    if not args.overlap<0 and not args.gap<0:
        df=df[(df["overlap"]<=args.overlap)&(df["gap"]<=args.gap)]
    elif not args.overlap<0:
        df=df[df["overlap"]<=args.overlap]
    elif not args.gap<0:
        df=df[df["gap"]<=args.gap]
    else:
        df=df

    # global reportDF
    dataLen=df["SEQ"].str.len()
    reportDF["number of reads after minimum overlap/gap cutoff"]=len(df)
    reportDF["read len mean after minimum overlap/gap cutoff"]=dataLen.mean()
    reportDF["read len std after minimum overlap/gap cutoff"]=dataLen.std()
    reportDF["read len min after minimum overlap/gap cutoff"]=dataLen.min()
    reportDF["read len max after minimum overlap/gap cutoff"]=dataLen.max()

    k=float(args.minLen)
    ssAl=args.steepSlopeAL

    df["PATHOGEN_AL_score"]=((df['PATHOGEN_AL']-k)/((ssAl+(df['PATHOGEN_AL']-k)**2.0)**0.5)+1.0)/2.0
    df["HOST_AL_score"]=((df['HOST_AL']-k)/((ssAl+(df['HOST_AL']-k)**2.0)**0.5)+1.0)/2.0

    df['jointEntropy']=((df['entropyScore_pathogen'] \
                    +df['entropyScore_host']) \
                    /(2))
    df['jointAlLen']=((df['HOST_AL_score'] \
                    +df['PATHOGEN_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.HOST,df.PATHOGEN,df.overlap,df.gap))))

    return df

def addSpan(data,dataPos):
    def testR1(row,dataPathogen):
        dataPathogenR1=dataPathogen[(dataPathogen['R1_PATHOGEN_RS']-row['PATHOGEN_RS']>-500)&(dataPathogen['R1_PATHOGEN_RS']-row['PATHOGEN_RS']<0)] # check that the start of the pathogen is before the start of the pathogen in the grouped dataframe
        dataPathogenR1=dataPathogenR1[(dataPathogenR1['R2_HOST_RE']-row['HOST_RE']<500)&(dataPathogenR1['R2_HOST_RE']-row['HOST_RE']>0)]
        if len(dataPathogenR1)>0:
            return [set(dataPathogenR1["QNAME"]),len(dataPathogenR1)] # return set of reads and count of reads
        return [{''},0]
        
    def testR2(row,dataPathogen):
        dataPathogenR2=dataPathogen[(dataPathogen['R2_PATHOGEN_RE']-row['PATHOGEN_RE']<500)&(dataPathogen['R2_PATHOGEN_RE']-row['PATHOGEN_RE']>0)] # check that the start of the pathogen is before the start of the pathogen in the grouped dataframe
        dataPathogenR2=dataPathogenR2[(dataPathogenR2['R1_HOST_RS']-row['HOST_RS']>-500)&(dataPathogenR2['R1_HOST_RS']-row['HOST_RS']<0)]
        if len(dataPathogenR2)>0:
            return [set(dataPathogenR2["QNAME"]),len(dataPathogenR2)] # return set of reads and count of reads
        return [{''},0]

    datapathogenR1=data[data['PATHOGEN'].str.contains('sepR2:1')]
    datapathogenR1=datapathogenR1[~(datapathogenR1['R1_PATHOGEN_ID']==0)] # check that the pathogen is on the same side
    dataPosR1=dataPos[dataPos['orient'].str.contains('pathogen:host')]
    if len(dataPosR1)>0:
        dataPosR1[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR1.apply(lambda row: testR1(row,datapathogenR1),axis=1)])

    dataPathogenR2=data[data['PATHOGEN'].str.contains('sepR1:2')]
    dataPathogenR2=dataPathogenR2[~(dataPathogenR2['R2_PATHOGEN_ID']==0)] # check that the pathogen is on the same side
    dataPosR2=dataPos[dataPos['orient'].str.contains('host:pathogen')]
    if len(dataPosR2)>0:
        dataPosR2[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR2.apply(lambda row: testR2(row,dataPathogenR2),axis=1)])

    frames=[dataPosR1,dataPosR2]
    df=pd.concat(frames).reset_index(drop=True)
    return df

# this function should do the following:
# 1. group data by split position
# 2. add a column of high confidence support reads from the first parameter DF
# 3. add a column of low confidence support reads from the second parameter DF
def findSupport(data,minLen,unpaired):
    dataPos=pd.DataFrame(data.groupby(by=["comb","split","HOST_ID","R","orient"])[["QNAME",
                                                                                    "HOST_RS",
                                                                                    "HOST_RE",
                                                                                    "PATHOGEN_RS",
                                                                                    "PATHOGEN_RE",
                                                                                    "HOST_AL",
                                                                                    "PATHOGEN_AL",
                                                                                    "READ_LEN",
                                                                                    "entropyScore_pathogen",
                                                                                    "meanQual_pathogen",
                                                                                    "entropyScore_host",
                                                                                    "meanQual_host",
                                                                                    "PATHOGEN_MAPQ",
                                                                                    "HOST_MAPQ",
                                                                                    "SEQ"]].agg(count=("QNAME", "count"),
                                                                                                reads=("QNAME", lambda x: set(x)),
                                                                                                seq=("SEQ", lambda x: (max(dict(list(x)), key=float),dict(list(x))[max(dict(list(x)), key=float)])),
                                                                                                HOST_RS=('HOST_RS', 'min'),
                                                                                                HOST_RE=('HOST_RE', 'max'),
                                                                                                PATHOGEN_RS=('PATHOGEN_RS','min'),
                                                                                                PATHOGEN_RE=('PATHOGEN_RE','max'),
                                                                                                PATHOGEN_AL=('PATHOGEN_AL','sum'),
                                                                                                HOST_AL=('HOST_AL','sum'),
                                                                                                READ_LEN=('READ_LEN','sum'),
                                                                                                entropyScore_pathogen=('entropyScore_pathogen','sum'),
                                                                                                entropyScore_host=('entropyScore_host','sum'),
                                                                                                meanQual_pathogen=('meanQual_pathogen','sum'),
                                                                                                meanQual_host=('meanQual_host','sum'),
                                                                                                PATHOGEN_MAPQ=('PATHOGEN_MAPQ','sum'),
                                                                                                HOST_MAPQ=('HOST_MAPQ','sum'))).reset_index()
    dataPos.rename(columns={'HOST_ID':'chr'}, inplace=True)            
    return dataPos

def score(dataPos,args):
    dataPos["READ_LEN"]=dataPos["READ_LEN"].astype(float)/dataPos["count"].astype(float)
    dataPos["PATHOGEN_MAPQ"]=dataPos["PATHOGEN_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["HOST_MAPQ"]=dataPos["HOST_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["PATHOGEN_AL"]=dataPos["PATHOGEN_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["HOST_AL"]=dataPos["HOST_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_pathogen"]=dataPos["entropyScore_pathogen"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_host"]=dataPos["entropyScore_host"].astype(float)/dataPos["count"].astype(float)
    dataPos["count"]=dataPos["count"].astype(float)

    dataPos["PATHOGEN_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["PATHOGEN_AL"]-dataPos["READ_LEN"]/2).abs())
    dataPos["HOST_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["HOST_AL"]-dataPos["READ_LEN"]/2).abs())

    k=float(args.minLen)
    ssAl=args.steepSlopeAL

    dataPos["PATHOGEN_AL_score"]=((dataPos['PATHOGEN_AL']-k)/((ssAl+(dataPos['PATHOGEN_AL']-k)**2.0)**0.5)+1.0)/2.0
    dataPos["HOST_AL_score"]=((dataPos['HOST_AL']-k)/((ssAl+(dataPos['HOST_AL']-k)**2.0)**0.5)+1.0)/2.0

    m=float(args.minCount)/2.0
    dataPos["count_score"]=((dataPos['count']-m)/((1.0+(dataPos['count']-m)**2.0)**0.5)+1.0)/2.0 \
                                *(1.0-args.maxCountPenalty) \
                                +args.maxCountPenalty # algebraic sigmoid function of read count score

    dataPos['jointEntropy']=((dataPos['entropyScore_pathogen'] \
                    +dataPos['entropyScore_host']) \
                    /(2))
    dataPos['jointAlLen']=((dataPos['HOST_AL_score'] \
                    +dataPos['PATHOGEN_AL_score']) \
                    /(2))
    dataPos["score"]=(dataPos['jointEntropy'] \
                    *dataPos['count_score'] \
                    *dataPos['jointAlLen'])

    dataPos.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    dataPos=dataPos.round({'score': 4})
    return dataPos

# the following function orders the integration events by distance and computes distance to the next site
def approxCloseness(data,args):
    data.sort_values(by=["HOST_RS","gene_name"],inplace=True)
    data["diff1"]=abs(data['HOST_RS']-data['HOST_RS'].shift(-1))
    data["diff2"]=abs(data['HOST_RS']-data['HOST_RS'].shift(1))
    data['t1']=data["diff1"]<args.close
    data['t2']=data["diff2"]<args.close
    data['t']=data['t1']|data['t2']

    l=list(data["t"])
    uid=0
    start=True
    newL=[]
    for x in range(len(l)):
        if start:
            newL.append(uid)
            start=False
        elif l[x]==l[x-1]:
            newL.append(uid)
        else:
            uid=uid+1
            newL.append(uid)
    data.reset_index(inplace=True,drop=True)
    data["uid"]=pd.Series(newL)

    return data.drop(['diff1','diff2','t1','t2','t'],axis=1).reset_index(drop=True)

def writeReadNames(outDir,row,fileNameR1,fileNameR2,baseName,dirPath):
    outD=os.path.abspath(outDir)
    tempF=outD+'/tmp/fq'
    stringReads=row['reads'].replace(';','|')

    outPath1_All="'"+outD+"/Positions/"+baseName+"/"+str(row['gene_name'].strip('\n')+"@"+row["R"]+"@"+str(row['comb'].split("@")[0]))+"_R1.all.fa'"
    outPath2_All="'"+outD+"/Positions/"+baseName+"/"+str(row['gene_name'].strip('\n')+"@"+row["R"]+"@"+str(row['comb'].split("@")[0]))+"_R2.all.fa'"
    cmdR1="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR1+" | seqtk seq -a - > "+outPath1_All
    cmdR2="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR2+" | seqtk seq -a - > "+outPath2_All
    os.system(cmdR1)
    os.system(cmdR2)

def writeReadNamesUnpaired(outDir,row,fileName,baseName,dirPath):
    outD=os.path.abspath(outDir)
    tempF=outD+'/tmp/fq'
    stringReads=row['reads'].replace(';','|')

    outPath="'"+outD+"/Positions/"+baseName+"/"+str(row['gene_name'].strip('\n')+"@"+str(row['comb'].split("@")[0]))+".fa'"
    cmd="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileName+" | seqtk seq -a - > "+outPath
    os.system(cmd)
    
# the following function is designed to combine the two dataframes together
def combineLocalFull(dataPos,dataPosFull,minLen):    
    data=pd.DataFrame([],columns=list(dataPosFull))

    setLocalPos=set(dataPos["comb"])
    setFullPos=set(dataPosFull["comb"])
    intersect=setFullPos.intersection(setLocalPos)
    diff=setFullPos.symmetric_difference(setLocalPos)

    if len(intersect)>0:
        for el in intersect:
            df2=pd.DataFrame(0, index=np.arange(1), columns=['orient',
                                                            'prim',
                                                            'chr',
                                                            'split',
                                                            'R',
                                                            'spanR1-R2',
                                                            'allReads',
                                                            'spanCount',
                                                            'comb',
                                                            'reads',
                                                            'count',
                                                            'score'])

            l1=set(list(dataPos[dataPos['comb']==el]['spanR1-R2'])[0])
            l2=set(list(dataPosFull[dataPosFull['comb']==el]['spanR1-R2'])[0])
            setL=l1.union(l2)
            l1a=set(list(dataPos[dataPos['comb']==el]['reads'])[0])
            l2a=set(list(dataPosFull[dataPosFull['comb']==el]['reads'])[0])
            setLa=l1a.union(l2a)
            df2[['comb',
                'orient',
                'prim',
                'chr',
                'split',
                'R',
                'spanR1-R2',
                'reads',
                'spanCount']]=[el,
                               dataPos[dataPos["comb"]==el]["orient"].iloc[0],
                               10,
                               dataPos[dataPos["comb"]==el]["chr"].iloc[0],
                               dataPos[dataPos["comb"]==el]["split"].iloc[0],
                               dataPos[dataPos["comb"]==el]["R"].iloc[0],
                               ";".join(list(setL)),
                               ";".join(list(setLa)),
                               len(setL)]

            l1=set(list(dataPos[dataPos['comb']==el]['reads'])[0])
            l2=set(list(dataPosFull[dataPosFull['comb']==el]['reads'])[0])
            setL=l1.union(l2)
            df2[['count','reads']]=[len(setL),";".join(list(setL))]
            df2['score']=list(dataPos[dataPos['comb']==el]['score'])[0]+list(dataPosFull[dataPosFull['comb']==el]['score'])[0]

            data=data.append(df2)
            
        data=data.reset_index(drop=True)

    dataLocalPosDiff=dataPos[(dataPos['comb'].isin(diff))]
    if len(dataLocalPosDiff)>0:
        dataLocalPosDiff["reads"]=dataLocalPosDiff['reads'].str.join(";")

        dataLocalPosDiff['spanR1-R2']=dataLocalPosDiff['spanR1-R2'].str.join(";")
    dataLocalPosDiff["prim"]=1

    dataFullPosDiff=dataPosFull[(dataPosFull['comb'].isin(diff))]
    if len(dataFullPosDiff)>0:
        dataFullPosDiff["reads"]=dataFullPosDiff['reads'].str.join(";")

        dataFullPosDiff['spanR1-R2']=dataFullPosDiff['spanR1-R2'].str.join(";")
    dataFullPosDiff["prim"]=0
    data=pd.concat([data,dataLocalPosDiff,dataFullPosDiff])
    data.reset_index(drop=True)

    data.replace("","-",inplace=True)
    return data

# this function will produce a dataframe with information grouped by the SpliceSites
# those reads that do not contain a valid spliceSite shall be discarded
def groupBySpliceSites(data):
    colNames=["orient"]
        
    data.drop(colNames,axis=1,inplace=True)
    
    dfg=pd.DataFrame(data.groupby(by=["gene_name",
                                      "chr",
                                      "R",
                                      "uid"])[["comb",
                                                "reads",
                                                "count",
                                                "spanCount",
                                                "spanR1-R2",
                                                "entropyScore_host",
                                                "entropyScore_pathogen",
                                                "PATHOGEN_AL",
                                                "HOST_AL",
                                                "READ_LEN",
                                                "PATHOGEN_MAPQ",
                                                "HOST_MAPQ",
                                                "seq"]].agg(groupsCount=('reads','count'),
                                                            reads=('reads',lambda x: ';'.join(set(x))),
                                                            seq=('seq',lambda x: dict(list(x))[max(dict(list(x)), key=float)]),
                                                            spanReads=('spanR1-R2',lambda x: ';'.join(list(set([el for el in x if not el=='-'])))),
                                                            spanCount=('spanCount','sum'),
                                                            count=('count','sum'),
                                                            comb=('comb',lambda x: ';'.join(set(x))),
                                                            entropyScore_host=('entropyScore_host','sum'),
                                                            entropyScore_pathogen=('entropyScore_pathogen','sum'),
                                                            PATHOGEN_AL=('PATHOGEN_AL','sum'),
                                                            HOST_AL=('HOST_AL','sum'),
                                                            READ_LEN=('READ_LEN','sum'),
                                                            PATHOGEN_MAPQ=('PATHOGEN_MAPQ','sum'),
                                                            HOST_MAPQ=('HOST_MAPQ','sum'))).reset_index()

    return dfg.reset_index(drop=True)

def groupBySpliceSitesUnpaired(data):
    colNames=["orient"]
        
    data.drop(colNames,axis=1,inplace=True)
    
    dfg=pd.DataFrame(data.groupby(by=["gene_name",
                                      "chr",
                                      "R",
                                      "uid"])[["comb",
                                                "reads",
                                                "count",
                                                "entropyScore_host",
                                                "entropyScore_pathogen",
                                                "PATHOGEN_AL",
                                                "HOST_AL",
                                                "READ_LEN",
                                                "PATHOGEN_MAPQ",
                                                "HOST_MAPQ",
                                                "seq"]].agg(groupsCount=("reads","count"),
                                                            reads=("reads",lambda x: ';'.join(set(x))),
                                                            seq=("seq",lambda x: dict(list(x))[max(dict(list(x)), key=float)]),
                                                            count=('count','sum'),
                                                            comb=("comb",lambda x: ';'.join(set(x))),
                                                            entropyScore_host=('entropyScore_host','sum'),
                                                            entropyScore_pathogen=('entropyScore_pathogen','sum'),
                                                            PATHOGEN_AL=('PATHOGEN_AL','sum'),
                                                            HOST_AL=('HOST_AL','sum'),
                                                            READ_LEN=('READ_LEN','sum'),
                                                            PATHOGEN_MAPQ=('PATHOGEN_MAPQ','sum'),
                                                            HOST_MAPQ=('HOST_MAPQ','sum'))).reset_index()

    return dfg.reset_index(drop=True)

def get_gene_name(attribute_str: str) -> str:
    gff = attribute_str.startswith("ID=") or attribute_str.startswith("Parent=")

    attrs = attribute_str.rstrip().rstrip(";").split(";")
    attrs = [x.strip() for x in attrs]
    attrs = [x.strip("\"") for x in attrs]
    attrs_dict = dict()
    sep = " \""
    if gff:
        sep = "="
    for at in attrs:
        k, v = at.split(sep)
        attrs_dict.setdefault(k.lower(), v)
        
    for attr in ["gene_name","gene","name","gene_id","id"]:
        if attr in attrs_dict:
            return attrs_dict[attr]
        
    return "-"

def annotate(dataBed,annPath,data):
    # extract gene names for the annotation
    ann_df = pd.read_csv(annPath,sep="\t",comment="#",names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
    ann_df = ann_df[ann_df["type"]=="gene"].reset_index(drop=True)
    assert len(ann_df)>0,"annotation file must have records with type \"gene\" present. You can add gene records to your file with gffread"
        
    ann_df["gene_name"] = ann_df.attributes.apply(lambda row: get_gene_name(row))
    ann_df = ann_df[~(ann_df["gene_name"]=="-")].reset_index(drop=True)

    annBed = ann_df[["seqid","start","end","gene_name","score","strand"]]

    # load both chimeric results and annotation into bedtools objects
    sites=BedTool.from_dataframe(dataBed)
    annotation=BedTool.from_dataframe(annBed)

    nearby=annotation.intersect(sites, wo=True)
    if nearby.count()==0:
        data['gene_name']="-"
        return data
    
    # read the results of the intersection into a dataframe
    df = pd.read_table(nearby.fn,names=["seqid","start","end","gene_name","score","strand","q_seqid","q_start","q_end","distance"])
    
    # select only the closest gene to each chimeric site
    df["q_coords"]=df["q_seqid"].astype(str)+":"+df["q_start"].astype(str)+":"+df['q_end'].astype(str)
    min_dist_idxs = df.groupby('q_coords', group_keys=False)['distance'].idxmin()
    finalBed = df.loc[min_dist_idxs]

    if len(finalBed)==0:
        data['gene_name']="-"
        return data

    finalBed = finalBed[["q_seqid","q_start","q_end","gene_name"]]
    finalBed.columns=['chr',
                    'HOST_RS',
                    'HOST_RE',
                    'gene_name']

    finalDF=pd.DataFrame(pd.merge(data,finalBed,on=['chr',
                                                    'HOST_RS',
                                                    'HOST_RE'],how='left')).reset_index(drop=True)
    finalDF["gene_name"] = np.where(finalDF["gene_name"].isna(),"-",finalDF["gene_name"])
    return finalDF.reset_index(drop=True)

def rest(dataPos,args,data,unpaired,baseName,outDir,dirPath,mate):
    if not unpaired:
        dataPos=addSpan(data,dataPos)
        dataPos.loc[dataPos['spanCount'].isnull(),['spanCount']]=dataPos.loc[dataPos['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
        dataPos.loc[dataPos['spanR1-R2'].isnull(),['spanR1-R2']]=dataPos.loc[dataPos['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())

    dataBed=dataPos[['chr','HOST_RS','HOST_RE']].drop_duplicates()
    if len(dataPos)>0:
        dataPos=annotate(dataBed,os.path.abspath(args.annotation),dataPos)
        if not dataPos is None:
            dataPos=approxCloseness(dataPos,args)
            dataPos["reads"]=dataPos["reads"].str.join(";")
            if not unpaired:
                dataPos["spanR1-R2"]=dataPos["spanR1-R2"].str.join(";")
                dataPos["spanR1-R2"].replace("","-",inplace=True)

            if unpaired:
                dataPos=groupBySpliceSitesUnpaired(dataPos)
            else:
                dataPos=groupBySpliceSites(dataPos)

            dataPos=score(dataPos,args)
            dataPos=dataPos.sort_values(by='score',ascending=False).reset_index(drop=True)
            dataPos[['seq','host_pos','drop','overlap','gap']]=dataPos['seq'].apply(pd.Series)
            dataPos.drop("drop",axis=1,inplace=True)

            dataPos["sample"] = args.pathogenR1.split("/")[-1].split(".")[0]
    
            dataPos.to_csv(os.path.abspath(args.out)+"."+mate+".full.csv",index=False)
            dataPosClean=dataPos[(dataPos['entropyScore_pathogen']>args.minEntropy) \
                                &(dataPos['entropyScore_host']>args.minEntropy) \
                                &(dataPos['score']>args.score)]

            colsOrder=["sample",
                       "gene_name",
                       "chr",
                       "host_pos",
                       "R",
                       "seq",
                       "count",
                       "score",
                       "overlap",
                       "gap",
                       "fileName"]
            if unpaired:
                colsOrder.remove("R")
            if args.writeReads is None:
                colsOrder.remove("fileName")

            if unpaired:
                dataPosClean["fileName"]="tmp"#dataPosClean['gene_name'].str.strip('\n')+"@"+dataPosClean['comb'].str.split("@",expand=True)[0]+".fa"
            else:
                dataPosClean["fileName"]="tmp"#dataPosClean['gene_name'].str.strip('\n')+"@"+dataPosClean["R"]+"@"+dataPosClean['comb'].str.split("@",expand=True)[0]+"_R1.all.fa"

            dataPosClean[colsOrder].to_csv(os.path.abspath(args.out)+"."+mate+".csv",index=False)

            # completely forgot that we can write reads from the sam file
            # should make it much faster
            # for each comb in the Pos.csv file
            # find corresponding positions in the original host file
            # save to csv with \n delimeter
            if args.writeReads:
                if not unpaired:
                    fileNameR1=baseName+"_R1.fastq"
                    fileNameR2=baseName+"_R2.fastq"
                    dataPosClean=dataPosClean[dataPosClean['score']>0.9]
                    dataPosClean.apply(lambda row: writeReadNames(os.path.abspath(outDir),row,fileNameR1,fileNameR2,baseName,dirPath),axis=1)
                else:
                    fileName=baseName+".fastq"
                    dataPosClean=dataPosClean[dataPosClean['score']>0.9]
                    dataPosClean.apply(lambda row: writeReadNamesUnpaired(os.path.abspath(outDir),row,fileName,baseName,dirPath),axis=1)

def wrapperSpan(outDir,baseName,dirPath,fileName,minLen,args): 

    def filterReads(dataHost_R1,dataPathogen_R1,dataHost_R2,dataPathogen_R2):
        #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
        dataHost_R1=dataHost_R1[(dataHost_R1["secondaryAlignment"]==0)&(dataHost_R1["PCRdup"]==0)&(dataHost_R1["suppAl"]==0)&(dataHost_R1["noPassFilter"]==0)]
        dataPathogen_R1=dataPathogen_R1[(dataPathogen_R1["secondaryAlignment"]==0)&(dataPathogen_R1["PCRdup"]==0)&(dataPathogen_R1["suppAl"]==0)&(dataPathogen_R1["noPassFilter"]==0)]
        dataHost_R2=dataHost_R2[(dataHost_R2["secondaryAlignment"]==0)&(dataHost_R2["PCRdup"]==0)&(dataHost_R2["suppAl"]==0)&(dataHost_R2["noPassFilter"]==0)]
        dataPathogen_R2=dataPathogen_R2[(dataPathogen_R2["secondaryAlignment"]==0)&(dataPathogen_R2["PCRdup"]==0)&(dataPathogen_R2["suppAl"]==0)&(dataPathogen_R2["noPassFilter"]==0)]
        return dataHost_R1.reset_index(drop=True),dataPathogen_R1.reset_index(drop=True),dataHost_R2.reset_index(drop=True),dataPathogen_R2.reset_index(drop=True)

    def createData(data,dataHost_R1,dataPathogen_R1,dataHost_R2,dataPathogen_R2):
        dataHostR1_PathogenR2=data.merge(dataHost_R1,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
        dataHostR1_PathogenR2=dataHostR1_PathogenR2.merge(dataPathogen_R2,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
        
        dataHostR2_PathogenR1=data.merge(dataHost_R2,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
        dataHostR2_PathogenR1=dataHostR2_PathogenR1.merge(dataPathogen_R1,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
        
        return data

    def processAligns(seqPathogen,seqHost,qual_pathogen,qual_host,i1_pathogen,i2_pathogen,i1_host,i2_host,len_pathogen,len_host,rPathogen,rHost):

        entropyScore_pathogen=0
        entropyScore_host=0
        meanQual_pathogen=0
        meanQual_host=0

        s_pathogen=""
        if rPathogen==True: # if reverse complemented take into account
            s_pathogen=seqPathogen[len_pathogen-i2_pathogen:len_pathogen-i1_pathogen]
        else:
            s_pathogen=seqPathogen[i1_pathogen:i2_pathogen]
        if not len(s_pathogen)==0:
            entropyScore_pathogen=topologicalNormalizedEntropy(s_pathogen)
            q_pathogen=qual_pathogen[i1_pathogen:i2_pathogen]
            if len(q_pathogen)>0:
                meanQual_pathogen=sum([ord(x)-33 for x in q_pathogen])/len(q_pathogen)
            else:
                meanQual_pathogen=0

        s_host=""
        if rHost==True: # if reverse complemented take into account
            s_host=seqHost[len_host-i2_host:len_host-i1_host]
        else:
            s_host=seqHost[i1_host:i2_host]
        if not len(s_host)==0:
            entropyScore_host=topologicalNormalizedEntropy(s_host)
            q_host=qual_host[i1_host:i2_host]
            if len(q_host)>0:
                meanQual_host=sum([ord(x)-33 for x in q_host])/len(q_host)
            else:
                meanQual_host=0

        return pd.Series({"entropyScore_pathogen":entropyScore_pathogen,
                            "meanQual_pathogen":meanQual_pathogen,
                            "entropyScore_host":entropyScore_host,
                            "meanQual_host":meanQual_host})

    def filterOverlapCombine(data):
        #right
        df=data.copy(deep=True)
        df['entropyScore_pathogen']=0
        df['meanQual_pathogen']=0
        df['mapQual_pathogen']=0
        df['entropyScore_host']=0
        df['meanQual_host']=0
        df['mapQual_host']=0
        df['HOST']=df['HOST_RE'].astype(int)
        df['PATHOGEN']=df['PATHOGEN_RE'].astype(int)
        df["comb"]=df.HOST.astype(str)+"@"+df.HOST_ID+":"+df.PATHOGEN_ID
        df["HOST_TS"]=df["HOST_TS"].astype(int)
        df["HOST_TE"]=df["HOST_TE"].astype(int)
        df["HOST_RS"]=df["HOST_RS"].astype(int)
        df["HOST_RE"]=df["HOST_RE"].astype(int)
        df["HOST_MAPQ"]=df["HOST_MAPQ"].astype(int)
        df["PATHOGEN_TS"]=df["PATHOGEN_TS"].astype(int)
        df["PATHOGEN_TE"]=df["PATHOGEN_TE"].astype(int)
        df["PATHOGEN_RS"]=df["PATHOGEN_RS"].astype(int)
        df["PATHOGEN_RE"]=df["PATHOGEN_RE"].astype(int)
        df["PATHOGEN_MAPQ"]=df["PATHOGEN_MAPQ"].astype(int)
        df["PATHOGEN_AL"]=df["PATHOGEN_TE"]-df["PATHOGEN_TS"]
        df["HOST_AL"]=df["HOST_TE"]-df["HOST_TS"]
        df=df[(df['PATHOGEN_AL']>minLen)&((df["PATHOGEN_LEN"]-df['PATHOGEN_AL'])<args.maxLenUnmapped)&(df['HOST_AL']>minLen)&((df["HOST_LEN"]-df['HOST_AL'])<args.maxLenUnmapped)].reset_index(drop=True)
        if len(df)>0:
            df[['entropyScore_pathogen',
                  'meanQual_pathogen',
                  'entropyScore_host',
                  'meanQual_host']]=df.merge(df.apply(lambda row: processAligns(row['PATHOGEN_SEQ'],
                                                                                    row['HOST_SEQ'],
                                                                                    row['PATHOGEN_QUAL'],
                                                                                    row['HOST_QUAL'],
                                                                                    int(row['PATHOGEN_TS']),
                                                                                    int(row['PATHOGEN_TE']),
                                                                                    int(row['HOST_TS']),
                                                                                    int(row['HOST_TE']),
                                                                                    int(row['PATHOGEN_LEN']),
                                                                                    int(row['HOST_LEN']),
                                                                                    row["PATHOGEN_reversedCurr"]==16,
                                                                                    row["HOST_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_pathogen_y',
                                                                                                                                                            'meanQual_pathogen_y',
                                                                                                                                                            'entropyScore_host_y',
                                                                                                                                                            'meanQual_host_y']]
            df = df.reset_index(drop=True)


        k=float(args.minLen)
        ssAl=args.steepSlopeAL
        df["PATHOGEN_AL_score"]=((df['PATHOGEN_AL']-k)/((ssAl+(df['PATHOGEN_AL']-k)**2.0)**0.5)+1.0)/2.0
        df["HOST_AL_score"]=((df['HOST_AL']-k)/((ssAl+(df['HOST_AL']-k)**2.0)**0.5)+1.0)/2.0

        df['jointEntropy']=((df['entropyScore_pathogen'] \
                        +df['entropyScore_host']) \
                        /(2))
        df['jointAlLen']=((df['HOST_AL_score'] \
                        +df['PATHOGEN_AL_score']) \
                        /(2))
        df["scorePrelim"]=(df['jointEntropy'] \
                        *df['jointAlLen'])

        df.drop(["jointEntropy",
                      "jointAlLen"],axis=1,inplace=True)
        df=df.round({'scorePrelim': 4})
        df=df[df["scorePrelim"]>=args.score]
        df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.HOST_SEQ,df.PATHOGEN_SEQ,df.HOST_RS.astype(int),df.HOST_RE.astype(int),df.PATHOGEN_RS.astype(int),df.PATHOGEN_RE.astype(int)))))

        return df.reset_index(drop=True)

    def findSupport(data,minLen,unpaired):
        
        dataPos=pd.DataFrame(data.groupby(by=["comb","HOST_ID"])[["QNAME",
                                                                "HOST_RS",
                                                                "HOST_RE",
                                                                "PATHOGEN_RS",
                                                                "PATHOGEN_RE",
                                                                "HOST_AL",
                                                                "PATHOGEN_AL",
                                                                "HOST_LEN",
                                                                "PATHOGEN_LEN",
                                                                "entropyScore_pathogen",
                                                                "meanQual_pathogen",
                                                                "entropyScore_host",
                                                                "meanQual_host",
                                                                "PATHOGEN_MAPQ",
                                                                "HOST_MAPQ",
                                                                "SEQ"]].agg(count=('QNAME','count'),
                                                                            reads=('QNAME', lambda x: set(x)),
                                                                            seq=('SEQ',lambda x: (max(dict(list(x)), key=float),dict(list(x))[max(dict(list(x)), key=float)])),
                                                                            HOST_RS=('HOST_RS','min'),
                                                                            HOST_RE=('HOST_RE','max'),
                                                                            PATHOGEN_RS=('PATHOGEN_RS','min'),
                                                                            PATHOGEN_RE=('PATHOGEN_RE','max'),
                                                                            PATHOGEN_AL=('PATHOGEN_AL','sum'),
                                                                            HOST_AL=('HOST_AL','sum'),
                                                                            HOST_LEN=('HOST_LEN','sum'),
                                                                            PATHOGEN_LEN=('PATHOGEN_LEN','sum'),
                                                                            entropyScore_pathogen=('entropyScore_pathogen','sum'),
                                                                            entropyScore_host=('entropyScore_host','sum'),
                                                                            meanQual_pathogen=('meanQual_pathogen','sum'),
                                                                            meanQual_host=('meanQual_host','sum'),
                                                                            PATHOGEN_MAPQ=('PATHOGEN_MAPQ','sum'),
                                                                            HOST_MAPQ=('HOST_MAPQ','sum'))).reset_index()
        dataPos.rename(columns={'HOST_ID':'chr'}, inplace=True)            
        return dataPos

    def groupBySpliceSites(data):
        
        dfg=pd.DataFrame(data.groupby(by=["gene_name",
                                          "chr",
                                          "uid"])[["comb",
                                                    "reads",
                                                    "count",
                                                    "entropyScore_host",
                                                    "entropyScore_pathogen",
                                                    "PATHOGEN_AL",
                                                    "HOST_AL",
                                                    "HOST_LEN",
                                                    "PATHOGEN_LEN",
                                                    "PATHOGEN_MAPQ",
                                                    "HOST_MAPQ",
                                                    "seq"]].agg(groupsCount=('reads','count'),
                                                                reads=('reads',lambda x: ';'.join(set(x))),
                                                                seq=('seq',lambda x: dict(list(x))[max(dict(list(x)), key=float)]),
                                                                count=('count','sum'),
                                                                comb=('comb',lambda x: ';'.join(set(x))),
                                                                entropyScore_host=('entropyScore_host','sum'),
                                                                entropyScore_pathogen=('entropyScore_pathogen','sum'),
                                                                PATHOGEN_AL=('PATHOGEN_AL','sum'),
                                                                HOST_AL=('HOST_AL','sum'),
                                                                HOST_LEN=('HOST_LEN','sum'),
                                                                PATHOGEN_LEN=('PATHOGEN_LEN','sum'),
                                                                PATHOGEN_MAPQ=('PATHOGEN_MAPQ','sum'),
                                                                HOST_MAPQ=('HOST_MAPQ','sum'))).reset_index()

        return dfg.reset_index(drop=True)

    def score(dataPos,args):
        dataPos["HOST_LEN"]=dataPos["HOST_LEN"].astype(float)/dataPos["count"].astype(float)
        dataPos["PATHOGEN_LEN"]=dataPos["PATHOGEN_LEN"].astype(float)/dataPos["count"].astype(float)
        dataPos["PATHOGEN_MAPQ"]=dataPos["PATHOGEN_MAPQ"].astype(float)/dataPos["count"].astype(float)
        dataPos["HOST_MAPQ"]=dataPos["HOST_MAPQ"].astype(float)/dataPos["count"].astype(float)
        dataPos["PATHOGEN_AL"]=dataPos["PATHOGEN_AL"].astype(float)/dataPos["count"].astype(float)
        dataPos["HOST_AL"]=dataPos["HOST_AL"].astype(float)/dataPos["count"].astype(float)
        dataPos["entropyScore_pathogen"]=dataPos["entropyScore_pathogen"].astype(float)/dataPos["count"].astype(float)
        dataPos["entropyScore_host"]=dataPos["entropyScore_host"].astype(float)/dataPos["count"].astype(float)
        dataPos["count"]=dataPos["count"].astype(float)

        k=float(args.minLen)
        ssAl=args.steepSlopeAL

        dataPos["PATHOGEN_AL_score"]=((dataPos['PATHOGEN_AL']-k)/((ssAl+(dataPos['PATHOGEN_AL']-k)**2.0)**0.5)+1.0)/2.0
        dataPos["HOST_AL_score"]=((dataPos['HOST_AL']-k)/((ssAl+(dataPos['HOST_AL']-k)**2.0)**0.5)+1.0)/2.0

        m=float(args.minCount)/2.0
        dataPos["count_score"]=((dataPos['count']-m)/((1.0+(dataPos['count']-m)**2.0)**0.5)+1.0)/2.0 \
                                    *(1.0-args.maxCountPenalty) \
                                    +args.maxCountPenalty # algebraic sigmoid function of read count score

        dataPos['jointEntropy']=((dataPos['entropyScore_pathogen'] \
                        +dataPos['entropyScore_host']) \
                        /(2))
        dataPos['jointAlLen']=((dataPos['HOST_AL_score'] \
                        +dataPos['PATHOGEN_AL_score']) \
                        /(2))
        dataPos["score"]=(dataPos['jointEntropy'] \
                        *dataPos['count_score'] \
                        *dataPos['jointAlLen'])

        dataPos.drop(["jointEntropy",
                      "jointAlLen"],axis=1,inplace=True)
        dataPos=dataPos.round({'score': 4})
        return dataPos
    
    global sam_colnames
    dataPathogen_R1=pd.read_csv(os.path.abspath(args.pathogenR1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)
    dataHost_R1=pd.read_csv(os.path.abspath(args.hostR1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)
    dataPathogen_R2=pd.read_csv(os.path.abspath(args.pathogenR2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)
    dataHost_R2=pd.read_csv(os.path.abspath(args.hostR2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)

    if ((len(dataPathogen_R1)==0 and len(dataHost_R2)==0) or (len(dataPathogen_R2)==0 or len(dataHost_R1)==0)): #exit if either alignment is empty
        print("incorrect")
        
    extractFlagBits(dataPathogen_R1)
    extractFlagBits(dataHost_R1)
    extractFlagBits(dataPathogen_R2)
    extractFlagBits(dataHost_R2)

    dataHost_R1,dataPathogen_R1,dataHost_R2,dataPathogen_R2=filterReads(dataHost_R1,dataPathogen_R1,dataHost_R2,dataPathogen_R2)

    dataPathogen_R1=extractStartEnd(dataPathogen_R1)
    dataPathogen_R2=extractStartEnd(dataPathogen_R2)
    dataHost_R1=extractStartEnd(dataHost_R1)
    dataHost_R2=extractStartEnd(dataHost_R2)

    data=pd.DataFrame(pd.Series(list(set(dataPathogen_R1["QNAME"]).union(set(dataPathogen_R2["QNAME"]))))).reset_index(drop=True)
    data.columns=["QNAME"]

    dataHostR1_PathogenR2=data.merge(dataHost_R1,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
    dataHostR1_PathogenR2=dataHostR1_PathogenR2[dataHostR1_PathogenR2["_merge"]=="both"].drop("_merge",axis=1).reset_index(drop=True)
    dataHostR1_PathogenR2=dataHostR1_PathogenR2.merge(dataPathogen_R2,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
    dataHostR1_PathogenR2=dataHostR1_PathogenR2[dataHostR1_PathogenR2["_merge"]=="both"].drop("_merge",axis=1).reset_index(drop=True)

    dataHostR2_PathogenR1=data.merge(dataHost_R2,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
    dataHostR2_PathogenR1=dataHostR2_PathogenR1[dataHostR2_PathogenR1["_merge"]=="both"].drop("_merge",axis=1).reset_index(drop=True)
    dataHostR2_PathogenR1=dataHostR2_PathogenR1.merge(dataPathogen_R1,left_on='QNAME',right_on='QNAME',how='inner',indicator=True)
    dataHostR2_PathogenR1=dataHostR2_PathogenR1[dataHostR2_PathogenR1["_merge"]=="both"].drop("_merge",axis=1).reset_index(drop=True)

    colnames = ['QNAME','RNAME_x','MAPQ_x','QUAL_x','SEQ_x','reversedCurr_x','Template_start_x','Template_end_x',
                'Reference_start_x','Reference_end_x','READ_LEN_x','RNAME_y','MAPQ_y','QUAL_y','SEQ_y','reversedCurr_y',
                'Template_start_y','Template_end_y','Reference_start_y','Reference_end_y','READ_LEN_y']
    dataHostR1_PathogenR2=dataHostR1_PathogenR2[colnames]
    dataHostR2_PathogenR1=dataHostR2_PathogenR1[colnames]

    colnames = ['QNAME','HOST_ID','HOST_MAPQ','HOST_QUAL','HOST_SEQ','HOST_reversedCurr','HOST_TS','HOST_TE','HOST_RS','HOST_RE','HOST_LEN',
                'PATHOGEN_ID','PATHOGEN_MAPQ','PATHOGEN_QUAL','PATHOGEN_SEQ','PATHOGEN_reversedCurr','PATHOGEN_TS','PATHOGEN_TE','PATHOGEN_RS','PATHOGEN_RE','PATHOGEN_LEN']
    dataHostR1_PathogenR2.columns=colnames
    dataHostR2_PathogenR1.columns=colnames

    dataHostR1_PathogenR2.replace('', np.nan,inplace=True)
    dataHostR1_PathogenR2.fillna(0,inplace=True)
    dataHostR2_PathogenR1.replace('', np.nan,inplace=True)
    dataHostR2_PathogenR1.fillna(0,inplace=True)

    dataHostR1_PathogenR2=dataHostR1_PathogenR2[~(dataHostR1_PathogenR2["HOST_SEQ"]==dataHostR1_PathogenR2["PATHOGEN_SEQ"])]
    dataHostR2_PathogenR1=dataHostR2_PathogenR1[~(dataHostR2_PathogenR1["HOST_SEQ"]==dataHostR2_PathogenR1["PATHOGEN_SEQ"])]

    dataPosHostR1_PathogenR2=pd.DataFrame([])
    dataPosHostR1_PathogenR2=filterOverlapCombine(dataHostR1_PathogenR2)

    dataPosHostR2_PathogenR1=pd.DataFrame([])
    dataPosHostR2_PathogenR1=filterOverlapCombine(dataHostR2_PathogenR1)

    dataPosHostR1_PathogenR2=findSupport(dataPosHostR1_PathogenR2,args.minLen,False)
    dataPosHostR2_PathogenR1=findSupport(dataPosHostR2_PathogenR1,args.minLen,False)

    dataBedHostR1_PathogenR2=dataPosHostR1_PathogenR2[['chr','HOST_RS','HOST_RE']].drop_duplicates().reset_index(drop=True)
    dataBedHostR2_PathogenR1=dataPosHostR2_PathogenR1[['chr','HOST_RS','HOST_RE']].drop_duplicates().reset_index(drop=True)

    if len(dataPosHostR2_PathogenR1)>0:
        dataPosHostR2_PathogenR1=annotate(dataBedHostR2_PathogenR1,os.path.abspath(args.annotation),dataPosHostR2_PathogenR1)
    if len(dataPosHostR1_PathogenR2)>0:
        dataPosHostR1_PathogenR2=annotate(dataBedHostR1_PathogenR2,os.path.abspath(args.annotation),dataPosHostR1_PathogenR2)

    if not dataPosHostR1_PathogenR2 is None and len(dataPosHostR1_PathogenR2)>0:
        dataPosHostR1_PathogenR2=approxCloseness(dataPosHostR1_PathogenR2,args)
        dataPosHostR1_PathogenR2["reads"]=dataPosHostR1_PathogenR2["reads"].str.join(";")
    if not dataPosHostR2_PathogenR1 is None and len(dataPosHostR2_PathogenR1)>0:
        dataPosHostR2_PathogenR1=approxCloseness(dataPosHostR2_PathogenR1,args)
        dataPosHostR2_PathogenR1["reads"]=dataPosHostR2_PathogenR1["reads"].str.join(";")

    if not dataPosHostR1_PathogenR2 is None and len(dataPosHostR1_PathogenR2)>0:
        dataPosHostR1_PathogenR2=groupBySpliceSites(dataPosHostR1_PathogenR2)
    if not dataPosHostR2_PathogenR1 is None and len(dataPosHostR2_PathogenR1)>0:
        dataPosHostR2_PathogenR1=groupBySpliceSites(dataPosHostR2_PathogenR1)

    if not dataPosHostR1_PathogenR2 is None and len(dataPosHostR1_PathogenR2)>0:
        dataPosHostR1_PathogenR2=score(dataPosHostR1_PathogenR2,args)
    if not dataPosHostR2_PathogenR1 is None and len(dataPosHostR2_PathogenR1)>0:
        dataPosHostR2_PathogenR1=score(dataPosHostR2_PathogenR1,args)

    if not dataPosHostR1_PathogenR2 is None and len(dataPosHostR1_PathogenR2)>0:
        dataPosHostR1_PathogenR2=dataPosHostR1_PathogenR2.sort_values(by='score',ascending=False).reset_index(drop=True)
        dataPosHostR1_PathogenR2[['host_seq','pathogen_seq','host_rs','host_re','pathogen_rs','pathogen_re']]=dataPosHostR1_PathogenR2['seq'].apply(pd.Series)
        
    if not dataPosHostR2_PathogenR1 is None and len(dataPosHostR2_PathogenR1)>0:
        dataPosHostR2_PathogenR1=dataPosHostR2_PathogenR1.sort_values(by='score',ascending=False).reset_index(drop=True)
        dataPosHostR2_PathogenR1[['host_seq','pathogen_seq','host_rs','host_re','pathogen_rs','pathogen_re']]=dataPosHostR2_PathogenR1['seq'].apply(pd.Series)

    dataPos=pd.concat([pd.DataFrame(dataPosHostR1_PathogenR2),pd.DataFrame(dataPosHostR2_PathogenR1)])
    if not len(dataPos)>0:
        return
    
    dataPos["sample"] = args.pathogenR1.split("/")[-1].split(".")[0]
    colsOrder=["sample",
               "gene_name",
               "chr",
               "host_rs",
               "host_re",
               "pathogen_rs",
               "pathogen_re",
               "host_seq",
               "pathogen_seq",
               "count",
               "score"]

    dataPos.to_csv(os.path.abspath(args.out)+".span.full.csv",index=False)
    dataPos_Clean=dataPos[(dataPos['entropyScore_pathogen']>args.minEntropy) \
                            &(dataPos['entropyScore_host']>args.minEntropy) \
                            &(dataPos['score']>args.score)].reset_index(drop=True)

    if not dataPos_Clean is None and len(dataPos_Clean)>0:
        dataPos_Clean[colsOrder].to_csv(os.path.abspath(args.out)+".span.csv",index=False)

def wrapper(outDir,baseName,dirPath,fileName,minLen,args,pathogen_fname,host_fname,s,mate):
    # load data from local alignments
    global sam_colnames
    dataPathogen=pd.read_csv(os.path.abspath(pathogen_fname),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)
    dataHost=pd.read_csv(os.path.abspath(host_fname),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)

    if (len(dataPathogen)==0 or len(dataHost)==0): #exit if either alignment is empty
        return
    
    outDirPOS=outDir+"/Positions/"
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
        os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

    global reportDF
    if len(dataPathogen)>0 and len(dataHost)>0:
        dataLen=dataHost["SEQ"].str.len()
        reportDF["number of reads before removing spliced reads"]=len(dataHost)
        reportDF["original read len mean"]=dataLen.mean()
        reportDF["original read len std"]=dataLen.std()
        reportDF["original read len min"]=dataLen.min()
        reportDF["original read len max"]=dataLen.max()
        # extract flag information
        extractFlagBits(dataPathogen)
        extractFlagBits(dataHost)
        unpaired=False
        if len(dataPathogen[dataPathogen["paired"]>0])==0 and len(dataHost[dataHost["paired"]>0])==0:
            unpaired=True
        
        # need to verify that both paired and unpaired work below.
        # i guess for paired we need to consider the case when splicing occurs on one mate, and thus we only want to remove that mate
        # but not the other mate, since that might still contain a chimeric site, provided it is not an end-to-tnd alignment
        if s is not None:
            if os.path.exists(os.path.abspath(s)):
                dataSpliced=pd.read_csv(os.path.abspath(s),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=sam_colnames)
                extractFlagBits(dataSpliced)
                dataSpliced["tid"]=dataSpliced['QNAME']+dataSpliced['firstRead'].astype(str)+dataSpliced['lastRead'].astype(str)
                dataHost["tid"]=dataHost['QNAME']+dataHost['firstRead'].astype(str)+dataHost['lastRead'].astype(str)
                dataPathogen["tid"]=dataPathogen['QNAME']+dataPathogen['firstRead'].astype(str)+dataPathogen['lastRead'].astype(str)
                dataPathogen=dataPathogen[~dataPathogen['tid'].isin(set(dataSpliced['tid']))]
                dataHost=dataHost[dataHost['tid'].isin(set(dataPathogen['tid']))]
                dataPathogen=dataPathogen[dataPathogen['tid'].isin(set(dataHost['tid']))]
                dataPathogen.drop(['tid'],axis=1,inplace=True)
                dataHost.drop(['tid'],axis=1,inplace=True)

                dataLen=dataHost["SEQ"].str.len()
                reportDF["number of reads shared between pathogen and host post splicing"]=len(dataHost)
                reportDF["read len mean shared between pathogen and host post splicing"]=dataLen.mean()
                reportDF["read len std shared between pathogen and host post splicing"]=dataLen.std()
                reportDF["read len min shared between pathogen and host post splicing"]=dataLen.min()
                reportDF["read len max shared between pathogen and host post splicing"]=dataLen.max()
            else:
                print("Spliced file does not exist")
        elif unpaired==False: # still remove all host which are not in pathogen
            dataHost["tid"]=dataHost['QNAME']+dataHost['firstRead'].astype(str)+dataHost['lastRead'].astype(str)
            dataPathogen["tid"]=dataPathogen['QNAME']+dataPathogen['firstRead'].astype(str)+dataPathogen['lastRead'].astype(str)
            dataPathogen=dataPathogen[dataPathogen['tid'].isin(set(dataHost['tid']))]
            dataHost=dataHost[dataHost['tid'].isin(set(dataPathogen['tid']))]
            dataPathogen.drop(['tid'],axis=1,inplace=True)
            dataHost.drop(['tid'],axis=1,inplace=True)
            dataLen=dataHost["SEQ"].str.len()
            reportDF["number of reads shared between pathogen and host post splicing"]=len(dataHost)
            reportDF["read len mean shared between pathogen and host post splicing"]=dataLen.mean()
            reportDF["read len std shared between pathogen and host post splicing"]=dataLen.std()
            reportDF["read len min shared between pathogen and host post splicing"]=dataLen.min()
            reportDF["read len max shared between pathogen and host post splicing"]=dataLen.max()
        else:
            dataHost["tid"]=dataHost['QNAME']
            dataPathogen["tid"]=dataPathogen['QNAME']
            dataPathogen=dataPathogen[dataPathogen['tid'].isin(set(dataHost['tid']))]
            dataHost=dataHost[dataHost['tid'].isin(set(dataPathogen['tid']))]
            dataPathogen.drop(['tid'],axis=1,inplace=True)
            dataHost.drop(['tid'],axis=1,inplace=True)
            dataLen=dataHost["SEQ"].str.len()
            reportDF["number of reads shared between pathogen and host post splicing"]=len(dataHost)
            reportDF["read len mean shared between pathogen and host post splicing"]=dataLen.mean()
            reportDF["read len std shared between pathogen and host post splicing"]=dataLen.std()
            reportDF["read len min shared between pathogen and host post splicing"]=dataLen.min()
            reportDF["read len max shared between pathogen and host post splicing"]=dataLen.max()

        # extract start and end for both template and reference
        dataPathogen=extractStartEnd(dataPathogen)
        dataHost=extractStartEnd(dataHost)

        dataHost,dataPathogen=filterReads(dataHost,dataPathogen)
        dataLen=dataHost["SEQ"].str.len()
        reportDF["number of reads after sam flag filtering"]=len(dataHost)
        reportDF["read len mean after sam flag filtering"]=dataLen.mean()
        reportDF["read len std after sam flag filtering"]=dataLen.std()
        reportDF["read len min after sam flag filtering"]=dataLen.min()
        reportDF["read len max after sam flag filtering"]=dataLen.max()
        if len(dataHost)==0:
            return

        data=pd.DataFrame(dataPathogen["QNAME"]).reset_index(drop=True)
        dataHost=dataHost.reset_index(drop=True)
        dataPathogen=dataPathogen.reset_index(drop=True)
        if unpaired:
            data=createDataUnpaired(data,dataHost,dataPathogen)
        else:
            data=createData(data,dataHost,dataPathogen)

        del dataHost
        del dataPathogen
        data.replace('', np.nan,inplace=True)
        data.fillna(0,inplace=True)
        if unpaired:
            data=leftRightUnpaired(data,minLen)
            data=overlapUnpaired(data)
        else:
            data=leftRight(data,minLen)
            data=overlap(data)
        data.drop_duplicates(inplace=True)

        dataPos=pd.DataFrame([])
        if unpaired:
            dataPos=filterOverlapCombineUnpaired(data,args)
        else:
            dataPos=filterOverlapCombine(data,args)
        dataPos=dataPos[(dataPos["entropyScore_host"]>args.minEntropy)&(dataPos["entropyScore_pathogen"]>args.minEntropy)]

        dataPos=findSupport(dataPos,minLen,unpaired)
        if len(dataPos)>0:
            rest(dataPos,args,data,unpaired,baseName,outDir,dirPath,mate)
    
    return 1

def main(args):
    outDir="/".join(os.path.abspath(args.out).split("/")[:-1])
    if not os.path.exists(outDir):
        print('output directory does not exist')
        return

    pathogen_r1=os.path.abspath(args.pathogenR1)
    host_r1=os.path.abspath(args.hostR1)
    pathogen_r2=os.path.abspath(args.pathogenR2)
    host_r2=os.path.abspath(args.hostR2)
    fullPath1r1=os.path.abspath(os.path.realpath(pathogen_r1))
    fullPath2r1=os.path.abspath(os.path.realpath(host_r1))
    fullPath1r2=os.path.abspath(os.path.realpath(pathogen_r2))
    fullPath2r2=os.path.abspath(os.path.realpath(host_r2))

    if os.path.exists(fullPath1r1) and os.path.exists(fullPath2r1) and os.path.exists(fullPath1r2) and os.path.exists(fullPath2r2):
        fileName=fullPath1r1.split('/')[-1]
        dirPath="/".join(fullPath1r1.split('/')[:-1])
        baseName=fileName.split(".")[0]
        ext=".".join(fileName.split(".")[1:-1])

        global reportDF
        reportDF["sample"]=args.pathogenR1.split("/")[-1].split(".")[0]
        resultsRow=wrapper(outDir,baseName,dirPath,fileName,args.minLen,args,args.pathogenR1,args.hostR1,args.splicedR1,"r1")
        if not args.quiet:
            print("report after wrapper r1")
            printReport()
        reportDF.to_csv(os.path.abspath(args.out)+".report",index=False)

        reportDF=pd.DataFrame([[]])
        reportDF["sample"]=args.pathogenR1.split("/")[-1].split(".")[0]
        resultsRow=wrapper(outDir,baseName,dirPath,fileName,args.minLen,args,args.pathogenR2,args.hostR2,args.splicedR2,"r2")
        if not args.quiet:
            print("report after wrapper r2")
            printReport()
        reportDF.to_csv(os.path.abspath(args.out)+".report",index=False)
        
        wrapperSpan(outDir,baseName,dirPath,fileName,args.minLen,args)

    else:
        print('not real path')
        return

def chimFinder(argv):

    parser=argparse.ArgumentParser(description='''Help Page''')

#==========================================
#==================STEP1===================
#==== Take two alignments and output ======
#=========== suggested chimeras ===========
#==========================================
    parser.add_argument('--pathogenR1',
                              required=True,
                              type=str,
                              help="Alignments of Mate 1 reads to pathogen genome")
    parser.add_argument('--pathogenR2',
                              required=True,
                              type=str,
                              help="Alignments of Mate 2 reads to pathogen genome")
    parser.add_argument('--hostR1',
                              required=True,
                              type=str,
                              help="Alignments of Mate 1 reads to host genome")
    parser.add_argument('--hostR2',
                              required=True,
                              type=str,
                              help="Alignments of Mate 2 reads to host genome")
    parser.add_argument('-o',
                              '--out',
                              required=False,
                              type=str,
                              default="./out",
                              help="output file")
    parser.add_argument('-t',
                              '--threads',
                              required=False,
                              default=1,
                              type=int,
                              help="the number of threads to use in the computation")
    parser.add_argument('--minCount',
                              required=False,
                              default=1,
                              type=float,
                              help="the minimum number of reads supporting an integration site.")
    parser.add_argument('--maxCountPenalty',
                              required=False,
                              default=0.85,
                              type=float,
                              help="the maximum penalty to give for minCount when hit.")
    parser.add_argument('--weightCount',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the number of reads score in the cumulative score equation")
    parser.add_argument('--minEntropy',
                              required=False,
                              default=0.8,
                              type=float,
                              help="minimum alignment entropy score.")
    parser.add_argument('--weightEntropy',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the entropy score in the cumulative score equation")
    parser.add_argument('--minQual',
                              required=False,
                              default=10,
                              type=float,
                              help="minimum mean mapping quality of sequenced read. Everything below this threshold will be reported as multialignment for which no assumption can be made from the annotation")
    parser.add_argument('--weightQual',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the mapping score in the cumulative score equation")
    parser.add_argument('--minLen',
                              required=False,
                              default=30,
                              type=float,
                              help="the minimum number of nucleotides in alignment to keep a read.")
    parser.add_argument('--maxLenUnmapped',
                              required=False,
                              default=30,
                              type=float,
                              help="minimum percent of the read to be aligned correctly in spanning read search")
    parser.add_argument('--weightLen',
                              required=False,
                              default=1.0,
                              type=float,
                              help="Weight of the alignment length score in the cumulative score equation")
    parser.add_argument('--steepSlopeAL',
                              required=False,
                              default=0.1,
                              type=float,
                              help="Weight of the alignment length score in the cumulative score equation")
    parser.add_argument('--maxAlLenPenalty',
                              required=False,
                              default=0.0,
                              type=float,
                              help="the maximum penalty to give for minLen when hit.")
    parser.add_argument('-s',
                              '--score',
                              required=False,
                              default=0.5,
                              type=float,
                              help="the minimum overall score to keep.")
    parser.add_argument('--overlap',
                              required=False,
                              default=-1,
                              type=int,
                              help="overlap threshold")
    parser.add_argument('--gap',
                              required=False,
                              default=-1,
                              type=int,
                              help="gap threshold")
    parser.add_argument('--close',
                              required=False,
                              default=30,
                              type=int,
                              help="distance between two integration sites to group together")
    parser.add_argument('-a',
                              '--annotation',
                              required=True,
                              type=str,
                              help="annotation for the host genome")
    parser.add_argument('-w',
                              '--writeReads',
                              action="store_true",
                              help="write reads to fasta files")
    parser.add_argument('-p',
                              '--plot',
                              action="store_true",
                              help="plot snapshots")
    parser.add_argument('-q',
                              '--quiet',
                              action="store_true",
                              help="do not print to std out. report will still be saved")
    parser.add_argument('--splicedR1',
                              required=False,
                              type=str,
                              help="spliced end-to-end alignment r1. Reads will be subtracted from the main alignments and will not be reported in the final report of integrations sites")
    parser.add_argument('--splicedR2',
                              required=False,
                              type=str,
                              help="spliced end-to-end alignment r2. Reads will be subtracted from the main alignments and will not be reported in the final report of integrations sites")
    parser.set_defaults(func=main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    chimFinder(sys.argv[1:])

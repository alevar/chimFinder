#!/usr/bin/env python

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

# the report function should report statistics about the each run
# The following fields can be filled in:
# 1. initial number of reads
# 2. read length mean prefiltering
# 3. read length std prefiltering
# 4. read length max prefiltering
# 5. read length min prefiltering
# 6. number of reads removed due to observed end-to-end splicing in human genome
# 7. number of reads filtered due to overlap/gap threshold
# 8. number of reads filtered due to insufficient alignment length
# 9. number of reads filtered due to complexity
# 10. number of reads filtered due to cumulative score
# 11. see if there are any other places where filtering might occur.
# 12. read length mean postfiltering
# 13. read length std postfiltering
# 14. read length max postfiltering
# 15. read length min postfiltering
reportDF=pd.DataFrame([[]])

def printReport():
    global reportDF
    for i in list(reportDF):
        print(i+str(": ")+str(reportDF[i].iloc[0]))

# Instead just write to a report file on separate lines
# those files can then be used to parse by the results.sh script

# Find chimeric reads (even if < 31nt alignment length to hiv or hum) from full alignments
def getFromSplits(dataG2,dataG1):
    cleanG1=dataG1[(dataG1["paired"]==1)&(dataG1["secondaryAlignment"]==0)&(dataG1["aligned2Mates"]==0)&(dataG1["unmappedCurr"]==0)&(dataG1["unmappedMate"]==8)]
    cleanG1_First=cleanG1[(cleanG1["firstRead"]==64)&(cleanG1["lastRead"]==0)]
    cleanG1_Last=cleanG1[(cleanG1["firstRead"]==0)&(cleanG1["lastRead"]==128)]

    cleanG2=dataG2[(dataG2["paired"]==1)&(dataG2["secondaryAlignment"]==0)&(dataG2["aligned2Mates"]==0)&(dataG2["unmappedCurr"]==0)&(dataG2["unmappedMate"]==8)]
    cleanG2_First=cleanG2[(cleanG2["firstRead"]==64)&(cleanG2["lastRead"]==0)]
    cleanG2_Last=cleanG2[(cleanG2["firstRead"]==0)&(cleanG2["lastRead"]==128)]
    #find mates within paired-end reads which aligned to both hiv and hum
    setG2_First=set(cleanG2_First["QNAME"])
    setG1_First=set(cleanG1_First["QNAME"])
    splitMates=setG2_First.intersection(setG1_First)

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

# mark reads that have G2 on the right side
def leftRight(data,minLen):
    data["G2"]=""
    data["G2_R1"]=""
    data["G2_R2"]=""
    data["sepR1"]=""
    data["sepR2"]=""

    data.loc[~(data["R1_G1_ID"]==0)&~(data["R1_G2_ID"]==0)&(data["R1_G2_TS"]-data["R1_G1_TS"]>minLen)&(data["R1_G2_TE"]-data["R1_G1_TE"]>minLen),"G2_R1"]="R1:r"
    data.loc[~(data["R1_G1_ID"]==0)&~(data["R1_G2_ID"]==0)&(data["R1_G1_TS"]-data["R1_G2_TS"]>minLen)&(data["R1_G1_TE"]-data["R1_G2_TE"]>minLen),"G2_R1"]="R1:l"

    data.loc[~(data["R2_G1_ID"]==0)&~(data["R2_G2_ID"]==0)&(data["R2_G2_TS"]-data["R2_G1_TS"]>minLen)&(data["R2_G2_TE"]-data["R2_G1_TE"]>minLen),"G2_R2"]="R2:r"
    data.loc[~(data["R2_G1_ID"]==0)&~(data["R2_G2_ID"]==0)&(data["R2_G1_TS"]-data["R2_G2_TS"]>minLen)&(data["R2_G1_TE"]-data["R2_G2_TE"]>minLen),"G2_R2"]="R2:l"

    data.loc[(~(data["R1_G1_ID"].astype(str)=="0")&~(data["R2_G2_ID"].astype(str)=="0")&(data["R2_G1_ID"].astype(str)=="0")), "sepR1"]="sepR1:2"
    data.loc[(~(data["R2_G1_ID"].astype(str)=="0")&~(data["R1_G2_ID"].astype(str)=="0")&(data["R1_G1_ID"].astype(str)=="0")), "sepR2"]="sepR2:1"

    data["G2"]=data["G2_R1"]+data["G2_R2"]+data["sepR1"]+data["sepR2"]
    data=data[data['G2'].str.len()>0]
    return data.drop(["G2_R1","G2_R2","sepR1","sepR2"],axis=1)

def leftRightUnpaired(data,minLen):
    data["G2"]=""

    data.loc[~(data["G1_ID"]==0)&~(data["G2_ID"]==0)&(data["G2_TS"]-data["G1_TS"]>minLen)&(data["G2_TE"]-data["G1_TE"]>minLen),"G2"]="r"
    data.loc[~(data["G1_ID"]==0)&~(data["G2_ID"]==0)&(data["G1_TS"]-data["G2_TS"]>minLen)&(data["G1_TE"]-data["G2_TE"]>minLen),"G2"]="l"

    data=data[data['G2'].str.len()>0]
    return data

def overlap(data):
    data.reset_index(inplace=True,drop=True)
    data["overlapR1"]=data[['R1_G1_TE','R1_G2_TE']].min(axis=1)-data[['R1_G1_TS','R1_G2_TS']].max(axis=1) #min end - max start
    data["gapR1"]=0
    data.loc[data['overlapR1']<0,['gapR1']]=data.loc[data['overlapR1']<0,['overlapR1']]["overlapR1"].abs()
    data.loc[data['overlapR1']<0,['overlapR1']]=0
    data["overlapR2"]=data[['R2_G1_TE','R2_G2_TE']].min(axis=1)-data[['R2_G1_TS','R2_G2_TS']].max(axis=1)
    data["gapR2"]=0
    data.loc[data['overlapR2']<0,['gapR2']]=data.loc[data['overlapR2']<0,['overlapR2']]["overlapR2"].abs()
    data.loc[data['overlapR2']<0,['overlapR2']]=0
    return data

def overlapUnpaired(data):
    data.reset_index(inplace=True,drop=True)
    data["overlap"]=data[['G1_TE','G2_TE']].min(axis=1)-data[['G1_TS','G2_TS']].max(axis=1) #min end - max start
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

# def extractStartEnd(data):
#     data["CIGAR"].replace("*",np.nan,inplace=True)
#     data.dropna(axis=0,inplace=True)

#     data["READ_LEN"]=data.SEQ.str.len()
#     data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$").replace(np.nan,0).astype(int)
#     data["END"]=data.READ_LEN-data.CIGAR_POST
#     data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[S]").replace(np.nan,0).astype(int)

#     data16=data[data["reversedCurr"]==16]
#     data0=data[data["reversedCurr"]==0]
#     data16["Template_start"]=data16.READ_LEN-data16.END
#     data16["Template_end"]=data16.READ_LEN-data16.CIGAR_PRE-1
#     data0["Template_start"]=data0.CIGAR_PRE
#     data0["Template_end"]=data0.END

#     data16["Reference_start"]=data16.READ_LEN-data16.END+data16.POS-data16.Template_start
#     data16["Reference_end"]=data16.READ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start
#     data0["Reference_start"]=data0.POS
#     data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE

#     data=pd.concat([data16,data0]).reset_index(drop=True)
#     data.drop(["CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
#     return data

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
def filterReads(dataG1,dataG2):
    #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
    dataG1=dataG1[(dataG1["secondaryAlignment"]==0)&(dataG1["PCRdup"]==0)&(dataG1["suppAl"]==0)&(dataG1["noPassFilter"]==0)]
    dataG2=dataG2[(dataG2["secondaryAlignment"]==0)&(dataG2["PCRdup"]==0)&(dataG2["suppAl"]==0)&(dataG2["noPassFilter"]==0)]
    return dataG1, dataG2

def createData(data,dataG1,dataG2):
    dataG1=dataG1[dataG1['QNAME'].isin(set(data['QNAME']))]
    dataG1R1=pd.DataFrame([])
    dataG1R1[["QNAME",
                "R1_G1_TS",
                "R1_G1_TE",
                "R1_G1_ID",
                "R1_G1_RS",
                "R1_G1_RE",
                "READ_LEN_R1",
                "R1_G1_SEQ",
                "QUAL_R1",
                "R1_G1_MAPQ",
                "R1_G1_reversedCurr"]]=dataG1[dataG1['firstRead']==64][["QNAME",
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
    dataG1R2=pd.DataFrame([])
    dataG1R2[["QNAME",
                "R2_G1_TS",
                "R2_G1_TE",
                "R2_G1_ID",
                "R2_G1_RS",
                "R2_G1_RE",
                "READ_LEN_R2",
                "R2_G1_SEQ",
                "QUAL_R2",
                "R2_G1_MAPQ",
                "R2_G1_reversedCurr"]]=dataG1[dataG1['lastRead']==128][["QNAME",
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
    dataG2R1=pd.DataFrame([])
    dataG2R1[["QNAME",
                "R1_G2_TS",
                "R1_G2_TE",
                "R1_G2_ID",
                "R1_G2_RS",
                "R1_G2_RE",
                "R1_G2_SEQ",
                "R1_G2_MAPQ",
                "R1_G2_reversedCurr"]]=dataG2[dataG2['firstRead']==64][["QNAME",
                                                                            "Template_start",
                                                                            "Template_end",
                                                                            "RNAME",
                                                                            "Reference_start",
                                                                            "Reference_end",
                                                                            "SEQ",
                                                                            "MAPQ",
                                                                            "reversedCurr"]]
    dataG2R2=pd.DataFrame([])
    dataG2R2[["QNAME",
                "R2_G2_TS",
                "R2_G2_TE",
                "R2_G2_ID",
                "R2_G2_RS",
                "R2_G2_RE",
                "R2_G2_SEQ",
                "R2_G2_MAPQ",
                "R2_G2_reversedCurr"]]=dataG2[dataG2['lastRead']==128][["QNAME",
                                                                        "Template_start",
                                                                        "Template_end",
                                                                        "RNAME",
                                                                        "Reference_start",
                                                                        "Reference_end",
                                                                        "SEQ",
                                                                        "MAPQ",
                                                                        "reversedCurr"]]
    data=data.merge(dataG1R1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataG1R2,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataG2R1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataG2R2,left_on='QNAME', right_on='QNAME', how='left')
    data[["R1_G1_ID",
        "R2_G1_ID",
        "R1_G2_ID",
        "R2_G2_ID",
        "R1_G1_MAPQ",
        "R2_G1_MAPQ",
        "R1_G2_MAPQ",
        "R2_G2_MAPQ",
        "R1_G1_reversedCurr",
        "R2_G1_reversedCurr",
        "R1_G2_reversedCurr",
        "R2_G2_reversedCurr"]]=data[["R1_G1_ID",
                                    "R2_G1_ID",
                                    "R1_G2_ID",
                                    "R2_G2_ID",
                                    "R1_G1_MAPQ",
                                    "R2_G1_MAPQ",
                                    "R1_G2_MAPQ",
                                    "R2_G2_MAPQ",
                                    "R1_G1_reversedCurr",
                                    "R2_G1_reversedCurr",
                                    "R1_G2_reversedCurr",
                                    "R2_G2_reversedCurr"]].fillna('',inplace=True)
    data.fillna(0,inplace=True)
    return data

def createDataUnpaired(data,dataG1,dataG2):
    dataG1=dataG1[dataG1['QNAME'].isin(set(data['QNAME']))].reset_index(drop=True)
    dfG1=pd.DataFrame([])
    dfG2=pd.DataFrame([])
    dfG1[["QNAME",
            "G1_TS",
            "G1_TE",
            "G1_ID",
            "G1_RS",
            "G1_RE",
            "READ_LEN",
            "G1_SEQ",
            "QUAL",
            "G1_MAPQ",
            "G1_reversedCurr"]]=dataG1[["QNAME",
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
    dfG2[["QNAME",
            "G2_TS",
            "G2_TE",
            "G2_ID",
            "G2_RS",
            "G2_RE",
            "G2_SEQ",
            "G2_MAPQ",
            "G2_reversedCurr"]]=dataG2[["QNAME",
                                            "Template_start",
                                            "Template_end",
                                            "RNAME",
                                            "Reference_start",
                                            "Reference_end",
                                            "SEQ",
                                            "MAPQ",
                                            "reversedCurr"]]
    data=data.merge(dfG1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dfG2,left_on='QNAME', right_on='QNAME', how='left')
    data[["G1_ID",
        "G2_ID",
        "G1_MAPQ",
        "G2_MAPQ",
        "G1_reversedCurr",
        "G2_reversedCurr"]]=data[["G1_ID",
                                    "G2_ID",
                                    "G1_MAPQ",
                                    "G2_MAPQ",
                                    "G1_reversedCurr",
                                    "G2_reversedCurr"]].fillna(value='')
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

    # expected=0.1#expectedEntropy(maxSubstringLen)
    # ent=bool(res>=expected)
    # res=res if ent else 0
    # temp="CACTTACATTGGGGAGTCAGGCTTCTCATCCACAGCCATGCCGTTCACACCCAGTCGCCGCCCCTCGCCTCTTGCTGTGCGCGCTTCAGCAAGCCGAGTCCTGCGTCGAGAGATCTCCTCTGGTTTTCCTTTCGCTTTCAGGTCCCTGCCG"
    # if s in temp:
    #     print(maxSubstringLen,expected,res)

    return res

# How to compute the mimimum entropy for a string of a given length and characters

#consider using technique described in https://academic.oup.com/bioinformatics/article/27/8/1061/227307/Topological-entropy-of-DNA-sequences
#another paper: Look at figure 2 https://www.nature.com/articles/srep19788#f2
#more https://www.xycoon.com/normalized_entropy.htm

def processAligns(seqG1,seqG2,qual,i1_g2,i2_g2,i1_g1,i2_g1,readLen,rG1,rG2):

    entropyScore_g2=0
    entropyScore_g1=0
    meanQual_g2=0
    meanQual_g1=0

    s_g2=""
    if rG2==True: # if reverse complemented take into account
        s_g2=seqG2[readLen-i2_g2:readLen-i1_g2]
    else:
        s_g2=seqG2[i1_g2:i2_g2]
    if not len(s_g2)==0:
        entropyScore_g2=topologicalNormalizedEntropy(s_g2)
        q_g2=qual[i1_g2:i2_g2]
        if len(q_g2)>0:
            meanQual_g2=sum([ord(x)-33 for x in q_g2])/len(q_g2)
        else:
            meanQual_g2=0

    s_g1=""
    if rG1==True: # if reverse complemented take into account
        s_g1=seqG1[readLen-i2_g1:readLen-i1_g1]
    else:
        s_g1=seqG1[i1_g1:i2_g1]
    if not len(s_g1)==0:
        entropyScore_g1=topologicalNormalizedEntropy(s_g1)
        q_g1=qual[i1_g1:i2_g1]
        if len(q_g1)>0:
            meanQual_g1=sum([ord(x)-33 for x in q_g1])/len(q_g1)
        else:
            meanQual_g1=0

    return pd.Series({"entropyScore_g2":entropyScore_g2,
                        "meanQual_g2":meanQual_g2,
                        "entropyScore_g1":entropyScore_g1,
                        "meanQual_g1":meanQual_g1})


# this function filters by overlap and flanking alignment length
# further it removes unnecessary data and combines into a single full dataframe with unified naming avoiding R1/R2 conventions
# output can be saved as .full.csv and then grouped by split position all at once

# this function allows identifying best reads
# However, since greater minLen values for G2 and G1 alignments will yield fewer reads but at higher confidence
# support reads could ignore the min len requirement as long as the split position is identical
# or perhaps the minLen requirement for the support reads should be lower
def filterOverlapCombine(data,args):
    dropList=["R1_G1_TS",
              "R1_G1_TE",
              "R1_G1_ID",
              "R1_G1_RS",
              "R1_G1_RE",
              "R2_G1_TS",
              "R2_G1_TE",
              "R2_G1_ID",
              "R2_G1_RS",
              "R2_G1_RE",
              "R1_G2_TS",
              "R1_G2_TE",
              "R1_G2_ID",
              "R1_G2_RS",
              "R1_G2_RE",
              "R2_G2_TS",
              "R2_G2_TE",
              "R2_G2_ID",
              "R2_G2_RS",
              "R2_G2_RE",
              "overlapR1",
              "overlapR2",
              "gapR1",
              "gapR2",
              "R1_G1_SEQ",
              "R2_G1_SEQ",
              "R1_G2_SEQ",
              "R2_G2_SEQ",
              "QUAL_R1",
              "QUAL_R2",
              "R1_G1_MAPQ",
              "R1_G2_MAPQ",
              "R2_G1_MAPQ",
              "R2_G2_MAPQ",
              "READ_LEN_R1",
              "READ_LEN_R2",
              "R1_G1_reversedCurr",
              "R2_G1_reversedCurr",
              "R1_G2_reversedCurr",
              "R2_G2_reversedCurr"]

    # R1 right
    data['entropyScore_g2']=0
    data['meanQual_g2']=0
    data['mapQual_g2']=0
    data['entropyScore_g1']=0
    data['meanQual_g1']=0
    data['mapQual_g1']=0
    dataR1Right=data[data["G2"].str.contains("R1:r")]
    dataR1Right=dataR1Right[~((dataR1Right["R1_G1_TS"]<dataR1Right["R1_G2_TS"])&(dataR1Right["R1_G1_TE"]>dataR1Right["R1_G2_TE"]))]
    dataR1Right=dataR1Right[~((dataR1Right["R1_G2_TS"]<dataR1Right["R1_G1_TS"])&(dataR1Right["R1_G2_TE"]>dataR1Right["R1_G1_TE"]))]
    dataR1Right["ins"]=dataR1Right["R1_G2_TS"]-dataR1Right["R1_G1_TE"]
    dataR1Right["split"]=dataR1Right['R1_G1_RE'].astype(int).astype(str)+":"+dataR1Right['R1_G2_RS'].astype(int).astype(str)
    dataR1Right['G1']=dataR1Right['R1_G1_RE'].astype(int)
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1_G1_ID+":"+dataR1Right.R1_G2_ID
    dataR1Right["orient"]="R1-g1:g2"
    dataR1Right["overlap"]=dataR1Right["overlapR1"]
    dataR1Right["gap"]=dataR1Right["gapR1"]
    dataR1Right["G1_TS"]=dataR1Right["R1_G1_TS"].astype(int)
    dataR1Right["G1_TE"]=dataR1Right["R1_G1_TE"].astype(int)
    dataR1Right["G1_RS"]=dataR1Right["R1_G1_RS"].astype(int)
    dataR1Right["G1_RE"]=dataR1Right["R1_G1_RE"].astype(int)
    dataR1Right["G2_TS"]=dataR1Right["R1_G2_TS"].astype(int)
    dataR1Right["G2_TE"]=dataR1Right["R1_G2_TE"].astype(int)
    dataR1Right["G2_RS"]=dataR1Right["R1_G2_RS"].astype(int)
    dataR1Right["G2_RE"]=dataR1Right["R1_G2_RE"].astype(int)
    dataR1Right["G1_ID"]=dataR1Right["R1_G1_ID"]
    dataR1Right["G2_ID"]=dataR1Right["R1_G2_ID"]
    dataR1Right["SEQ"]=dataR1Right["R1_G1_SEQ"]
#     dataR1Right["SEQ"]=dataR1Right[dataR1Right['R1_G1_reversedCurr']==0]["R1_G1_SEQ"]
#     dataR1Right["R_SEQ"]=dataR1Right[dataR1Right['R1_G1_reversedCurr']==16]["R1_G1_SEQ"]
    dataR1Right["G2_MAPQ"]=dataR1Right["R1_G2_MAPQ"].astype(int)
    dataR1Right["G1_MAPQ"]=dataR1Right["R1_G1_MAPQ"].astype(int)
    dataR1Right["G2_AL"]=dataR1Right["G2_TE"]-dataR1Right["G2_TS"]-dataR1Right["overlap"]
    dataR1Right["G1_AL"]=dataR1Right["G1_TE"]-dataR1Right["G1_TS"]-dataR1Right["overlap"]
    dataR1Right=dataR1Right[(dataR1Right['G2_AL']>args.minLen)&(dataR1Right['G1_AL']>args.minLen)]
    if len(dataR1Right>0):
        dataR1Right[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataR1Right.merge(dataR1Right.apply(lambda row: processAligns(row['R1_G1_SEQ'],
                                                                                                    row['R1_G2_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['G2_TS']+row['overlap']),
                                                                                                    int(row['G2_TE']),
                                                                                                    int(row['G1_TS']),
                                                                                                    int(row['G1_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1_G1_reversedCurr"]==16,
                                                                                                    row["R1_G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                                   'meanQual_g1_y']]
    dataR1Right["READ_LEN"]=dataR1Right["READ_LEN_R1"]
    dataR1Right["R"]="R1"
    dataR1Right.drop(dropList,axis=1,inplace=True)

    # R1 left
    dataR1Left=data[data["G2"].str.contains("R1:l")]
    dataR1Left=dataR1Left[~((dataR1Left["R1_G1_TS"]<dataR1Left["R1_G2_TS"])&(dataR1Left["R1_G1_TE"]>dataR1Left["R1_G2_TE"]))]
    dataR1Left=dataR1Left[~((dataR1Left["R1_G2_TS"]<dataR1Left["R1_G1_TS"])&(dataR1Left["R1_G2_TE"]>dataR1Left["R1_G1_TE"]))]
    dataR1Left["ins"]=dataR1Left["R1_G1_TS"]-dataR1Left["R1_G2_TE"]
    dataR1Left["split"]=dataR1Left['R1_G2_RE'].astype(int).astype(str)+":"+dataR1Left['R1_G1_RS'].astype(int).astype(str)
    dataR1Left['G1']=dataR1Left['R1_G1_RS'].astype(int)
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1_G2_ID+":"+dataR1Left.R1_G1_ID
    dataR1Left["orient"]="R1-g2:g1"
    dataR1Left["overlap"]=dataR1Left["overlapR1"]
    dataR1Left["gap"]=dataR1Left["gapR1"]
    dataR1Left["G1_TS"]=dataR1Left["R1_G1_TS"].astype(int)
    dataR1Left["G1_TE"]=dataR1Left["R1_G1_TE"].astype(int)
    dataR1Left["G1_RS"]=dataR1Left["R1_G1_RS"].astype(int)
    dataR1Left["G1_RE"]=dataR1Left["R1_G1_RE"].astype(int)
    dataR1Left["G2_TS"]=dataR1Left["R1_G2_TS"].astype(int)
    dataR1Left["G2_TE"]=dataR1Left["R1_G2_TE"].astype(int)
    dataR1Left["G2_RS"]=dataR1Left["R1_G2_RS"].astype(int)
    dataR1Left["G2_RE"]=dataR1Left["R1_G2_RE"].astype(int)
    dataR1Left["G1_ID"]=dataR1Left["R1_G1_ID"]
    dataR1Left["G2_ID"]=dataR1Left["R1_G2_ID"]
    dataR1Left["SEQ"]=dataR1Left["R1_G1_SEQ"]
#     dataR1Left["SEQ"]=dataR1Left[dataR1Left['R1_G1_reversedCurr']==0]["R1_G1_SEQ"]
#     dataR1Left["R_SEQ"]=dataR1Left[dataR1Left['R1_G1_reversedCurr']==16]["R1_G1_SEQ"]
    dataR1Left["G2_MAPQ"]=dataR1Left["R1_G2_MAPQ"].astype(int)
    dataR1Left["G1_MAPQ"]=dataR1Left["R1_G1_MAPQ"].astype(int)
    dataR1Left["G2_AL"]=dataR1Left["G2_TE"]-dataR1Left["G2_TS"]-dataR1Left["overlap"]
    dataR1Left["G1_AL"]=dataR1Left["G1_TE"]-dataR1Left["G1_TS"]-dataR1Left["overlap"]
    dataR1Left=dataR1Left[(dataR1Left["G2_AL"]>args.minLen)&(dataR1Left["G1_AL"]>args.minLen)]
    if len(dataR1Left)>0:
        dataR1Left[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataR1Left.merge(dataR1Left.apply(lambda row: processAligns(row['R1_G1_SEQ'],
                                                                                                    row['R1_G2_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['G2_TS']),
                                                                                                    int(row['G2_TE']-row['overlap']),
                                                                                                    int(row['G1_TS']+row['overlap']),
                                                                                                    int(row['G1_TE']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1_G1_reversedCurr"]==16,
                                                                                                    row["R1_G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                   'meanQual_g1_y']]

    dataR1Left["READ_LEN"]=dataR1Left["READ_LEN_R1"]
    dataR1Left["R"]="R1"
    dataR1Left.drop(dropList,axis=1,inplace=True)

    # R2 right
    dataR2Right=data[data["G2"].str.contains("R2:r")]
    dataR2Right=dataR2Right[~((dataR2Right["R2_G1_TS"]<dataR2Right["R2_G2_TS"])&(dataR2Right["R2_G1_TE"]>dataR2Right["R2_G2_TE"]))]
    dataR2Right=dataR2Right[~((dataR2Right["R2_G2_TS"]<dataR2Right["R2_G1_TS"])&(dataR2Right["R2_G2_TE"]>dataR2Right["R2_G1_TE"]))]
    dataR2Right["ins"]=dataR2Right["R2_G2_TS"]-dataR2Right["R2_G1_TE"]
    dataR2Right["split"]=dataR2Right['R2_G1_RE'].astype(int).astype(str)+":"+dataR2Right['R2_G2_RS'].astype(int).astype(str)
    dataR2Right['G1']=dataR2Right['R2_G1_RE'].astype(int)
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2_G1_ID+":"+dataR2Right.R2_G2_ID
    dataR2Right["orient"]="R2-g1:g2"
    dataR2Right["overlap"]=dataR2Right["overlapR2"]
    dataR2Right["gap"]=dataR2Right["gapR2"]
    dataR2Right["G1_TS"]=dataR2Right["R2_G1_TS"].astype(int)
    dataR2Right["G1_TE"]=dataR2Right["R2_G1_TE"].astype(int)
    dataR2Right["G1_RS"]=dataR2Right["R2_G1_RS"].astype(int)
    dataR2Right["G1_RE"]=dataR2Right["R2_G1_RE"].astype(int)
    dataR2Right["G2_TS"]=dataR2Right["R2_G2_TS"].astype(int)
    dataR2Right["G2_TE"]=dataR2Right["R2_G2_TE"].astype(int)
    dataR2Right["G2_RS"]=dataR2Right["R2_G2_RS"].astype(int)
    dataR2Right["G2_RE"]=dataR2Right["R2_G2_RE"].astype(int)
    dataR2Right["G1_ID"]=dataR2Right["R2_G1_ID"]
    dataR2Right["G2_ID"]=dataR2Right["R2_G2_ID"]
    dataR2Right["SEQ"]=dataR2Right["R2_G1_SEQ"]
#     dataR2Right["SEQ"]=dataR2Right[dataR2Right['R2_G1_reversedCurr']==0]["R2_G1_SEQ"]
#     dataR2Right["R_SEQ"]=dataR2Right[dataR2Right['R2_G1_reversedCurr']==16]["R2_G1_SEQ"]
    dataR2Right["G2_MAPQ"]=dataR2Right["R2_G2_MAPQ"].astype(int)
    dataR2Right["G1_MAPQ"]=dataR2Right["R2_G1_MAPQ"].astype(int)
    dataR2Right["G2_AL"]=dataR2Right["G2_TE"]-dataR2Right["G2_TS"]-dataR2Right["overlap"]
    dataR2Right["G1_AL"]=dataR2Right["G1_TE"]-dataR2Right["G1_TS"]-dataR2Right["overlap"]
    dataR2Right=dataR2Right[(dataR2Right["G2_AL"]>args.minLen)&(dataR2Right["G1_AL"]>args.minLen)]
    if len(dataR2Right)>0:
        dataR2Right[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataR2Right.merge(dataR2Right.apply(lambda row: processAligns(row['R2_G1_SEQ'],
                                                                                                    row['R2_G2_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['G2_TS']+row['overlap']),
                                                                                                    int(row['G2_TE']),
                                                                                                    int(row['G1_TS']),
                                                                                                    int(row['G1_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2_G1_reversedCurr"]==16,
                                                                                                    row["R2_G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                                   'meanQual_g1_y']]
    dataR2Right["READ_LEN"]=dataR2Right["READ_LEN_R2"]
    dataR2Right["R"]="R2"
    dataR2Right.drop(dropList,axis=1,inplace=True)

    # R2 left
    dataR2Left=data[data["G2"].str.contains("R2:l")]
    dataR2Left=dataR2Left[~((dataR2Left["R2_G1_TS"]<dataR2Left["R2_G2_TS"])&(dataR2Left["R2_G1_TE"]>dataR2Left["R2_G2_TE"]))]
    dataR2Left=dataR2Left[~((dataR2Left["R2_G2_TS"]<dataR2Left["R2_G1_TS"])&(dataR2Left["R2_G2_TE"]>dataR2Left["R2_G1_TE"]))]
    dataR2Left["ins"]=dataR2Left["R2_G1_TS"]-dataR2Left["R2_G2_TE"]
    dataR2Left["split"]=dataR2Left['R2_G2_RE'].astype(int).astype(str)+":"+dataR2Left['R2_G1_RS'].astype(int).astype(str)
    dataR2Left['G1']=dataR2Left['R2_G1_RS'].astype(int)
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2_G2_ID+":"+dataR2Left.R2_G1_ID
    dataR2Left["orient"]="R2-g2:g1"
    dataR2Left["overlap"]=dataR2Left["overlapR2"]
    dataR2Left["gap"]=dataR2Left["gapR2"]
    dataR2Left["G1_TS"]=dataR2Left["R2_G1_TS"].astype(int)
    dataR2Left["G1_TE"]=dataR2Left["R2_G1_TE"].astype(int)
    dataR2Left["G1_RS"]=dataR2Left["R2_G1_RS"].astype(int)
    dataR2Left["G1_RE"]=dataR2Left["R2_G1_RE"].astype(int)
    dataR2Left["G2_TS"]=dataR2Left["R2_G2_TS"].astype(int)
    dataR2Left["G2_TE"]=dataR2Left["R2_G2_TE"].astype(int)
    dataR2Left["G2_RS"]=dataR2Left["R2_G2_RS"].astype(int)
    dataR2Left["G2_RE"]=dataR2Left["R2_G2_RE"].astype(int)
    dataR2Left["G1_ID"]=dataR2Left["R2_G1_ID"]
    dataR2Left["G2_ID"]=dataR2Left["R2_G2_ID"]
    dataR2Left["SEQ"]=dataR2Left["R2_G1_SEQ"]
#     dataR2Left["SEQ"]=dataR2Left[dataR2Left['R2_G1_reversedCurr']==0]["R2_G1_SEQ"]
#     dataR2Left["R_SEQ"]=dataR2Left[dataR2Left['R2_G1_reversedCurr']==16]["R2_G1_SEQ"]
    dataR2Left["G2_MAPQ"]=dataR2Left["R2_G2_MAPQ"].astype(int)
    dataR2Left["G1_MAPQ"]=dataR2Left["R2_G1_MAPQ"].astype(int)
    dataR2Left["G2_AL"]=dataR2Left["G2_TE"]-dataR2Left["G2_TS"]-dataR2Left["overlap"]
    dataR2Left["G1_AL"]=dataR2Left["G1_TE"]-dataR2Left["G1_TS"]-dataR2Left["overlap"]
    dataR2Left=dataR2Left[(dataR2Left["G2_AL"]>args.minLen)&(dataR2Left["G1_AL"]>args.minLen)]
    if len(dataR2Left)>0:
        dataR2Left[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataR2Left.merge(dataR2Left.apply(lambda row: processAligns(row['R2_G1_SEQ'],
                                                                                                    row['R2_G2_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['G2_TS']),
                                                                                                    int(row['G2_TE']-row['overlap']),
                                                                                                    int(row['G1_TS']+row['overlap']),
                                                                                                    int(row['G1_TE']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2_G1_reversedCurr"]==16,
                                                                                                    row["R2_G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                   'meanQual_g1_y']]

    
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
    
    # df["G2_AL"]=(df["READ_LEN"]/2)-((df["G2_AL"]-df["READ_LEN"]/2).abs())
    # df["G1_AL"]=(df["READ_LEN"]/2)-((df["G1_AL"]-df["READ_LEN"]/2).abs())

    k=float(args.minLen)
    ssAl=args.steepSlopeAL
    # df["G2_AL_score"]=(((((df['G2_AL']-k)/k) \
    #                     /(((ssAl)+((df['G2_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    df["G2_AL_score"]=((df['G2_AL']-k)/((ssAl+(df['G2_AL']-k)**2.0)**0.5)+1.0)/2.0
    df["G1_AL_score"]=((df['G1_AL']-k)/((ssAl+(df['G1_AL']-k)**2.0)**0.5)+1.0)/2.0

    # df["G1_AL_score"]=(((((df['G1_AL']-k)/k) \
    #                     /(((ssAl)+((df['G1_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of G2 alignment length score

    df['jointEntropy']=((df['entropyScore_g2'] \
                    +df['entropyScore_g1']) \
                    /(2))
    df['jointAlLen']=((df['G1_AL_score'] \
                    +df['G2_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.G1,df.G2,df.overlap,df.gap))))

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
    data['entropyScore_g2']=0
    data['meanQual_g2']=0
    data['mapQual_g2']=0
    data['entropyScore_g1']=0
    data['meanQual_g1']=0
    data['mapQual_g1']=0
    dataRight=data[data["G2"].str.contains("r")]
    dataRight=dataRight[~((dataRight["G1_TS"]<dataRight["G2_TS"])&(dataRight["G1_TE"]>dataRight["G2_TE"]))]
    dataRight=dataRight[~((dataRight["G2_TS"]<dataRight["G1_TS"])&(dataRight["G2_TE"]>dataRight["G1_TE"]))]
    dataRight["ins"]=dataRight["G2_TS"]-dataRight["G1_TE"]
    dataRight["split"]=dataRight['G1_RE'].astype(int).astype(str)+":"+dataRight['G2_RS'].astype(int).astype(str)
    dataRight['G1']=dataRight['G1_RE'].astype(int)
    dataRight["comb"]=dataRight.split+"@"+dataRight.G1_ID+":"+dataRight.G2_ID
    dataRight["orient"]="`"
    dataRight["overlap"]=dataRight["overlap"]
    dataRight["gap"]=dataRight["gap"]
    dataRight["G1_TS"]=dataRight["G1_TS"].astype(int)
    dataRight["G1_TE"]=dataRight["G1_TE"].astype(int)
    dataRight["G1_RS"]=dataRight["G1_RS"].astype(int)
    dataRight["G1_RE"]=dataRight["G1_RE"].astype(int)
    dataRight["G2_TS"]=dataRight["G2_TS"].astype(int)
    dataRight["G2_TE"]=dataRight["G2_TE"].astype(int)
    dataRight["G2_RS"]=dataRight["G2_RS"].astype(int)
    dataRight["G2_RE"]=dataRight["G2_RE"].astype(int)
    dataRight["G2_MAPQ"]=dataRight["G2_MAPQ"].astype(int)
    dataRight["G1_MAPQ"]=dataRight["G1_MAPQ"].astype(int)
    dataRight["SEQ"]=dataRight["G1_SEQ"]
#     dataRight["SEQ"]=dataRight[dataRight['G1_reversedCurr']==0]["G1_SEQ"]
#     dataRight["R_SEQ"]=dataRight[dataRight['G1_reversedCurr']==16]["G1_SEQ"]
    dataRight["G2_AL"]=dataRight["G2_TE"]-dataRight["G2_TS"]-dataRight["overlap"]
    dataRight["G1_AL"]=dataRight["G1_TE"]-dataRight["G1_TS"]-dataRight["overlap"]
    dataRight=dataRight[(dataRight['G2_AL']>args.minLen)&(dataRight['G1_AL']>args.minLen)]
    if len(dataRight>0):
        dataRight[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataRight.merge(dataRight.apply(lambda row: processAligns(row['G1_SEQ'],
                                                                                                    row['G2_SEQ'],
                                                                                                    row['QUAL'],
                                                                                                    int(row['G2_TS']+row['overlap']),
                                                                                                    int(row['G2_TE']),
                                                                                                    int(row['G1_TS']),
                                                                                                    int(row['G1_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN']),
                                                                                                    row["G1_reversedCurr"]==16,
                                                                                                    row["G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                                   'meanQual_g1_y']]
    dataRight["READ_LEN"]=dataRight["READ_LEN"]
    dataRight["R"]=0

    # left
    dataLeft=data[data["G2"].str.contains("l")]
    dataLeft=dataLeft[~((dataLeft["G1_TS"]<dataLeft["G2_TS"])&(dataLeft["G1_TE"]>dataLeft["G2_TE"]))]
    dataLeft=dataLeft[~((dataLeft["G2_TS"]<dataLeft["G1_TS"])&(dataLeft["G2_TE"]>dataLeft["G1_TE"]))]
    dataLeft["ins"]=dataLeft["G1_TS"]-dataLeft["G2_TE"]
    dataLeft["split"]=dataLeft['G2_RE'].astype(int).astype(str)+":"+dataLeft['G1_RS'].astype(int).astype(str)
    dataLeft['G1']=dataLeft['G1_RS'].astype(int)
    dataLeft["comb"]=dataLeft.split+"@"+dataLeft.G2_ID+":"+dataLeft.G1_ID
    dataLeft["orient"]="g2:g1"
    dataLeft["overlap"]=dataLeft["overlap"]
    dataLeft["gap"]=dataLeft["gap"]
    dataLeft["G1_TS"]=dataLeft["G1_TS"].astype(int)
    dataLeft["G1_TE"]=dataLeft["G1_TE"].astype(int)
    dataLeft["G1_RS"]=dataLeft["G1_RS"].astype(int)
    dataLeft["G1_RE"]=dataLeft["G1_RE"].astype(int)
    dataLeft["G2_TS"]=dataLeft["G2_TS"].astype(int)
    dataLeft["G2_TE"]=dataLeft["G2_TE"].astype(int)
    dataLeft["G2_RS"]=dataLeft["G2_RS"].astype(int)
    dataLeft["G2_RE"]=dataLeft["G2_RE"].astype(int)
    dataLeft["G2_MAPQ"]=dataLeft["G2_MAPQ"].astype(int)
    dataLeft["G1_MAPQ"]=dataLeft["G1_MAPQ"].astype(int)
    dataLeft["SEQ"]=dataLeft["G1_SEQ"]
#     dataLeft["SEQ"]=dataLeft[dataLeft['G1_reversedCurr']==0]["G1_SEQ"]
#     dataLeft["R_SEQ"]=dataLeft[dataLeft['G1_reversedCurr']==16]["G1_SEQ"]
    dataLeft["G2_AL"]=dataLeft["G2_TE"]-dataLeft["G2_TS"]-dataLeft["overlap"]
    dataLeft["G1_AL"]=dataLeft["G1_TE"]-dataLeft["G1_TS"]-dataLeft["overlap"]
    dataLeft=dataLeft[(dataLeft["G2_AL"]>args.minLen)&(dataLeft["G1_AL"]>args.minLen)]
    if len(dataLeft)>0:
        dataLeft[['entropyScore_g2',
                     'meanQual_g2',
                     'entropyScore_g1',
                     'meanQual_g1']]=dataLeft.merge(dataLeft.apply(lambda row: processAligns(row['G1_SEQ'],
                                                                                            row['G2_SEQ'],
                                                                                            row['QUAL'],
                                                                                            int(row['G2_TS']),
                                                                                            int(row['G2_TE']-row['overlap']),
                                                                                            int(row['G1_TS']+row['overlap']),
                                                                                            int(row['G1_TE']),
                                                                                            int(row['READ_LEN']),
                                                                                            row["G1_reversedCurr"]==16,
                                                                                            row["G2_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_g2_y',
                                                                                                                                                                   'meanQual_g2_y',
                                                                                                                                                                   'entropyScore_g1_y',
                                                                                                                                                                   'meanQual_g1_y']]

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
    # df["G2_AL_score"]=(((((df['G2_AL']-k)/k) \
    #                     /(((ssAl)+((df['G2_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    # df["G1_AL_score"]=(((((df['G1_AL']-k)/k) \
    #                     /(((ssAl)+((df['G1_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of G2 alignment length score

    df["G2_AL_score"]=((df['G2_AL']-k)/((ssAl+(df['G2_AL']-k)**2.0)**0.5)+1.0)/2.0
    df["G1_AL_score"]=((df['G1_AL']-k)/((ssAl+(df['G1_AL']-k)**2.0)**0.5)+1.0)/2.0

    df['jointEntropy']=((df['entropyScore_g2'] \
                    +df['entropyScore_g1']) \
                    /(2))
    df['jointAlLen']=((df['G1_AL_score'] \
                    +df['G2_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.G1,df.G2,df.overlap,df.gap))))

    # global reportDF
    # dataLen=df["SEQ"].str.len()
    # reportDF["number of reads after intermediate score cutoff"]=len(df)
    # reportDF["read len mean after intermediate score cutoff"]=dataLen.mean()
    # reportDF["read len std after intermediate score cutoff"]=dataLen.std()
    # reportDF["read len min after intermediate score cutoff"]=dataLen.min()
    # reportDF["read len max after intermediate score cutoff"]=dataLen.max()

    return df

def addSpan(data,dataPos):
    def testR1(row,dataG2):
        dataG2R1=dataG2[(dataG2['R1_G2_RS']-row['G2_RS']>-500)&(dataG2['R1_G2_RS']-row['G2_RS']<0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataG2R1=dataG2R1[(dataG2R1['R2_G1_RE']-row['G1_RE']<500)&(dataG2R1['R2_G1_RE']-row['G1_RE']>0)]
        if len(dataG2R1)>0:
            return [set(dataG2R1["QNAME"]),len(dataG2R1)] # return set of reads and count of reads
        return [{''},0]
        
    def testR2(row,dataG2):
        dataG2R2=dataG2[(dataG2['R2_G2_RE']-row['G2_RE']<500)&(dataG2['R2_G2_RE']-row['G2_RE']>0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataG2R2=dataG2R2[(dataG2R2['R1_G1_RS']-row['G1_RS']>-500)&(dataG2R2['R1_G1_RS']-row['G1_RS']<0)]
        if len(dataG2R2)>0:
            return [set(dataG2R2["QNAME"]),len(dataG2R2)] # return set of reads and count of reads
        return [{''},0]

    dataG2R1=data[data['G2'].str.contains('sepR2:1')]
    dataG2R1=dataG2R1[~(dataG2R1['R1_G2_ID']==0)] # check that the hiv is on the same side
    dataPosR1=dataPos[dataPos['orient'].str.contains('g2:g1')]
    if len(dataPosR1)>0:
        dataPosR1[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR1.apply(lambda row: testR1(row,dataG2R1),axis=1)])

    dataG2R2=data[data['G2'].str.contains('sepR1:2')]
    dataG2R2=dataG2R2[~(dataG2R2['R2_G2_ID']==0)] # check that the hiv is on the same side
    dataPosR2=dataPos[dataPos['orient'].str.contains('g1:g2')]
    if len(dataPosR2)>0:
        dataPosR2[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR2.apply(lambda row: testR2(row,dataG2R2),axis=1)])

    frames=[dataPosR1,dataPosR2]
    df=pd.concat(frames).reset_index(drop=True)
    return df

# this function should do the following:
# 1. group data by split position
# 2. add a column of high confidence support reads from the first parameter DF
# 3. add a column of low confidence support reads from the second parameter DF
def findSupport(data,minLen,unpaired):
    aggregations={
        'QNAME':{
            'count':'count',
            'reads': lambda x: set(x)
        },
        'SEQ':{
            'seq': lambda x: (max(dict(list(x)), key=int),dict(list(x))[max(dict(list(x)), key=int)])
        },
        'G1_RS':{
            'G1_RS':'min'
        },
        'G1_RE':{
            'G1_RE':'max'
        },
        'G2_RS':{
            'G2_RS':'min'
        },
        'G2_RE':{
            'G2_RE':'max'
        },
        'G2_AL':{
            'G2_AL':'sum'
        },
        'G1_AL':{
            'G1_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'entropyScore_g2':{
            'entropyScore_g2':'sum'
        },
        'entropyScore_g1':{
            'entropyScore_g1':'sum'
        },
        'meanQual_g2':{
            'meanQual_g2':'sum'
        },
        'meanQual_g1':{
            'meanQual_g1':'sum'
        },
        'G2_MAPQ':{
            'G2_MAPQ':'sum'
        },
        'G1_MAPQ':{
            'G1_MAPQ':'sum'
        }
    }

    dataPos=pd.DataFrame(data.groupby(by=["comb","split","G1_ID","R","orient"])[["QNAME",
                                                                                    "G1_RS",
                                                                                    "G1_RE",
                                                                                    "G2_RS",
                                                                                    "G2_RE",
                                                                                    "G1_AL",
                                                                                    "G2_AL",
                                                                                    "READ_LEN",
                                                                                    "entropyScore_g2",
                                                                                    "meanQual_g2",
                                                                                    "entropyScore_g1",
                                                                                    "meanQual_g1",
                                                                                    "G2_MAPQ",
                                                                                    "G1_MAPQ",
                                                                                    "SEQ"]].agg(aggregations)).reset_index()
    dataPos.rename(columns={'G1_ID':'chr'}, inplace=True)            
    return dataPos

def score(dataPos,args,minLen):
    dataPos["READ_LEN"]=dataPos["READ_LEN"].astype(float)/dataPos["count"].astype(float)
    dataPos["G2_MAPQ"]=dataPos["G2_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["G1_MAPQ"]=dataPos["G1_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["G2_AL"]=dataPos["G2_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["G1_AL"]=dataPos["G1_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_g2"]=dataPos["entropyScore_g2"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_g1"]=dataPos["entropyScore_g1"].astype(float)/dataPos["count"].astype(float)
    dataPos["count"]=dataPos["count"].astype(float)

    dataPos["G2_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["G2_AL"]-dataPos["READ_LEN"]/2).abs())
    dataPos["G1_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["G1_AL"]-dataPos["READ_LEN"]/2).abs())

    k=float(args.minLen)
    ssAl=args.steepSlopeAL
    # dataPos["G2_AL_score"]=(((((dataPos['G2_AL']-k)) \
    #                     /(((ssAl)+((dataPos['G2_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    # dataPos["G1_AL_score"]=(((((dataPos['G1_AL']-k)) \
    #                     /(((ssAl)+((dataPos['G1_AL']-k))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxAlLenPenalty))) \
    #                             +args.maxAlLenPenalty # algebraic sigmoid function of G2 alignment length score

    dataPos["G2_AL_score"]=((dataPos['G2_AL']-k)/((ssAl+(dataPos['G2_AL']-k)**2.0)**0.5)+1.0)/2.0
    dataPos["G1_AL_score"]=((dataPos['G1_AL']-k)/((ssAl+(dataPos['G1_AL']-k)**2.0)**0.5)+1.0)/2.0

    m=float(args.minCount)/2.0
    dataPos["count_score"]=((dataPos['count']-m)/((1.0+(dataPos['count']-m)**2.0)**0.5)+1.0)/2.0 \
                                *(1.0-args.maxCountPenalty) \
                                +args.maxCountPenalty # algebraic sigmoid function of read count score

    # dataPos["count_score"]=(((((dataPos['count']-m)) \
    #                     /(((1.0)+((dataPos['count']-m))**2.0)**0.5)) \
    #                             /2.0 \
    #                             +0.5) \
    #                             /(1.0/(1.0-args.maxCountPenalty))) \
    #                             +args.maxCountPenalty # algebraic sigmoid function of read count score
    
    # dataPos["HIV_MAPQ_score"]=1-(10**(-(dataPos["G2_MAPQ"]/10)))
    # dataPos["HUM_MAPQ_score"]=1-(10**(-(dataPos["G1_MAPQ"]/10)))

    dataPos['jointEntropy']=((dataPos['entropyScore_g2'] \
                    +dataPos['entropyScore_g1']) \
                    /(2))
    dataPos['jointAlLen']=((dataPos['G1_AL_score'] \
                    +dataPos['G2_AL_score']) \
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
    data.sort_values(by=["G1_RS","hum_nearest_SS"],inplace=True)
    data["diff1"]=abs(data['G1_RS']-data['G1_RS'].shift(-1))
    data["diff2"]=abs(data['G1_RS']-data['G1_RS'].shift(1))
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

def getSet(row):
    tmpLocalQNAME=list(row["reads"])
    return tmpLocalQNAME

def writeReadNames(outDir,row,fileNameR1,fileNameR2,baseName,dirPath):
    outD=os.path.abspath(outDir)
    tempF=outD+'/tmp/fq'
    stringReads=row['reads'].replace(';','|')

    outPath1_All="'"+outD+"/Positions/"+baseName+"/"+str(row['hum_nearest_SS'].strip('\n')+"@"+row["R"]+"@"+str(row['comb'].split("@")[0]))+"_R1.all.fa'"
    outPath2_All="'"+outD+"/Positions/"+baseName+"/"+str(row['hum_nearest_SS'].strip('\n')+"@"+row["R"]+"@"+str(row['comb'].split("@")[0]))+"_R2.all.fa'"
    cmdR1="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR1+" | seqtk seq -a - > "+outPath1_All
    cmdR2="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR2+" | seqtk seq -a - > "+outPath2_All
    os.system(cmdR1)
    os.system(cmdR2)

def writeReadNamesUnpaired(outDir,row,fileName,baseName,dirPath):
    outD=os.path.abspath(outDir)
    tempF=outD+'/tmp/fq'
    stringReads=row['reads'].replace(';','|')

    outPath="'"+outD+"/Positions/"+baseName+"/"+str(row['hum_nearest_SS'].strip('\n')+"@"+str(row['comb'].split("@")[0]))+".fa'"
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
    aggregations={
        'reads':{
            'groupsCount':'count',
            'reads': lambda x: ';'.join(set(x))
        },
        'seq':{
            'seq': lambda x: dict(list(x))[max(dict(list(x)), key=int)]
        },
        'spanR1-R2':{
            'spanReads':lambda x: ';'.join(list(set([el for el in x if not el=='-'])))
        },
        'spanCount':{
            'spanCount':'sum'
        },
        'count':{
            'count':'sum'
        },
        'comb':{
            'comb': lambda x: ';'.join(set(x))
        },
        'entropyScore_g1':{
            'entropyScore_g1':'sum'
        },
        'entropyScore_g2':{
            'entropyScore_g2':'sum'
        },
        'G2_AL':{
            'G2_AL':'sum'
        },
        'G1_AL':{
            'G1_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'G2_MAPQ':{
            'G2_MAPQ':'sum'
        },
        'G1_MAPQ':{
            'G1_MAPQ':'sum'
        }
    }
    
    dfg=pd.DataFrame(data.groupby(by=["hum_nearest_SS",
                                      "chr",
                                      "R",
                                      "uid"])[["comb",
                                                "reads",
                                                "count",
                                                "spanCount",
                                                "spanR1-R2",
                                                "entropyScore_g1",
                                                "entropyScore_g2",
                                                "G2_AL",
                                                "G1_AL",
                                                "READ_LEN",
                                                "G2_MAPQ",
                                                "G1_MAPQ",
                                                "seq"]].agg(aggregations)).reset_index()

    return dfg.reset_index(drop=True)

def groupBySpliceSitesUnpaired(data):
    colNames=["orient"]
        
    data.drop(colNames,axis=1,inplace=True)
    aggregations={
        'reads':{
            'groupsCount':'count',
            'reads': lambda x: ';'.join(set(x))
        },
        'seq':{
            'seq': lambda x: dict(list(x))[max(dict(list(x)), key=int)]
        },
        'count':{
            'count':'sum'
        },
        'comb':{
            'comb': lambda x: ';'.join(set(x))
        },
        'entropyScore_g1':{
            'entropyScore_g1':'sum'
        },
        'entropyScore_g2':{
            'entropyScore_g2':'sum'
        },
        'G2_AL':{
            'G2_AL':'sum'
        },
        'G1_AL':{
            'G1_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'G2_MAPQ':{
            'G2_MAPQ':'sum'
        },
        'G1_MAPQ':{
            'G1_MAPQ':'sum'
        }
    }
    
    dfg=pd.DataFrame(data.groupby(by=["hum_nearest_SS",
                                      "chr",
                                      "R",
                                      "uid"])[["comb",
                                                "reads",
                                                "count",
                                                "entropyScore_g1",
                                                "entropyScore_g2",
                                                "G2_AL",
                                                "G1_AL",
                                                "READ_LEN",
                                                "G2_MAPQ",
                                                "G1_MAPQ",
                                                "seq"]].agg(aggregations)).reset_index()

    return dfg.reset_index(drop=True)

def annotate(dataBed,annPath,data):
    data.reset_index(drop=True,inplace=True)
    sites=BedTool.from_dataframe(dataBed)
    annotation=BedTool(annPath)
    nearby=annotation.intersect(sites, wo=True)
    df=pd.read_table(nearby.fn,names=["chr",
                                      "remove1",
                                      "ord",
                                      "startRegionPos",
                                      "endRegionPos",
                                      "something1",
                                      "something2",
                                      "something3",
                                      "useful_info",
                                      "remove2",
                                      "startQuery",
                                      "endQuery",
                                      "distance"])
    if len(df)==0:
        return None

    order={"exon":0,
           "CDS":1,
           "mRNA":2,
           "gene":3,
           "ncRNA":4,
           "transcript":5,
           "primary_transcript":6,
           "rRNA":7,
           "tRNA":8,
           "locus":9}
    reverseOrder={0:"exon",
                  1:"CDS",
                  2:"mRNA",
                  3:"gene",
                  4:"ncRNA",
                  5:"transcript",
                  6:"primary_transcript",
                  7:"rRNA",
                  8:"tRNA",
                  9:"locus"}

    df=df.replace({'ord':order})
    df["queryPair"]=df["startQuery"].astype(str)+":"+df['endQuery'].astype(str)

    df[["parent","notID"]]=df["useful_info"].str.extract('ID=(.+?)\;(.*)',expand=True)
    parentDF=df.dropna(axis=0) # this shall be used for pulling information from links

    df=df.groupby('queryPair', group_keys=False).apply(lambda x: x.ix[x.ord.idxmin()]).reset_index(drop=True) # this might contain links and non links. these need to be separated and processed separately. For links pull information from Non-links and for others just calculate from there
    linksDF=pd.DataFrame([])
    if len(df)>0:
        df["link"]=df["useful_info"].str.extract('Parent=(.*)',expand=False)
        linksDF=df.dropna(subset=["link"],axis=0)# contains all which need to be linked

    readyDF=df[~df.index.isin(linksDF.index)] # contains all which can do not need anything else done to them except extracting info
    readyDF=readyDF.replace({'ord':reverseOrder})
    if len(readyDF)>0:
        readyDF["gene_name"]=readyDF["notID"].str.extract('gene_name=(.+?)\;',expand=False)
        readyDF['hum_nearest_SS']=readyDF['ord']+":"+readyDF['gene_name']
        readyDF.drop(['remove1',
                      'ord',
                      'useful_info',
                      'remove2',
                      'queryPair',
                      'parent',
                      'notID',
                      'link',
                      'gene_name'],inplace=True,axis=1)

    # now lets work on the links
    setLink=set(linksDF['link'])
    parentDF=parentDF[parentDF["parent"].isin(setLink)]
    parentDF=parentDF.groupby('queryPair', group_keys=False).apply(lambda x: x.ix[x.ord.idxmin()]).reset_index(drop=True)
    parentDF=parentDF.replace({'ord':reverseOrder})
    if len(parentDF)>0:
        parentDF["gene_name"]=parentDF["notID"].str.extract('gene_name=(.+?)\;',expand=False)
        parentDF['hum_nearest_SS']=parentDF['ord']+":"+parentDF['gene_name']
        parentDF.drop(['remove1',
                      'ord',
                      'useful_info',
                      'remove2',
                      'queryPair',
                      'parent',
                      'notID',
                      'gene_name'],inplace=True,axis=1)

    finalBed=pd.DataFrame([])
    if len(parentDF)>0 and len(readyDF)>0:
        finalBed=pd.concat([readyDF,parentDF]).reset_index(drop=True)
    if len(parentDF)>0 and len(readyDF)==0:
        finalBed=parentDF.copy(deep=True)
    if len(parentDF)==0 and len(readyDF)>0:
        finalBed=readyDF.copy(deep=True)

    if len(finalBed)==0:
        return

    finalBed.drop(["startRegionPos",
                   "endRegionPos",
                   "something1",
                   "something2",
                   "something3",
                   "distance"],axis=1,inplace=True)
    finalBed.columns=['chr',
                      'G1_RS',
                      'G1_RE',
                      'hum_nearest_SS']

    finalDF=pd.DataFrame(pd.merge(data,finalBed,on=['chr',
                                                      'G1_RS',
                                                      'G1_RE'],how='inner')).reset_index(drop=True)
    # data1=pd.DataFrame([])
    data=data[~(data['comb'].isin(set(finalDF["comb"])))].copy()
    data.reset_index(drop=True,inplace=True)
    data1=data.copy()
    data1['hum_nearest_SS']="-"

    return pd.concat([finalDF,data1]).reset_index(drop=True)

#this function should find the closest known acceptor/donor for each of the final records
# we could also implement passing the known donor acceptor lists as parameters
# we could also score, if required by how close the alignment is to the known donor/acceptor

# should be done right before the results cleanup

# also makes sense to keep G2 alignment information
# needs to be kept along with the best sequence information in the tuple

# What needs to be done:
# this has to be done separately based on which side of the sequence the G2 is aligned
# If aligned first half - search for donors and acceptors to the end of the hiv alignment
# if aligned second half of the read - search for donors and acceptors to the beginning of the hiv

# need to read what is meant by donor and acceptor
def donorAcceptor(dataPos):
    donors=[("D1b",725),
            ("D1",742),
            ("D1c",746),
            ("D1a",4720),
            ("D2",4961),
            ("D3",5462),
            ("D4",6043),
            ("D5",6729),
            ("D6",8422),
            ("End",9713)]
    acceptors=[("Beg",1),
            ("A1a",4542),
            ("A1",4912),
            ("A2",5389),
            ("A3",5776),
            ("A4c",5935),
            ("A4a",5953),
            ("A4b",5959),
            ("A5",5975),
            ("A5a",5979),
            ("A5b",5982),
            ("A6",6602),
            ("A6a",6654),
            ("A7a",8334),
            ("A7b",8340),
            ("A7c",8344),
            ("A7d",8362),
            ("A7",8368),
            ("A7e",8381),
            ("A7f",8485),
            ("A8d",9128),
            ("A8a",9144),
            ("A8c",9156),
            ("A8",9164),
            ("A8b",9177),
            ("A8e",9183),
            ("A8f",9221)]

    donorsDF=pd.DataFrame(donors)
    acceptorsDF=pd.DataFrame(acceptors)

    def findClosest(val,df):
        df['diff']=abs(df[1]-val)
        return df.sort_values(by="diff",ascending=True).reset_index(drop=True).iloc[0][0]

    dataPos.reset_index(drop=True,inplace=True)
    dataPos["nearestDonor"]=dataPos.apply(lambda row: findClosest(row["G2"]),axis=1)
    dataPos["nearestAcceptor"]=dataPos.apply(lambda row: findClosest(row["G2"]),axis=1)

    return dataPos

def rest(dataPos,args,data,unpaired,baseName,outDir,dirPath):
    if not unpaired:
        dataPos=addSpan(data,dataPos)
        dataPos.loc[dataPos['spanCount'].isnull(),['spanCount']]=dataPos.loc[dataPos['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
        dataPos.loc[dataPos['spanR1-R2'].isnull(),['spanR1-R2']]=dataPos.loc[dataPos['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())

    dataBed=dataPos[['chr','G1_RS','G1_RE']].drop_duplicates()
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

            dataPos=score(dataPos,args,args.minLen)
            dataPos=dataPos.sort_values(by='score',ascending=False).reset_index(drop=True)
            dataPos[['seq','hum_pos','drop','overlap','gap']]=dataPos['seq'].apply(pd.Series)
            dataPos.drop("drop",axis=1,inplace=True)

            dataPos.to_csv(os.path.abspath(args.out)+".full.csv",index=False)
            dataPosClean=dataPos[(dataPos['entropyScore_g2']>args.minEntropy) \
                                &(dataPos['entropyScore_g1']>args.minEntropy) \
                                &(dataPos['score']>args.score)]

            colsOrder=["hum_nearest_SS",
                       "chr",
                       "hum_pos",
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
                dataPosClean["fileName"]=dataPosClean['hum_nearest_SS'].str.strip('\n')+"@"+dataPosClean['comb'].str.split("@",expand=True)[0]+".fa"
            else:
                dataPosClean["fileName"]=dataPosClean['hum_nearest_SS'].str.strip('\n')+"@"+dataPosClean["R"]+"@"+dataPosClean['comb'].str.split("@",expand=True)[0]+"_R1.all.fa"

            dataPosClean[colsOrder].to_csv(os.path.abspath(args.out)+".csv",index=False)

            # completely forgot that we can write reads from the sam file
            # should make it much faster
            # for each comb in the Pos.csv file
            # find corresponding positions in the original human file
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

# def child(path,outDir,covRange,numReps,sequenceRef,annotationRef,threads,cont):
#     baseDirName = path.split("/")[-1].split(".")[:-1][0]
#     finDir = outDir+"/"+baseDirName+"/"+baseDirName
#     if not cont==False:
#         # If this method with the script works - consider writing the config file from here rather than passing parameters
#         for scalefactor in xfrange(cont,covRange[1],covRange[2]):
#             randSeed = random.sample(range(1,10000),numReps) # Need to think where to place this line
#             for rep in range(numReps):
#                 scriptCMD = "./rnaseq_al_pipe.sh "+path+" "+finDir+" "+str(randSeed[rep]+scalefactor)+" "+str(rep)+" "+sequenceRef+" "+annotationRef+" "+str(threads)
#                 os.system(scriptCMD)
#     else:
#         # If this method with the script works - consider writing the config file from here rather than passing parameters
#         for scalefactor in xfrange(covRange[0],covRange[1],covRange[2]):
#             randSeed = random.sample(range(1,10000),numReps) # Need to think where to place this line
#             for rep in range(numReps):
#                 scriptCMD = "./rnaseq_al_pipe.sh "+path+" "+finDir+" "+str(randSeed[rep]+scalefactor)+" "+str(rep)+" "+sequenceRef+" "+annotationRef+" "+str(threads)
#                 os.system(scriptCMD)
#     os.remove(finDir+".sam")
#     os._exit(0)

# # multithreaded version of the algorithm
# def wrapperThreaded(outDir,baseName,dirPath,fileName,minLen,end,args):
#     dataG2=pd.read_csv(os.path.abspath(args.input1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
#                                                                                                                     'FLAG',
#                                                                                                                     'RNAME',
#                                                                                                                     'POS',
#                                                                                                                     'MAPQ',
#                                                                                                                     'CIGAR',
#                                                                                                                     'RNEXT',
#                                                                                                                     'PNEXT',
#                                                                                                                     'TLEN',
#                                                                                                                     'SEQ',
#                                                                                                                     'QUAL'])
#     dataG1=pd.read_csv(os.path.abspath(args.input2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
#                                                                                                                     'FLAG',
#                                                                                                                     'RNAME',
#                                                                                                                     'POS',
#                                                                                                                     'MAPQ',
#                                                                                                                     'CIGAR',
#                                                                                                                     'RNEXT',
#                                                                                                                     'PNEXT',
#                                                                                                                     'TLEN',
#                                                                                                                     'SEQ',
#                                                                                                                     'QUAL'])

#     if (len(dataG2)==0 or len(dataG1)==0): #exit if either alignment is empty
#         return
    
#     outDirPOS=outDir+"/Positions/"
#     if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
#         os.mkdir(os.path.abspath(outDir+"/Positions/"))
#     if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
#         os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

#     if len(dataG2)>0 and len(dataG1)>0:
#         # first calculate how many parts the dataFrame will have
#         # also, should preprocess dataG2 as much as possible before parallelization
#         # dont forget about unpaired option.
#         extractFlagBits(dataG2)
#         extractFlagBits(dataG1)

#         childPIDS = []
#         def parent():
#             for dfPart in dfParts:
#                 if len(childPIDS) >= tf[0]:
#                         childPIDS[0].join()
#                         childPIDS.remove(childPIDS[0])
#                 else:
#                     p = multiprocessing.Process(target=child, args=(dfPart,))
#                     childPIDS.append(p)
#                     p.start()

#             while(len(childPIDS) > 0):
#                 childPIDS[-1].join()
#                 childPIDS.remove(childPIDS[-1])

def wrapper(outDir,baseName,dirPath,fileName,minLen,end,args):
    # load data from local alignments
    dataG2=pd.read_csv(os.path.abspath(args.input1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
                                                                                                                    'FLAG',
                                                                                                                    'RNAME',
                                                                                                                    'POS',
                                                                                                                    'MAPQ',
                                                                                                                    'CIGAR',
                                                                                                                    'RNEXT',
                                                                                                                    'PNEXT',
                                                                                                                    'TLEN',
                                                                                                                    'SEQ',
                                                                                                                    'QUAL'])
    dataG1=pd.read_csv(os.path.abspath(args.input2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
                                                                                                                    'FLAG',
                                                                                                                    'RNAME',
                                                                                                                    'POS',
                                                                                                                    'MAPQ',
                                                                                                                    'CIGAR',
                                                                                                                    'RNEXT',
                                                                                                                    'PNEXT',
                                                                                                                    'TLEN',
                                                                                                                    'SEQ',
                                                                                                                    'QUAL'])

    if (len(dataG2)==0 or len(dataG1)==0): #exit if either alignment is empty
        return
    
    outDirPOS=outDir+"/Positions/"
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
        os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

    if len(dataG2)>0 and len(dataG1)>0:
        global reportDF
        dataLen=dataG1["SEQ"].str.len()
        reportDF["number of reads before removing spliced reads"]=len(dataG1)
        reportDF["original read len mean"]=dataLen.mean()
        reportDF["original read len std"]=dataLen.std()
        reportDF["original read len min"]=dataLen.min()
        reportDF["original read len max"]=dataLen.max()
        # extract flag information
        extractFlagBits(dataG2)
        extractFlagBits(dataG1)
        unpaired=False
        if len(dataG2[dataG2["paired"]>0])==0 and len(dataG1[dataG1["paired"]>0])==0:
            unpaired=True
        
        # need to verify that both paired and unpaired work below.
        # i guess for paired we need to consider the case when splicing occurs on one mate, and thus we only want to remove that mate
        # but not the other mate, since that might still contain a chimeric site, provided it is not an end-to-tnd alignment

        if args.spliced is not None:
            if os.path.exists(os.path.abspath(args.spliced)):
                dataSpliced=pd.read_csv(os.path.abspath(args.spliced),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
                                                                                                                                   'FLAG',
                                                                                                                                   'RNAME',
                                                                                                                                   'POS',
                                                                                                                                   'MAPQ',
                                                                                                                                   'CIGAR',
                                                                                                                                   'RNEXT',
                                                                                                                                   'PNEXT',
                                                                                                                                   'TLEN',
                                                                                                                                   'SEQ',
                                                                                                                                   'QUAL'])
                # dataSpliced=dataSpliced[dataSpliced['CIGAR'].str.contains("N")]
                extractFlagBits(dataSpliced)
                dataSpliced["tid"]=dataSpliced['QNAME']+dataSpliced['firstRead'].astype(str)+dataSpliced['lastRead'].astype(str)
                dataG1["tid"]=dataG1['QNAME']+dataG1['firstRead'].astype(str)+dataG1['lastRead'].astype(str)
                dataG2["tid"]=dataG2['QNAME']+dataG2['firstRead'].astype(str)+dataG2['lastRead'].astype(str)
                dataG2=dataG2[~dataG2['tid'].isin(set(dataSpliced['tid']))]
                dataG1=dataG1[dataG1['tid'].isin(set(dataG2['tid']))]
                dataG2=dataG2[dataG2['tid'].isin(set(dataG1['tid']))]
                dataG2.drop(['tid'],axis=1,inplace=True)
                dataG1.drop(['tid'],axis=1,inplace=True)

                # global reportDF
                dataLen=dataG1["SEQ"].str.len()
                reportDF["number of reads shared between G2 and human post splicing"]=len(dataG1)
                reportDF["read len mean shared between G2 and human post splicing"]=dataLen.mean()
                reportDF["read len std shared between G2 and human post splicing"]=dataLen.std()
                reportDF["read len min shared between G2 and human post splicing"]=dataLen.min()
                reportDF["read len max shared between G2 and human post splicing"]=dataLen.max()
            else:
                print("Spliced file does not exist")
        elif unpaired==False: # still remove all human which are not in hiv
            dataG1["tid"]=dataG1['QNAME']+dataG1['firstRead'].astype(str)+dataG1['lastRead'].astype(str)
            dataG2["tid"]=dataG2['QNAME']+dataG2['firstRead'].astype(str)+dataG2['lastRead'].astype(str)
            dataG2=dataG2[dataG2['tid'].isin(set(dataG1['tid']))]
            dataG1=dataG1[dataG1['tid'].isin(set(dataG2['tid']))]
            dataG2.drop(['tid'],axis=1,inplace=True)
            dataG1.drop(['tid'],axis=1,inplace=True)
            # global reportDF
            dataLen=dataG1["SEQ"].str.len()
            reportDF["number of reads shared between G2 and human post splicing"]=len(dataG1)
            reportDF["read len mean shared between G2 and human post splicing"]=dataLen.mean()
            reportDF["read len std shared between G2 and human post splicing"]=dataLen.std()
            reportDF["read len min shared between G2 and human post splicing"]=dataLen.min()
            reportDF["read len max shared between G2 and human post splicing"]=dataLen.max()
        else:
            dataG1["tid"]=dataG1['QNAME']
            dataG2["tid"]=dataG2['QNAME']
            dataG2=dataG2[dataG2['tid'].isin(set(dataG1['tid']))]
            dataG1=dataG1[dataG1['tid'].isin(set(dataG2['tid']))]
            dataG2.drop(['tid'],axis=1,inplace=True)
            dataG1.drop(['tid'],axis=1,inplace=True)
            # global reportDF
            dataLen=dataG1["SEQ"].str.len()
            reportDF["number of reads shared between G2 and human post splicing"]=len(dataG1)
            reportDF["read len mean shared between G2 and human post splicing"]=dataLen.mean()
            reportDF["read len std shared between G2 and human post splicing"]=dataLen.std()
            reportDF["read len min shared between G2 and human post splicing"]=dataLen.min()
            reportDF["read len max shared between G2 and human post splicing"]=dataLen.max()

        # extract start and end for both template and reference
        dataG2=extractStartEnd(dataG2)
        dataG1=extractStartEnd(dataG1)

        dataG1,dataG2=filterReads(dataG1,dataG2)
        # global reportDF
        dataLen=dataG1["SEQ"].str.len()
        reportDF["number of reads after sam flag filtering"]=len(dataG1)
        reportDF["read len mean after sam flag filtering"]=dataLen.mean()
        reportDF["read len std after sam flag filtering"]=dataLen.std()
        reportDF["read len min after sam flag filtering"]=dataLen.min()
        reportDF["read len max after sam flag filtering"]=dataLen.max()
        if len(dataG1)==0:
            return

        data=pd.DataFrame(dataG2["QNAME"]).reset_index(drop=True)
        dataG1=dataG1.reset_index(drop=True)
        dataG2=dataG2.reset_index(drop=True)
        if unpaired:
            data=createDataUnpaired(data,dataG1,dataG2)
        else:
            data=createData(data,dataG1,dataG2)

        del dataG1
        del dataG2
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
            # print(d[d["QNAME"]=="NB501749:128:HFHJVBGX3:2:13304:6468:2639"])
        else:
            dataPos=filterOverlapCombine(data,args)
        dataPos=dataPos[(dataPos["entropyScore_g1"]>args.minEntropy)&(dataPos["entropyScore_g2"]>args.minEntropy)]
        # global reportDF
        # dataLen=dataPos["SEQ"].str.len()
        # reportDF["number of reads after minimum entropy cutoff"]=len(dataPos)
        # reportDF["read len mean after minimum entropy cutoff"]=dataLen.mean()
        # reportDF["read len std after minimum entropy cutoff"]=dataLen.std()
        # reportDF["read len min after minimum entropy cutoff"]=dataLen.min()
        # reportDF["read len max after minimum entropy cutoff"]=dataLen.max()
        dataPos=findSupport(dataPos,minLen,unpaired)
        if len(dataPos)>0:
            rest(dataPos,args,data,unpaired,baseName,outDir,dirPath)
    
    return 1

def main(args):
    end=""
    if args.end==True:
        end='.no_dup'

    outDir="/".join(os.path.abspath(args.out).split("/")[:-1])
    if not os.path.exists(outDir):
        print('output directory does not exist')
        return

    al1=os.path.abspath(args.input1)
    al2=os.path.abspath(args.input2)
    fullPath1=os.path.abspath(os.path.realpath(al1))
    fullPath2=os.path.abspath(os.path.realpath(al2))

    if os.path.exists(fullPath1) and os.path.exists(fullPath2):
        fileName=fullPath1.split('/')[-1]
        dirPath="/".join(fullPath1.split('/')[:-1])
        baseName=fileName.split(".")[0]
        ext=".".join(fileName.split(".")[1:-1])

        resultsRow=wrapper(outDir,baseName,dirPath,fileName,args.minLen,end,args)
        if args.quiet:
            printReport()

        global reportDF
        reportDF.to_csv(os.path.abspath(args.out)+".report",index=False)

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
#./hiv.py -i1 hiv.sam -i2 hum.sam -o ${outputDir}_R2/${sample}${baseEnd} -t 12 --minLen 30 -a ${annotation} --overlap 5 --gap 5
    parser.add_argument('-i1',
                              '--input1',
                              required=True,
                              type=str,
                              help="first alignment")
    parser.add_argument('-i2',
                              '--input2',
                              required=True,
                              type=str,
                              help="second alignment")
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
                              default=0.75,
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
                              help="annotation for the human genome")
    parser.add_argument('-e',
                              '--end',
                              default="",
                              type=str,
                              help="suffix to append to the end of the ouput name")
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
    parser.add_argument('--spliced',
                              required=False,
                              type=str,
                              help="spliced end-to-end alignment. Reads will be subtracted from the main alignments and will not be reported in the final report of integrations sites")
    parser.set_defaults(func=main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    chimFinder(sys.argv[1:])

# at the very end need to build python environment and requirements.txt
# since the code should work with both python2 and python3 need to do this for both

# implement writeReads function without grepping

# What really needs to be implemented is the following function
# The function will take as an input all weights and scoring arguments passed
# and it will calculate the minimum overall score for the desired combination
# unless the minimum score is specified, in which case it should override the calculation

# In the end the fasta records will have to be written from the dataframe,
# since the raw sequence data is not supplied to the application and for that matter should not be supplied
# This should not be a priority however at the moment.

# Should we consider the number of insertions/deletions/mismatches in the score calculations?

# Another idea for span reads:
# 1. Add only those that are fully hiv on one side and fully human on the other side
# should not contain a fragmented G2 HUMan read on either side
# it could be fragmented if it appears to indicate the same splice site

# Supply human splicing annotation to hisat and do not include novel splice sites
# use the same as the one used by the UCSC BLAT

# question: when provided with known splice sites, does hisat2 still report novel ones or only align in accordance with those provided

# perhaps the ITGAE gene is in the 1B6 sample, but is not detected in the annotation,
# because we are only annotating the 3'SS. Try adding Ellas script to the output once again,
# to see if anything can be revealed through that script
# for now will check manually
# the manual analysis has shown that the position is nowhere to be found in the results page

# It would be much more efficient to not hold read names and sequences in the memory and to not add them to the final results
# Instead It makes sense to output a bed file with positions
# bedtools intersect could then be used to show the actual reads in the alignments
# An option could be preserved to report read names and sequences as they are now

# workflow for building the graph when begin developing C++ code for the software:
# 0. Perhaps a network flow would work well for the task?
# 1. create a disjoint graph of all reads from G2 (whichever alignment has fewer reads)
# 2. when parsing through the second alignment, populate corresponding nodes with information
# 3. when performing collapsing - create edges between nodes

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Could it be that the issue which manifests itself in that one alignment is a subset of another
# like the majority of outHIV_HCV results
# could it be that it is due to not taking into account CIGAR string deletions and insertions properly?
# Need to check and figure this out with confidence once and for all
# Either it is something else or it further increases sensitivity of the results/method!
# Also try viewing alignments in IGV
# This might help understand the interpretation of the cigar string better
# Perhaps this will uncover what is missing from the extractStartEnd procedure
# I really need to figure this one out!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Also remove all combs from comb column
# replace them with the central number
# verify how this approximation works

# To do so, extract the final reads from the sam files in HIV_HCV experiment
# This will yield a much smaller set of data to work with
# Analyze the CIGAR string composition of these records

#PARALLELIZATION!!!!!!!!!!!!!!!!!!!!!!!!
# Should be very easy to achive
# simply split the human dataframe into n chunks
# run the rest of the analysis in parallel on each chunk
# merge resulting dataframes together
# sort by score
# save

# specify donors/acceptors as a command line argument
# calculate and report the offset from the nearest donor/acceptor on G2

# Perhaps would be better to read sam/bam in lines and perform all initial calculations and filtering then
# For instance
# Begin parsing G2 BAM/SAM (note that by using pysam we can work with both BAM and SAM files)
# have a cascade of filters (increasing order of computational complexity:
# if end-to-end - pass to the next read
# if len alignment less than threshold (30bp) - pass
# etc
# Then begin parsing human in the same manner, this time:
# if not in G2 set - pass
# if end-to-end - remove from G2-set and pass
# etc.
# The rest of the calculations can be done afterwards

# So here is the problem with the report.
# After data is grouped the number of reads retained must be computed not as length of the dataframe, but rather the sum of counts in the dataFrame
# additionally, sequence length must be estimated from the tuples that belong to each record
# that is going to complicate things a bit

# try adopting a sliding window approach to entropy calculation
# calculate entropy for each kmer of size (minimum alignment length) - report highest score
# also compute entropy of the full sequence and if higher than that of individual kmers then report full
# this will ensure that at least part of the sequence conforms to the desired entropy
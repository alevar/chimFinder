#!/usr/bin/env python

import pandas as pd
import numpy as np
import numba
import math
import os
import subprocess
import multiprocessing
import signal
import shutil
import sys
import glob
import itertools
import argparse
import scipy
import warnings
from pybedtools import BedTool
warnings.filterwarnings('ignore')

# Find chimeric reads (even if < 31nt alignment length to hiv or hum) from full alignments
def getFromSplits(dataHIV,dataHUM):
    cleanHUM=dataHUM[(dataHUM["paired"]==1)&(dataHUM["secondaryAlignment"]==0)&(dataHUM["aligned2Mates"]==0)&(dataHUM["unmappedCurr"]==0)&(dataHUM["unmappedMate"]==8)]
    cleanHUM_First=cleanHUM[(cleanHUM["firstRead"]==64)&(cleanHUM["lastRead"]==0)]
    cleanHUM_Last=cleanHUM[(cleanHUM["firstRead"]==0)&(cleanHUM["lastRead"]==128)]

    cleanHIV=dataHIV[(dataHIV["paired"]==1)&(dataHIV["secondaryAlignment"]==0)&(dataHIV["aligned2Mates"]==0)&(dataHIV["unmappedCurr"]==0)&(dataHIV["unmappedMate"]==8)]
    cleanHIV_First=cleanHIV[(cleanHIV["firstRead"]==64)&(cleanHIV["lastRead"]==0)]
    cleanHIV_Last=cleanHIV[(cleanHIV["firstRead"]==0)&(cleanHIV["lastRead"]==128)]
    #find mates within paired-end reads which aligned to both hiv and hum
    setHIV_First=set(cleanHIV_First["QNAME"])
    setHUM_First=set(cleanHUM_First["QNAME"])
    splitMates=setHIV_First.intersection(setHUM_First)

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

# mark reads that have HIV on the right side
def leftRight(data,minLen):
    data["HIV"]=""
    data["hivR1"]=""
    data["hivR2"]=""
    data["sepR1"]=""
    data["sepR2"]=""

    data.loc[~(data["R1HUM_ID"]==0)&~(data["R1HIV_ID"]==0)&(data["R1HIV_TS"]-data["R1HUM_TS"]>minLen)&(data["R1HIV_TE"]-data["R1HUM_TE"]>minLen),"hivR1"]="R1:r"
    data.loc[~(data["R1HUM_ID"]==0)&~(data["R1HIV_ID"]==0)&(data["R1HUM_TS"]-data["R1HIV_TS"]>minLen)&(data["R1HUM_TE"]-data["R1HIV_TE"]>minLen),"hivR1"]="R1:l"

    data.loc[~(data["R2HUM_ID"]==0)&~(data["R2HIV_ID"]==0)&(data["R2HIV_TS"]-data["R2HUM_TS"]>minLen)&(data["R2HIV_TE"]-data["R2HUM_TE"]>minLen),"hivR2"]="R2:r"
    data.loc[~(data["R2HUM_ID"]==0)&~(data["R2HIV_ID"]==0)&(data["R2HUM_TS"]-data["R2HIV_TS"]>minLen)&(data["R2HUM_TE"]-data["R2HIV_TE"]>minLen),"hivR2"]="R2:l"

    data.loc[(~(data["R1HUM_ID"].astype(str)=="0")&~(data["R2HIV_ID"].astype(str)=="0")&(data["R2HUM_ID"].astype(str)=="0")), "sepR1"]="sepR1:2"
    data.loc[(~(data["R2HUM_ID"].astype(str)=="0")&~(data["R1HIV_ID"].astype(str)=="0")&(data["R1HUM_ID"].astype(str)=="0")), "sepR2"]="sepR2:1"

    data["HIV"]=data["hivR1"]+data["hivR2"]+data["sepR1"]+data["sepR2"]
    data=data[data['HIV'].str.len()>0]
    return data.drop(["hivR1","hivR2","sepR1","sepR2"],axis=1)

def leftRightUnpaired(data,minLen):
    data["HIV"]=""

    data.loc[~(data["HUM_ID"]==0)&~(data["HIV_ID"]==0)&(data["HIV_TS"]-data["HUM_TS"]>minLen)&(data["HIV_TE"]-data["HUM_TE"]>minLen),"HIV"]="r"
    data.loc[~(data["HUM_ID"]==0)&~(data["HIV_ID"]==0)&(data["HUM_TS"]-data["HIV_TS"]>minLen)&(data["HUM_TE"]-data["HIV_TE"]>minLen),"HIV"]="l"

    data=data[data['HIV'].str.len()>0]
    return data

def overlap(data):
    data.reset_index(inplace=True,drop=True)
    data["overlapR1"]=data[['R1HUM_TE','R1HIV_TE']].min(axis=1)-data[['R1HUM_TS','R1HIV_TS']].max(axis=1) #min end - max start
    data["gapR1"]=0
    data.loc[data['overlapR1']<0,['gapR1']]=data.loc[data['overlapR1']<0,['overlapR1']]["overlapR1"].abs()
    data.loc[data['overlapR1']<0,['overlapR1']]=0
    data["overlapR2"]=data[['R2HUM_TE','R2HIV_TE']].min(axis=1)-data[['R2HUM_TS','R2HIV_TS']].max(axis=1)
    data["gapR2"]=0
    data.loc[data['overlapR2']<0,['gapR2']]=data.loc[data['overlapR2']<0,['overlapR2']]["overlapR2"].abs()
    data.loc[data['overlapR2']<0,['overlapR2']]=0
    return data

def overlapUnpaired(data):
    data.reset_index(inplace=True,drop=True)
    data["overlap"]=data[['HUM_TE','HIV_TE']].min(axis=1)-data[['HUM_TS','HIV_TS']].max(axis=1) #min end - max start
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

    data["READ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$").replace(np.nan,0).astype(int)
    data["END"]=data.READ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[S]").replace(np.nan,0).astype(int)

    data16=data[data["reversedCurr"]==16]
    data0=data[data["reversedCurr"]==0]
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
def filterReads(dataHUM,dataHIV):
    #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
    dataHUM=dataHUM[(dataHUM["secondaryAlignment"]==0)&(dataHUM["PCRdup"]==0)&(dataHUM["suppAl"]==0)&(dataHUM["noPassFilter"]==0)]
    dataHIV=dataHIV[(dataHIV["secondaryAlignment"]==0)&(dataHIV["PCRdup"]==0)&(dataHIV["suppAl"]==0)&(dataHIV["noPassFilter"]==0)]
    return dataHUM, dataHIV

def createData(data,dataHUM,dataHIV):
    dataHUM=dataHUM[dataHUM['QNAME'].isin(set(data['QNAME']))]
    dataHUMR1=pd.DataFrame([])
    dataHUMR1[["QNAME",
                "R1HUM_TS",
                "R1HUM_TE",
                "R1HUM_ID",
                "R1HUM_RS",
                "R1HUM_RE",
                "READ_LEN_R1",
                "R1HUM_SEQ",
                "QUAL_R1",
                "R1HUM_MAPQ",
                "R1HUM_reversedCurr"]]=dataHUM[dataHUM['firstRead']==64][["QNAME",
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
    dataHUMR2=pd.DataFrame([])
    dataHUMR2[["QNAME",
                "R2HUM_TS",
                "R2HUM_TE",
                "R2HUM_ID",
                "R2HUM_RS",
                "R2HUM_RE",
                "READ_LEN_R2",
                "R2HUM_SEQ",
                "QUAL_R2",
                "R2HUM_MAPQ",
                "R2HUM_reversedCurr"]]=dataHUM[dataHUM['lastRead']==128][["QNAME",
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
    dataHIVR1=pd.DataFrame([])
    dataHIVR1[["QNAME",
                "R1HIV_TS",
                "R1HIV_TE",
                "R1HIV_ID",
                "R1HIV_RS",
                "R1HIV_RE",
                "R1HIV_SEQ",
                "R1HIV_MAPQ",
                "R1HIV_reversedCurr"]]=dataHIV[dataHIV['firstRead']==64][["QNAME",
                                                                            "Template_start",
                                                                            "Template_end",
                                                                            "RNAME",
                                                                            "Reference_start",
                                                                            "Reference_end",
                                                                            "SEQ",
                                                                            "MAPQ",
                                                                            "reversedCurr"]]
    dataHIVR2=pd.DataFrame([])
    dataHIVR2[["QNAME",
                "R2HIV_TS",
                "R2HIV_TE",
                "R2HIV_ID",
                "R2HIV_RS",
                "R2HIV_RE",
                "R2HIV_SEQ",
                "R2HIV_MAPQ",
                "R2HIV_reversedCurr"]]=dataHIV[dataHIV['lastRead']==128][["QNAME",
                                                                        "Template_start",
                                                                        "Template_end",
                                                                        "RNAME",
                                                                        "Reference_start",
                                                                        "Reference_end",
                                                                        "SEQ",
                                                                        "MAPQ",
                                                                        "reversedCurr"]]
    data=data.merge(dataHUMR1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataHUMR2,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataHIVR1,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dataHIVR2,left_on='QNAME', right_on='QNAME', how='left')
    data[["R1HUM_ID",
        "R2HUM_ID",
        "R1HIV_ID",
        "R2HIV_ID",
        "R1HUM_MAPQ",
        "R2HUM_MAPQ",
        "R1HIV_MAPQ",
        "R2HIV_MAPQ",
        "R1HUM_reversedCurr",
        "R2HUM_reversedCurr",
        "R1HIV_reversedCurr",
        "R2HIV_reversedCurr"]].fillna('',inplace=True)
    data.fillna(0,inplace=True)
    return data

def createDataUnpaired(data,dataHUM,dataHIV):
    dataHUM=dataHUM[dataHUM['QNAME'].isin(set(data['QNAME']))]
    dfHUM=pd.DataFrame([])
    dfHIV=pd.DataFrame([])
    dfHUM[["QNAME",
            "HUM_TS",
            "HUM_TE",
            "HUM_ID",
            "HUM_RS",
            "HUM_RE",
            "READ_LEN",
            "HUM_SEQ",
            "QUAL",
            "HUM_MAPQ",
            "HUM_reversedCurr"]]=dataHUM[["QNAME",
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
    dfHIV[["QNAME",
            "HIV_TS",
            "HIV_TE",
            "HIV_ID",
            "HIV_RS",
            "HIV_RE",
            "HIV_SEQ",
            "HIV_MAPQ",
            "HIV_reversedCurr"]]=dataHIV[["QNAME",
                                            "Template_start",
                                            "Template_end",
                                            "RNAME",
                                            "Reference_start",
                                            "Reference_end",
                                            "SEQ",
                                            "MAPQ",
                                            "reversedCurr"]]
    data=data.merge(dfHUM,left_on='QNAME', right_on='QNAME', how='left')
    data=data.merge(dfHIV,left_on='QNAME', right_on='QNAME', how='left')
    
    data[["HUM_ID",
        "HIV_ID",
        "HUM_MAPQ",
        "HIV_MAPQ",
        "HUM_reversedCurr",
        "HIV_reversedCurr"]].fillna('',inplace=True)
    data.fillna(0,inplace=True)
    return data

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

    maxSubstringLen=findCF(len(s))
    n=countSubString(s,maxSubstringLen)
    return math.log(n,4)/maxSubstringLen

# calculate the expected entropy
def expectedEntropy(n):
    fourN=float(4**n)
    logTerm=(fourN-(fourN*((1.0-(1.0/fourN))**fourN)))
    return math.log(logTerm,4)/float(n)

# How to compute the mimimum entropy for a string of a given length and characters

#consider using technique described in https://academic.oup.com/bioinformatics/article/27/8/1061/227307/Topological-entropy-of-DNA-sequences
#another paper: Look at figure 2 https://www.nature.com/articles/srep19788#f2
#more https://www.xycoon.com/normalized_entropy.htm

def processAligns(seqHum,seqHIV,qual,i1_hiv,i2_hiv,i1_hum,i2_hum,readLen,rHUM,rHIV):

    entropyScore_hiv=0
    entropyScore_hum=0
    meanQual_hiv=0
    meanQual_hum=0

    s_hiv=""
    if rHIV==True: # if reverse complemented take into account
        s_hiv=seqHIV[readLen-i2_hiv:readLen-i1_hiv]
    else:
        s_hiv=seqHIV[i1_hiv:i2_hiv]
    if not len(s_hiv)==0:
        entropyScore_hiv=topologicalNormalizedEntropy(s_hiv)
        q_hiv=qual[i1_hiv:i2_hiv]
        if len(q_hiv)>0:
            meanQual_hiv=sum([ord(x)-33 for x in q_hiv])/len(q_hiv)
        else:
            meanQual_hiv=0

    s_hiv=""
    if rHUM==True: # if reverse complemented take into account
        s_hum=seqHum[readLen-i2_hum:readLen-i1_hum]
    else:
        s_hum=seqHum[i1_hum:i2_hum]
    if not len(s_hum)==0:
        entropyScore_hum=topologicalNormalizedEntropy(s_hum)
        q_hum=qual[i1_hum:i2_hum]
        if len(q_hum)>0:
            meanQual_hum=sum([ord(x)-33 for x in q_hum])/len(q_hum)
        else:
            meanQual_hum=0

    return pd.Series({"entropyScore_hiv":entropyScore_hiv,
                        "meanQual_hiv":meanQual_hiv,
                        "entropyScore_hum":entropyScore_hum,
                        "meanQual_hum":meanQual_hum})


# this function filters by overlap and flanking alignment length
# further it removes unnecessary data and combines into a single full dataframe with unified naming avoiding R1/R2 conventions
# output can be saved as .full.csv and then grouped by split position all at once

# this function allows identifying best reads
# However, since greater minLen values for HIV and HUM alignments will yield fewer reads but at higher confidence
# support reads could ignore the min len requirement as long as the split position is identical
# or perhaps the minLen requirement for the support reads should be lower
def filterOverlapCombine(data,args):
    dropList=["R1HUM_TS",
              "R1HUM_TE",
              "R1HUM_ID",
              "R1HUM_RS",
              "R1HUM_RE",
              "R2HUM_TS",
              "R2HUM_TE",
              "R2HUM_ID",
              "R2HUM_RS",
              "R2HUM_RE",
              "R1HIV_TS",
              "R1HIV_TE",
              "R1HIV_ID",
              "R1HIV_RS",
              "R1HIV_RE",
              "R2HIV_TS",
              "R2HIV_TE",
              "R2HIV_ID",
              "R2HIV_RS",
              "R2HIV_RE",
              "overlapR1",
              "overlapR2",
              "gapR1",
              "gapR2",
              "R1HUM_SEQ",
              "R2HUM_SEQ",
              "R1HIV_SEQ",
              "R2HIV_SEQ",
              "QUAL_R1",
              "QUAL_R2",
              "R1HUM_MAPQ",
              "R1HIV_MAPQ",
              "R2HUM_MAPQ",
              "R2HIV_MAPQ",
              "READ_LEN_R1",
              "READ_LEN_R2",
              "R1HUM_reversedCurr",
              "R2HUM_reversedCurr",
              "R1HIV_reversedCurr",
              "R2HIV_reversedCurr"]

    # R1 right
    data['entropyScore_hiv']=0
    data['meanQual_hiv']=0
    data['mapQual_hiv']=0
    data['entropyScore_hum']=0
    data['meanQual_hum']=0
    data['mapQual_hum']=0
    dataR1Right=data[data["HIV"].str.contains("R1:r")]
    dataR1Right=dataR1Right[~((dataR1Right["R1HUM_TS"]<dataR1Right["R1HIV_TS"])&(dataR1Right["R1HUM_TE"]>dataR1Right["R1HIV_TE"]))]
    dataR1Right=dataR1Right[~((dataR1Right["R1HIV_TS"]<dataR1Right["R1HUM_TS"])&(dataR1Right["R1HIV_TE"]>dataR1Right["R1HUM_TE"]))]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(int).astype(str)+":"+dataR1Right['R1HIV_RS'].astype(int).astype(str)
    dataR1Right['HUM']=dataR1Right['R1HUM_RE'].astype(int)
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1HUM_ID+":"+dataR1Right.R1HIV_ID
    dataR1Right["orient"]="R1-hum:hiv"
    dataR1Right["overlap"]=dataR1Right["overlapR1"]
    dataR1Right["gap"]=dataR1Right["gapR1"]
    dataR1Right["HUM_TS"]=dataR1Right["R1HUM_TS"].astype(int)
    dataR1Right["HUM_TE"]=dataR1Right["R1HUM_TE"].astype(int)
    dataR1Right["HUM_RS"]=dataR1Right["R1HUM_RS"].astype(int)
    dataR1Right["HUM_RE"]=dataR1Right["R1HUM_RE"].astype(int)
    dataR1Right["HIV_TS"]=dataR1Right["R1HIV_TS"].astype(int)
    dataR1Right["HIV_TE"]=dataR1Right["R1HIV_TE"].astype(int)
    dataR1Right["HIV_RS"]=dataR1Right["R1HIV_RS"].astype(int)
    dataR1Right["HIV_RE"]=dataR1Right["R1HIV_RE"].astype(int)
    dataR1Right["HUM_ID"]=dataR1Right["R1HUM_ID"]
    dataR1Right["HIV_ID"]=dataR1Right["R1HIV_ID"]
    dataR1Right["SEQ"]=dataR1Right["R1HUM_SEQ"]
#     dataR1Right["SEQ"]=dataR1Right[dataR1Right['R1HUM_reversedCurr']==0]["R1HUM_SEQ"]
#     dataR1Right["R_SEQ"]=dataR1Right[dataR1Right['R1HUM_reversedCurr']==16]["R1HUM_SEQ"]
    dataR1Right["HIV_MAPQ"]=dataR1Right["R1HIV_MAPQ"].astype(int)
    dataR1Right["HUM_MAPQ"]=dataR1Right["R1HUM_MAPQ"].astype(int)
    dataR1Right["HIV_AL"]=dataR1Right["HIV_TE"]-dataR1Right["HIV_TS"]-dataR1Right["overlap"]
    dataR1Right["HUM_AL"]=dataR1Right["HUM_TE"]-dataR1Right["HUM_TS"]-dataR1Right["overlap"]
    dataR1Right=dataR1Right[(dataR1Right['HIV_AL']>args.minLen)&(dataR1Right['HUM_AL']>args.minLen)]
    if len(dataR1Right>0):
        dataR1Right[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR1Right.merge(dataR1Right.apply(lambda row: processAligns(row['R1HUM_SEQ'],
                                                                                                    row['R1HIV_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1HUM_reversedCurr"]==16,
                                                                                                    row["R1HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                                   'meanQual_hum_y']]
    dataR1Right["READ_LEN"]=dataR1Right["READ_LEN_R1"]
    dataR1Right["R"]="R1"
    dataR1Right.drop(dropList,axis=1,inplace=True)

    # R1 left
    dataR1Left=data[data["HIV"].str.contains("R1:l")]
    dataR1Left=dataR1Left[~((dataR1Left["R1HUM_TS"]<dataR1Left["R1HIV_TS"])&(dataR1Left["R1HUM_TE"]>dataR1Left["R1HIV_TE"]))]
    dataR1Left=dataR1Left[~((dataR1Left["R1HIV_TS"]<dataR1Left["R1HUM_TS"])&(dataR1Left["R1HIV_TE"]>dataR1Left["R1HUM_TE"]))]
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(int).astype(str)+":"+dataR1Left['R1HUM_RS'].astype(int).astype(str)
    dataR1Left['HUM']=dataR1Left['R1HUM_RS'].astype(int)
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1HIV_ID+":"+dataR1Left.R1HUM_ID
    dataR1Left["orient"]="R1-hiv:hum"
    dataR1Left["overlap"]=dataR1Left["overlapR1"]
    dataR1Left["gap"]=dataR1Left["gapR1"]
    dataR1Left["HUM_TS"]=dataR1Left["R1HUM_TS"].astype(int)
    dataR1Left["HUM_TE"]=dataR1Left["R1HUM_TE"].astype(int)
    dataR1Left["HUM_RS"]=dataR1Left["R1HUM_RS"].astype(int)
    dataR1Left["HUM_RE"]=dataR1Left["R1HUM_RE"].astype(int)
    dataR1Left["HIV_TS"]=dataR1Left["R1HIV_TS"].astype(int)
    dataR1Left["HIV_TE"]=dataR1Left["R1HIV_TE"].astype(int)
    dataR1Left["HIV_RS"]=dataR1Left["R1HIV_RS"].astype(int)
    dataR1Left["HIV_RE"]=dataR1Left["R1HIV_RE"].astype(int)
    dataR1Left["HUM_ID"]=dataR1Left["R1HUM_ID"]
    dataR1Left["HIV_ID"]=dataR1Left["R1HIV_ID"]
    dataR1Left["SEQ"]=dataR1Left["R1HUM_SEQ"]
#     dataR1Left["SEQ"]=dataR1Left[dataR1Left['R1HUM_reversedCurr']==0]["R1HUM_SEQ"]
#     dataR1Left["R_SEQ"]=dataR1Left[dataR1Left['R1HUM_reversedCurr']==16]["R1HUM_SEQ"]
    dataR1Left["HIV_MAPQ"]=dataR1Left["R1HIV_MAPQ"].astype(int)
    dataR1Left["HUM_MAPQ"]=dataR1Left["R1HUM_MAPQ"].astype(int)
    dataR1Left["HIV_AL"]=dataR1Left["HIV_TE"]-dataR1Left["HIV_TS"]-dataR1Left["overlap"]
    dataR1Left["HUM_AL"]=dataR1Left["HUM_TE"]-dataR1Left["HUM_TS"]-dataR1Left["overlap"]
    dataR1Left=dataR1Left[(dataR1Left["HIV_AL"]>args.minLen)&(dataR1Left["HUM_AL"]>args.minLen)]
    if len(dataR1Left)>0:
        dataR1Left[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR1Left.merge(dataR1Left.apply(lambda row: processAligns(row['R1HUM_SEQ'],
                                                                                                    row['R1HIV_SEQ'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['HIV_TS']),
                                                                                                    int(row['HIV_TE']-row['overlap']),
                                                                                                    int(row['HUM_TS']+row['overlap']),
                                                                                                    int(row['HUM_TE']),
                                                                                                    int(row['READ_LEN_R1']),
                                                                                                    row["R1HUM_reversedCurr"]==16,
                                                                                                    row["R1HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                   'meanQual_hum_y']]

    dataR1Left["READ_LEN"]=dataR1Left["READ_LEN_R1"]
    dataR1Left["R"]="R1"
    dataR1Left.drop(dropList,axis=1,inplace=True)

    # R2 right
    dataR2Right=data[data["HIV"].str.contains("R2:r")]
    dataR2Right=dataR2Right[~((dataR2Right["R2HUM_TS"]<dataR2Right["R2HIV_TS"])&(dataR2Right["R2HUM_TE"]>dataR2Right["R2HIV_TE"]))]
    dataR2Right=dataR2Right[~((dataR2Right["R2HIV_TS"]<dataR2Right["R2HUM_TS"])&(dataR2Right["R2HIV_TE"]>dataR2Right["R2HUM_TE"]))]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(int).astype(str)+":"+dataR2Right['R2HIV_RS'].astype(int).astype(str)
    dataR2Right['HUM']=dataR2Right['R2HUM_RE'].astype(int)
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2HUM_ID+":"+dataR2Right.R2HIV_ID
    dataR2Right["orient"]="R2-hum:hiv"
    dataR2Right["overlap"]=dataR2Right["overlapR2"]
    dataR2Right["gap"]=dataR2Right["gapR2"]
    dataR2Right["HUM_TS"]=dataR2Right["R2HUM_TS"].astype(int)
    dataR2Right["HUM_TE"]=dataR2Right["R2HUM_TE"].astype(int)
    dataR2Right["HUM_RS"]=dataR2Right["R2HUM_RS"].astype(int)
    dataR2Right["HUM_RE"]=dataR2Right["R2HUM_RE"].astype(int)
    dataR2Right["HIV_TS"]=dataR2Right["R2HIV_TS"].astype(int)
    dataR2Right["HIV_TE"]=dataR2Right["R2HIV_TE"].astype(int)
    dataR2Right["HIV_RS"]=dataR2Right["R2HIV_RS"].astype(int)
    dataR2Right["HIV_RE"]=dataR2Right["R2HIV_RE"].astype(int)
    dataR2Right["HUM_ID"]=dataR2Right["R2HUM_ID"]
    dataR2Right["HIV_ID"]=dataR2Right["R2HIV_ID"]
    dataR2Right["SEQ"]=dataR2Right["R2HUM_SEQ"]
#     dataR2Right["SEQ"]=dataR2Right[dataR2Right['R2HUM_reversedCurr']==0]["R2HUM_SEQ"]
#     dataR2Right["R_SEQ"]=dataR2Right[dataR2Right['R2HUM_reversedCurr']==16]["R2HUM_SEQ"]
    dataR2Right["HIV_MAPQ"]=dataR2Right["R2HIV_MAPQ"].astype(int)
    dataR2Right["HUM_MAPQ"]=dataR2Right["R2HUM_MAPQ"].astype(int)
    dataR2Right["HIV_AL"]=dataR2Right["HIV_TE"]-dataR2Right["HIV_TS"]-dataR2Right["overlap"]
    dataR2Right["HUM_AL"]=dataR2Right["HUM_TE"]-dataR2Right["HUM_TS"]-dataR2Right["overlap"]
    dataR2Right=dataR2Right[(dataR2Right["HIV_AL"]>args.minLen)&(dataR2Right["HUM_AL"]>args.minLen)]
    if len(dataR2Right)>0:
        dataR2Right[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR2Right.merge(dataR2Right.apply(lambda row: processAligns(row['R2HUM_SEQ'],
                                                                                                    row['R2HIV_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2HUM_reversedCurr"]==16,
                                                                                                    row["R2HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                                   'meanQual_hum_y']]
    dataR2Right["READ_LEN"]=dataR2Right["READ_LEN_R2"]
    dataR2Right["R"]="R2"
    dataR2Right.drop(dropList,axis=1,inplace=True)

    # R2 left
    dataR2Left=data[data["HIV"].str.contains("R2:l")]
    dataR2Left=dataR2Left[~((dataR2Left["R2HUM_TS"]<dataR2Left["R2HIV_TS"])&(dataR2Left["R2HUM_TE"]>dataR2Left["R2HIV_TE"]))]
    dataR2Left=dataR2Left[~((dataR2Left["R2HIV_TS"]<dataR2Left["R2HUM_TS"])&(dataR2Left["R2HIV_TE"]>dataR2Left["R2HUM_TE"]))]
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(int).astype(str)+":"+dataR2Left['R2HUM_RS'].astype(int).astype(str)
    dataR2Left['HUM']=dataR2Left['R2HUM_RS'].astype(int)
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2HIV_ID+":"+dataR2Left.R2HUM_ID
    dataR2Left["orient"]="R2-hiv:hum"
    dataR2Left["overlap"]=dataR2Left["overlapR2"]
    dataR2Left["gap"]=dataR2Left["gapR2"]
    dataR2Left["HUM_TS"]=dataR2Left["R2HUM_TS"].astype(int)
    dataR2Left["HUM_TE"]=dataR2Left["R2HUM_TE"].astype(int)
    dataR2Left["HUM_RS"]=dataR2Left["R2HUM_RS"].astype(int)
    dataR2Left["HUM_RE"]=dataR2Left["R2HUM_RE"].astype(int)
    dataR2Left["HIV_TS"]=dataR2Left["R2HIV_TS"].astype(int)
    dataR2Left["HIV_TE"]=dataR2Left["R2HIV_TE"].astype(int)
    dataR2Left["HIV_RS"]=dataR2Left["R2HIV_RS"].astype(int)
    dataR2Left["HIV_RE"]=dataR2Left["R2HIV_RE"].astype(int)
    dataR2Left["HUM_ID"]=dataR2Left["R2HUM_ID"]
    dataR2Left["HIV_ID"]=dataR2Left["R2HIV_ID"]
    dataR2Left["SEQ"]=dataR2Left["R2HUM_SEQ"]
#     dataR2Left["SEQ"]=dataR2Left[dataR2Left['R2HUM_reversedCurr']==0]["R2HUM_SEQ"]
#     dataR2Left["R_SEQ"]=dataR2Left[dataR2Left['R2HUM_reversedCurr']==16]["R2HUM_SEQ"]
    dataR2Left["HIV_MAPQ"]=dataR2Left["R2HIV_MAPQ"].astype(int)
    dataR2Left["HUM_MAPQ"]=dataR2Left["R2HUM_MAPQ"].astype(int)
    dataR2Left["HIV_AL"]=dataR2Left["HIV_TE"]-dataR2Left["HIV_TS"]-dataR2Left["overlap"]
    dataR2Left["HUM_AL"]=dataR2Left["HUM_TE"]-dataR2Left["HUM_TS"]-dataR2Left["overlap"]
    dataR2Left=dataR2Left[(dataR2Left["HIV_AL"]>args.minLen)&(dataR2Left["HUM_AL"]>args.minLen)]
    if len(dataR2Left)>0:
        dataR2Left[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR2Left.merge(dataR2Left.apply(lambda row: processAligns(row['R2HUM_SEQ'],
                                                                                                    row['R2HIV_SEQ'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['HIV_TS']),
                                                                                                    int(row['HIV_TE']-row['overlap']),
                                                                                                    int(row['HUM_TS']+row['overlap']),
                                                                                                    int(row['HUM_TE']),
                                                                                                    int(row['READ_LEN_R2']),
                                                                                                    row["R2HUM_reversedCurr"]==16,
                                                                                                    row["R2HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                   'meanQual_hum_y']]

    
    dataR2Left["READ_LEN"]=dataR2Left["READ_LEN_R2"]
    dataR2Left["R"]="R2"
    dataR2Left.drop(dropList,axis=1,inplace=True)

    frames=[dataR1Right,dataR1Left,dataR2Right,dataR2Left]
    df=pd.concat(frames).reset_index(drop=True)
    if not args.overlap<0 and not args.gap<0:
        df=df[(df["overlap"]<=args.overlap)&(df["gap"]<=args.gap)]
    elif not args.overlap<0:
        df=df[df["overlap"]<=args.overlap]
    elif not args.gap<0:
        df=df[df["gap"]<=args.gap]
    else:
        df=df
    
    # df["HIV_AL"]=(df["READ_LEN"]/2)-((df["HIV_AL"]-df["READ_LEN"]/2).abs())
    # df["HUM_AL"]=(df["READ_LEN"]/2)-((df["HUM_AL"]-df["READ_LEN"]/2).abs())

    k=float(args.minLen)
    ssAl=args.steepSlopeAL
    df["HIV_AL_score"]=(((((df['HIV_AL']-k)/k) \
                        /(((ssAl/k)+((df['HIV_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    df["HUM_AL_score"]=(((((df['HUM_AL']-k)/k) \
                        /(((ssAl/k)+((df['HUM_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of HIV alignment length score

    df['jointEntropy']=((df['entropyScore_hiv'] \
                    +df['entropyScore_hum']) \
                    /(2))
    df['jointAlLen']=((df['HUM_AL_score'] \
                    +df['HIV_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.HUM,df.HIV,df.overlap,df.gap))))
    
    return df

def filterOverlapCombineUnpaired(data,args):
    #right
    data['entropyScore_hiv']=0
    data['meanQual_hiv']=0
    data['mapQual_hiv']=0
    data['entropyScore_hum']=0
    data['meanQual_hum']=0
    data['mapQual_hum']=0
    dataRight=data[data["HIV"].str.contains("r")]
    dataRight=dataRight[~((dataRight["HUM_TS"]<dataRight["HIV_TS"])&(dataRight["HUM_TE"]>dataRight["HIV_TE"]))]
    dataRight=dataRight[~((dataRight["HIV_TS"]<dataRight["HUM_TS"])&(dataRight["HIV_TE"]>dataRight["HUM_TE"]))]
    dataRight["ins"]=dataRight["HIV_TS"]-dataRight["HUM_TE"]
    dataRight["split"]=dataRight['HUM_RE'].astype(int).astype(str)+":"+dataRight['HIV_RS'].astype(int).astype(str)
    dataRight['HUM']=dataRight['HUM_RE'].astype(int)
    dataRight["comb"]=dataRight.split+"@"+dataRight.HUM_ID+":"+dataRight.HIV_ID
    dataRight["orient"]="hum:hiv"
    dataRight["overlap"]=dataRight["overlap"]
    dataRight["gap"]=dataRight["gap"]
    dataRight["HUM_TS"]=dataRight["HUM_TS"].astype(int)
    dataRight["HUM_TE"]=dataRight["HUM_TE"].astype(int)
    dataRight["HUM_RS"]=dataRight["HUM_RS"].astype(int)
    dataRight["HUM_RE"]=dataRight["HUM_RE"].astype(int)
    dataRight["HIV_TS"]=dataRight["HIV_TS"].astype(int)
    dataRight["HIV_TE"]=dataRight["HIV_TE"].astype(int)
    dataRight["HIV_RS"]=dataRight["HIV_RS"].astype(int)
    dataRight["HIV_RE"]=dataRight["HIV_RE"].astype(int)
    dataRight["HIV_MAPQ"]=dataRight["HIV_MAPQ"].astype(int)
    dataRight["HUM_MAPQ"]=dataRight["HUM_MAPQ"].astype(int)
    dataRight["SEQ"]=dataRight["HUM_SEQ"]
#     dataRight["SEQ"]=dataRight[dataRight['HUM_reversedCurr']==0]["HUM_SEQ"]
#     dataRight["R_SEQ"]=dataRight[dataRight['HUM_reversedCurr']==16]["HUM_SEQ"]
    dataRight["HIV_AL"]=dataRight["HIV_TE"]-dataRight["HIV_TS"]-dataRight["overlap"]
    dataRight["HUM_AL"]=dataRight["HUM_TE"]-dataRight["HUM_TS"]-dataRight["overlap"]
    dataRight=dataRight[(dataRight['HIV_AL']>args.minLen)&(dataRight['HUM_AL']>args.minLen)]
    if len(dataRight>0):
        dataRight[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataRight.merge(dataRight.apply(lambda row: processAligns(row['HUM_SEQ'],
                                                                                                    row['HIV_SEQ'],
                                                                                                    row['QUAL'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap']),
                                                                                                    int(row['READ_LEN']),
                                                                                                    row["HUM_reversedCurr"]==16,
                                                                                                    row["HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                                   'meanQual_hum_y']]
    dataRight["READ_LEN"]=dataRight["READ_LEN"]
    dataRight["R"]=0

    # left
    dataLeft=data[data["HIV"].str.contains("l")]
    dataLeft=dataLeft[~((dataLeft["HUM_TS"]<dataLeft["HIV_TS"])&(dataLeft["HUM_TE"]>dataLeft["HIV_TE"]))]
    dataLeft=dataLeft[~((dataLeft["HIV_TS"]<dataLeft["HUM_TS"])&(dataLeft["HIV_TE"]>dataLeft["HUM_TE"]))]
    dataLeft["ins"]=dataLeft["HUM_TS"]-dataLeft["HIV_TE"]
    dataLeft["split"]=dataLeft['HIV_RE'].astype(int).astype(str)+":"+dataLeft['HUM_RS'].astype(int).astype(str)
    dataLeft['HUM']=dataLeft['HUM_RS'].astype(int)
    dataLeft["comb"]=dataLeft.split+"@"+dataLeft.HIV_ID+":"+dataLeft.HUM_ID
    dataLeft["orient"]="hiv:hum"
    dataLeft["overlap"]=dataLeft["overlap"]
    dataLeft["gap"]=dataLeft["gap"]
    dataLeft["HUM_TS"]=dataLeft["HUM_TS"].astype(int)
    dataLeft["HUM_TE"]=dataLeft["HUM_TE"].astype(int)
    dataLeft["HUM_RS"]=dataLeft["HUM_RS"].astype(int)
    dataLeft["HUM_RE"]=dataLeft["HUM_RE"].astype(int)
    dataLeft["HIV_TS"]=dataLeft["HIV_TS"].astype(int)
    dataLeft["HIV_TE"]=dataLeft["HIV_TE"].astype(int)
    dataLeft["HIV_RS"]=dataLeft["HIV_RS"].astype(int)
    dataLeft["HIV_RE"]=dataLeft["HIV_RE"].astype(int)
    dataLeft["HIV_MAPQ"]=dataLeft["HIV_MAPQ"].astype(int)
    dataLeft["HUM_MAPQ"]=dataLeft["HUM_MAPQ"].astype(int)
    dataLeft["SEQ"]=dataLeft["HUM_SEQ"]
#     dataLeft["SEQ"]=dataLeft[dataLeft['HUM_reversedCurr']==0]["HUM_SEQ"]
#     dataLeft["R_SEQ"]=dataLeft[dataLeft['HUM_reversedCurr']==16]["HUM_SEQ"]
    dataLeft["HIV_AL"]=dataLeft["HIV_TE"]-dataLeft["HIV_TS"]-dataLeft["overlap"]
    dataLeft["HUM_AL"]=dataLeft["HUM_TE"]-dataLeft["HUM_TS"]-dataLeft["overlap"]
    dataLeft=dataLeft[(dataLeft["HIV_AL"]>args.minLen)&(dataLeft["HUM_AL"]>args.minLen)]
    if len(dataLeft)>0:
        dataLeft[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataLeft.merge(dataLeft.apply(lambda row: processAligns(row['HUM_SEQ'],
                                                                                            row['HIV_SEQ'],
                                                                                            row['QUAL'],
                                                                                            int(row['HIV_TS']),
                                                                                            int(row['HIV_TE']-row['overlap']),
                                                                                            int(row['HUM_TS']+row['overlap']),
                                                                                            int(row['HUM_TE']),
                                                                                            int(row['READ_LEN']),
                                                                                            row["HUM_reversedCurr"]==16,
                                                                                            row["HIV_reversedCurr"]==16),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                   'meanQual_hum_y']]

    dataLeft["READ_LEN"]=dataLeft["READ_LEN"]
    dataLeft["R"]=0

    frames=[dataRight,dataLeft]
    df=pd.concat(frames).reset_index(drop=True)
    if not args.overlap<0 and not args.gap<0:
        df=df[(df["overlap"]<=args.overlap)&(df["gap"]<=args.gap)]
    elif not args.overlap<0:
        df=df[df["overlap"]<=args.overlap]
    elif not args.gap<0:
        df=df[df["gap"]<=args.gap]
    else:
        df=df

    k=float(args.minLen)
    ssAl=args.steepSlopeAL
    df["HIV_AL_score"]=(((((df['HIV_AL']-k)/k) \
                        /(((ssAl/k)+((df['HIV_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    df["HUM_AL_score"]=(((((df['HUM_AL']-k)/k) \
                        /(((ssAl/k)+((df['HUM_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of HIV alignment length score

    df['jointEntropy']=((df['entropyScore_hiv'] \
                    +df['entropyScore_hum']) \
                    /(2))
    df['jointAlLen']=((df['HUM_AL_score'] \
                    +df['HIV_AL_score']) \
                    /(2))
    df["scorePrelim"]=(df['jointEntropy'] \
                    *df['jointAlLen'])

    df.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    df=df.round({'scorePrelim': 4})
    df=df[df["scorePrelim"]>=args.score]
    df["SEQ"]=list(zip(df.scorePrelim,list(zip(df.SEQ,df.HUM,df.HIV,df.overlap,df.gap))))

    return df

def addSpan(data,dataPos):
    def testR1(row,dataHIV):
        dataHIVR1=dataHIV[(dataHIV['R1HIV_RS']-row['HIV_RS']>-500)&(dataHIV['R1HIV_RS']-row['HIV_RS']<0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataHIVR1=dataHIVR1[(dataHIVR1['R2HUM_RE']-row['HUM_RE']<500)&(dataHIVR1['R2HUM_RE']-row['HUM_RE']>0)]
        if len(dataHIVR1)>0:
            return [set(dataHIVR1["QNAME"]),len(dataHIVR1)] # return set of reads and count of reads
        return [{''},0]
        
    def testR2(row,dataHIV):
        dataHIVR2=dataHIV[(dataHIV['R2HIV_RE']-row['HIV_RE']<500)&(dataHIV['R2HIV_RE']-row['HIV_RE']>0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataHIVR2=dataHIVR2[(dataHIVR2['R1HUM_RS']-row['HUM_RS']>-500)&(dataHIVR2['R1HUM_RS']-row['HUM_RS']<0)]
        if len(dataHIVR2)>0:
            return [set(dataHIVR2["QNAME"]),len(dataHIVR2)] # return set of reads and count of reads
        return [{''},0]

    dataHIVR1=data[data['HIV'].str.contains('sepR2:1')]
    dataHIVR1=dataHIVR1[~(dataHIVR1['R1HIV_ID']==0)] # check that the hiv is on the same side
    dataPosR1=dataPos[dataPos['orient'].str.contains('hiv:hum')]
    if len(dataPosR1)>0:
        dataPosR1[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR1.apply(lambda row: testR1(row,dataHIVR1),axis=1)])

    dataHIVR2=data[data['HIV'].str.contains('sepR1:2')]
    dataHIVR2=dataHIVR2[~(dataHIVR2['R2HIV_ID']==0)] # check that the hiv is on the same side
    dataPosR2=dataPos[dataPos['orient'].str.contains('hum:hiv')]
    if len(dataPosR2)>0:
        dataPosR2[['spanR1-R2',
                    'spanCount']]=pd.DataFrame([x for x in dataPosR2.apply(lambda row: testR2(row,dataHIVR2),axis=1)])

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
        'HUM_RS':{
            'HUM_RS':'min'
        },
        'HUM_RE':{
            'HUM_RE':'max'
        },
        'HIV_RS':{
            'HIV_RS':'min'
        },
        'HIV_RE':{
            'HIV_RE':'max'
        },
        'HIV_AL':{
            'HIV_AL':'sum'
        },
        'HUM_AL':{
            'HUM_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'entropyScore_hiv':{
            'entropyScore_hiv':'sum'
        },
        'entropyScore_hum':{
            'entropyScore_hum':'sum'
        },
        'meanQual_hiv':{
            'meanQual_hiv':'sum'
        },
        'meanQual_hum':{
            'meanQual_hum':'sum'
        },
        'HIV_MAPQ':{
            'HIV_MAPQ':'sum'
        },
        'HUM_MAPQ':{
            'HUM_MAPQ':'sum'
        }
    }

    dataPos=pd.DataFrame(data.groupby(by=["comb","split","HUM_ID","R","orient"])[["QNAME",
                                                                                    "HUM_RS",
                                                                                    "HUM_RE",
                                                                                    "HIV_RS",
                                                                                    "HIV_RE",
                                                                                    "HUM_AL",
                                                                                    "HIV_AL",
                                                                                    "READ_LEN",
                                                                                    "entropyScore_hiv",
                                                                                    "meanQual_hiv",
                                                                                    "entropyScore_hum",
                                                                                    "meanQual_hum",
                                                                                    "HIV_MAPQ",
                                                                                    "HUM_MAPQ",
                                                                                    "SEQ"]].agg(aggregations)).reset_index()
    dataPos.rename(columns={'HUM_ID':'chr'}, inplace=True)            
    return dataPos

def score(dataPos,args,minLen):
    dataPos["READ_LEN"]=dataPos["READ_LEN"].astype(float)/dataPos["count"].astype(float)
    dataPos["HIV_MAPQ"]=dataPos["HIV_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["HUM_MAPQ"]=dataPos["HUM_MAPQ"].astype(float)/dataPos["count"].astype(float)
    dataPos["HIV_AL"]=dataPos["HIV_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["HUM_AL"]=dataPos["HUM_AL"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_hiv"]=dataPos["entropyScore_hiv"].astype(float)/dataPos["count"].astype(float)
    dataPos["entropyScore_hum"]=dataPos["entropyScore_hum"].astype(float)/dataPos["count"].astype(float)
    dataPos["count"]=dataPos["count"].astype(float)

    dataPos["HIV_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["HIV_AL"]-dataPos["READ_LEN"]/2).abs())
    dataPos["HUM_AL"]=(dataPos["READ_LEN"]/2)-((dataPos["HUM_AL"]-dataPos["READ_LEN"]/2).abs())

    k=float(args.minLen)
    ssAl=args.steepSlopeAL
    dataPos["HIV_AL_score"]=(((((dataPos['HIV_AL']-k)/k) \
                        /(((ssAl/k)+((dataPos['HIV_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of human alignment length score

    dataPos["HUM_AL_score"]=(((((dataPos['HUM_AL']-k)/k) \
                        /(((ssAl/k)+((dataPos['HUM_AL']-k)/k)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxAlLenPenalty))) \
                                +args.maxAlLenPenalty # algebraic sigmoid function of HIV alignment length score

    m=float(args.minCount)/2.0
    dataPos["count_score"]=(((((dataPos['count']-m)/m) \
                        /(((1.0/m)+((dataPos['count']-m)/m)**2.0)**0.5)) \
                                /2.0 \
                                +0.5) \
                                /(1.0/(1.0-args.maxCountPenalty))) \
                                +args.maxCountPenalty # algebraic sigmoid function of read count score
    
    # dataPos["HIV_MAPQ_score"]=1-(10**(-(dataPos["HIV_MAPQ"]/10)))
    # dataPos["HUM_MAPQ_score"]=1-(10**(-(dataPos["HUM_MAPQ"]/10)))

    dataPos['jointEntropy']=((dataPos['entropyScore_hiv'] \
                    +dataPos['entropyScore_hum']) \
                    /(2))
    dataPos['jointAlLen']=((dataPos['HUM_AL_score'] \
                    +dataPos['HIV_AL_score']) \
                    /(2))
    dataPos["score"]=(dataPos['jointEntropy'] \
                    *dataPos['count_score'] \
                    *dataPos['jointAlLen'])

    dataPos.drop(["jointEntropy",
                  "jointAlLen"],axis=1,inplace=True)
    dataPos=dataPos.round({'score': 4})
    return dataPos

def approxCloseness(data,args):
    data.sort_values(by=["HUM_RS","hum_nearest_SS"],inplace=True)
    data["diff1"]=abs(data['HUM_RS']-data['HUM_RS'].shift(-1))
    data["diff2"]=abs(data['HUM_RS']-data['HUM_RS'].shift(1))
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
        'entropyScore_hum':{
            'entropyScore_hum':'sum'
        },
        'entropyScore_hiv':{
            'entropyScore_hiv':'sum'
        },
        'HIV_AL':{
            'HIV_AL':'sum'
        },
        'HUM_AL':{
            'HUM_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'HIV_MAPQ':{
            'HIV_MAPQ':'sum'
        },
        'HUM_MAPQ':{
            'HUM_MAPQ':'sum'
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
                                                "entropyScore_hum",
                                                "entropyScore_hiv",
                                                "HIV_AL",
                                                "HUM_AL",
                                                "READ_LEN",
                                                "HIV_MAPQ",
                                                "HUM_MAPQ",
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
        'entropyScore_hum':{
            'entropyScore_hum':'sum'
        },
        'entropyScore_hiv':{
            'entropyScore_hiv':'sum'
        },
        'HIV_AL':{
            'HIV_AL':'sum'
        },
        'HUM_AL':{
            'HUM_AL':'sum'
        },
        'READ_LEN':{
            'READ_LEN':'sum'
        },
        'HIV_MAPQ':{
            'HIV_MAPQ':'sum'
        },
        'HUM_MAPQ':{
            'HUM_MAPQ':'sum'
        }
    }
    
    dfg=pd.DataFrame(data.groupby(by=["hum_nearest_SS",
                                      "chr",
                                      "R",
                                      "uid"])[["comb",
                                                "reads",
                                                "count",
                                                "entropyScore_hum",
                                                "entropyScore_hiv",
                                                "HIV_AL",
                                                "HUM_AL",
                                                "READ_LEN",
                                                "HIV_MAPQ",
                                                "HUM_MAPQ",
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
                      'HUM_RS',
                      'HUM_RE',
                      'hum_nearest_SS']

    finalDF=pd.DataFrame(pd.merge(data,finalBed,on=['chr',
                                                      'HUM_RS',
                                                      'HUM_RE'],how='inner'))
    data=data[~(data['comb'].isin(set(finalDF["comb"])))]
    data['hum_nearest_SS']="-"

    return pd.concat([finalDF,data]).reset_index(drop=True)

#this function should find the closest known acceptor/donor for each of the final records
# we could also implement passing the known donor acceptor lists as parameters
# we could also score, if required by how close the alignment is to the known donor/acceptor

# should be done right before the results cleanup

# also makes sense to keep HIV alignment information
# needs to be kept along with the best sequence information in the tuple

# What needs to be done:
# this has to be done separately based on which side of the sequence the HIV is aligned
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
    dataPos["nearestDonor"]=dataPos.apply(lambda row: findClosest(row["HIV"]),axis=1)
    dataPos["nearestAcceptor"]=dataPos.apply(lambda row: findClosest(row["HIV"]),axis=1)

    return dataPos

def rest(dataPos,args,data,ext,unpaired,baseName,outDir,dirPath):
    if not unpaired:
        dataPos=addSpan(data,dataPos)
        dataPos.loc[dataPos['spanCount'].isnull(),['spanCount']]=dataPos.loc[dataPos['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
        dataPos.loc[dataPos['spanR1-R2'].isnull(),['spanR1-R2']]=dataPos.loc[dataPos['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())

    dataBed=dataPos[['chr','HUM_RS','HUM_RE']].drop_duplicates()
    if len(dataPos)>0:
        dataPos=annotate(dataBed,os.path.abspath(args.annotation),dataPos)
        if not len(dataPos)==0:
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

            dataPos.to_csv(os.path.abspath(args.out)+"_Pos"+ext+".csv",index=False)
            dataPosClean=dataPos[(dataPos['entropyScore_hiv']>args.minEntropy) \
                                &(dataPos['entropyScore_hum']>args.minEntropy) \
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

            dataPosClean[colsOrder].to_csv(os.path.abspath(args.out)+"_Pos"+ext+".clean.csv",index=False)

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
#     dataHIV=pd.read_csv(os.path.abspath(args.input1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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
#     dataHUM=pd.read_csv(os.path.abspath(args.input2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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

#     if (len(dataHIV)==0 or len(dataHUM)==0): #exit if either alignment is empty
#         return
    
#     outDirPOS=outDir+"/Positions/"
#     if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
#         os.mkdir(os.path.abspath(outDir+"/Positions/"))
#     if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
#         os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

#     if len(dataHIV)>0 and len(dataHUM)>0:
#         # first calculate how many parts the dataFrame will have
#         # also, should preprocess dataHIV as much as possible before parallelization
#         # dont forget about unpaired option.
#         extractFlagBits(dataHIV)
#         extractFlagBits(dataHUM)

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
    dataHIV=pd.read_csv(os.path.abspath(args.input1),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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
    dataHUM=pd.read_csv(os.path.abspath(args.input2),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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

    if (len(dataHIV)==0 or len(dataHUM)==0): #exit if either alignment is empty
        return
    
    outDirPOS=outDir+"/Positions/"
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
        os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

    if len(dataHIV)>0 and len(dataHUM)>0:
        # extract flag information
        extractFlagBits(dataHIV)
        extractFlagBits(dataHUM)
        unpaired=False
        if len(dataHIV[dataHIV["paired"]>0])==0 and len(dataHUM[dataHUM["paired"]>0])==0:
            unpaired=True
        
        if args.spliced is not None and not unpaired:
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
                dataSpliced=dataSpliced[dataSpliced['CIGAR'].str.contains("N")]
                extractFlagBits(dataSpliced)
                dataSpliced["tid"]=dataSpliced['QNAME']+dataSpliced['firstRead'].astype(str)+dataSpliced['lastRead'].astype(str)
                dataHUM["tid"]=dataHUM['QNAME']+dataHUM['firstRead'].astype(str)+dataHUM['lastRead'].astype(str)
                dataHIV["tid"]=dataHIV['QNAME']+dataHIV['firstRead'].astype(str)+dataHIV['lastRead'].astype(str)
                dataHIV=dataHIV[~dataHIV['tid'].isin(set(dataSpliced['tid']))]
                dataHUM=dataHUM[dataHUM['tid'].isin(set(dataHIV['tid']))]
                dataHIV=dataHIV[dataHIV['tid'].isin(set(dataHUM['tid']))]
                dataHIV.drop(['tid'],axis=1,inplace=True)
                dataHUM.drop(['tid'],axis=1,inplace=True)
            
            else:
                print("Spliced file does not exist")
        elif unpaired==False: # still remove all human which are not in hiv
            dataHUM["tid"]=dataHUM['QNAME']+dataHUM['firstRead'].astype(str)+dataHUM['lastRead'].astype(str)
            dataHIV["tid"]=dataHIV['QNAME']+dataHIV['firstRead'].astype(str)+dataHIV['lastRead'].astype(str)
            dataHIV=dataHIV[dataHIV['tid'].isin(set(dataHUM['tid']))]
            dataHUM=dataHUM[dataHUM['tid'].isin(set(dataHIV['tid']))]
            dataHIV.drop(['tid'],axis=1,inplace=True)
            dataHUM.drop(['tid'],axis=1,inplace=True)
        else:
            dataHUM["tid"]=dataHUM['QNAME']
            dataHIV["tid"]=dataHIV['QNAME']
            dataHIV=dataHIV[dataHIV['tid'].isin(set(dataHUM['tid']))]
            dataHUM=dataHUM[dataHUM['tid'].isin(set(dataHIV['tid']))]
            dataHIV.drop(['tid'],axis=1,inplace=True)
            dataHUM.drop(['tid'],axis=1,inplace=True)

        # extract start and end for both template and reference
        dataHIV=extractStartEnd(dataHIV)
        dataHUM=extractStartEnd(dataHUM)

        dataHUM,dataHIV=filterReads(dataHUM,dataHIV)
        if len(dataHUM)==0:
            return

        data=pd.DataFrame(dataHIV["QNAME"]).reset_index(drop=True)
        dataHUM=dataHUM.reset_index(drop=True)
        dataHIV=dataHIV.reset_index(drop=True)
        if unpaired:
            data=createDataUnpaired(data,dataHUM,dataHIV)
        else:
            data=createData(data,dataHUM,dataHIV)

        del dataHUM
        del dataHIV
        data.replace('', np.nan,inplace=True)
        data.fillna(0,inplace=True)
        if unpaired:
            data=leftRightUnpaired(data,minLen)
            data=overlapUnpaired(data)
        else:
            data=leftRight(data,minLen)
            data=overlap(data)
        data.drop_duplicates(inplace=True)

        d=pd.DataFrame([])
        if unpaired:
            d=filterOverlapCombineUnpaired(data,args)
        else:
            d=filterOverlapCombine(data,args)
        d=d[(d["entropyScore_hum"]>args.minEntropy)&(d["entropyScore_hiv"]>args.minEntropy)]
        dataPos=d.copy(deep=True)
        dataPos=findSupport(dataPos,minLen,unpaired)
        if len(dataPos)>0:
            rest(dataPos,args,data,"",unpaired,baseName,outDir,dirPath)
    
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

    else:
        print('not real path')
        return

def hiv(argv):

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
                              default=5,
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
                              default=0.6,
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
    parser.add_argument('--spliced',
                              required=False,
                              type=str,
                              help="spliced end-to-end alignment. Reads will be subtracted from the main alignments and will not be reported in the final report of integrations sites")
    parser.set_defaults(func=main)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    hiv(sys.argv[1:])
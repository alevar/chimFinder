import pandas as pd
import numpy as np
import numba
import math
import os
import subprocess
import sys
import glob
import itertools
import argparse
import scipy
import warnings
from pybedtools import BedTool
import matplotlib as mpl
mpl.use('Agg')
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
import seaborn as sbn
from pandas.tools.plotting import scatter_matrix
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
        if (row['reversedCurr']==16)and(reference==False):
            return str(len(row["SEQ"])-1-end)+":"+str(len(row["SEQ"])-1-start)
        else:
            return str(start)+":"+str(end)
    else:
        pass

def alMap(row,dataHIV,dataHUM):
    dataTMP_1=dataHIV[dataHIV['QNAME']==row["QNAME"]]
    dataTMP_2=dataHUM[dataHUM['QNAME']==row["QNAME"]]
    flag=0
    if(dataTMP_1["firstRead"].iloc[0]==64): # hum aligns on the first read
        flag|=2
        if(dataTMP_2["firstRead"].iloc[0]==64): #hum and hiv on the same first read
            flag|=4

            meanTMP_1=(dataTMP_1["Template_start"].iloc[0]+dataTMP_1["Template_end"].iloc[0])/2
            meanTMP_2=(dataTMP_2["Template_start"].iloc[0]+dataTMP_2["Template_end"].iloc[0])/2
            if(meanTMP_1<meanTMP_2): # hiv to the right
                flag|=8

    if(dataTMP_1["lastRead"].iloc[0]==128): #hum aligns on last read
        flag|=16
        if(dataTMP_2["lastRead"].iloc[0]==128): #hum and hiv on the same second read
            flag|=32

            meanTMP_1=(dataTMP_1["Template_start"].iloc[0]+dataTMP_1["Template_end"].iloc[0])/2
            meanTMP_2=(dataTMP_2["Template_start"].iloc[0]+dataTMP_2["Template_end"].iloc[0])/2
            if(meanTMP_1<meanTMP_2): # hiv to the right
                flag|=64

    return flag

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
    data.loc[data['overlapR1']<0,['overlapR1']]=0
    data["overlapR2"]=data[['R2HUM_TE','R2HIV_TE']].min(axis=1)-data[['R2HUM_TS','R2HIV_TS']].max(axis=1)
    data.loc[data['overlapR2']<0,['overlapR2']]=0
    return data

def overlapUnpaired(data):
    data.reset_index(inplace=True,drop=True)
    data["overlap"]=data[['HUM_TE','HIV_TE']].min(axis=1)-data[['HUM_TS','HIV_TS']].max(axis=1) #min end - max start
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

# extract start and end for both template and reference
def extractStartEndOld(data):
    data["CIGAR"].replace("*",np.nan,inplace=True)
    data.dropna(axis=0,inplace=True)
    t=data.values[:,5].astype("S")
    tx=t.view(np.uint8).reshape(-1, t.itemsize)
    ty=data.reversedCurr.values.astype("I")
    data['lenSEQ']=data["SEQ"].str.len()
    tz=data.lenSEQ.values.astype("I")
    tu=data.POS.values.astype("I")

    data["Template_start"]=calcAlignmentTemplateStart(tx,ty,tz)
    data["Template_end"]=calcAlignmentTemplateEnd(tx,ty,tz)
    data["Reference_start"]=calcAlignmentReferenceStart(tx,ty,tz,tu)
    data["Reference_end"]=calcAlignmentReferenceEnd(tx,ty,tz,tu)


def extractStartEnd(data):
    data["CIGAR"].replace("*",np.nan,inplace=True)
    data.dropna(axis=0,inplace=True)

    data["READ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$").replace(np.nan,0).astype(int)
    data["END"]=data.READ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[S]").replace(np.nan,0).astype(int)

    dataHIV6=data[data["reversedCurr"]==16]
    data0=data[data["reversedCurr"]==0]
    dataHIV6["Template_start"]=dataHIV6.READ_LEN-dataHIV6.END
    dataHIV6["Template_end"]=dataHIV6.READ_LEN-dataHIV6.CIGAR_PRE-1
    data0["Template_start"]=data0.CIGAR_PRE
    data0["Template_end"]=data0.END

    dataHIV6["Reference_start"]=dataHIV6.READ_LEN-dataHIV6.END+dataHIV6.POS-dataHIV6.Template_start
    dataHIV6["Reference_end"]=dataHIV6.READ_LEN-dataHIV6.CIGAR_PRE-1+dataHIV6.POS-dataHIV6.Template_start
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE

    data=pd.concat([dataHIV6,data0]).reset_index(drop=True)
    data.drop(["CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
    return data

# filtering the reads based on the flags:
def filterReads(dataHUM,dataHIV):
    #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
    dataHUM=dataHUM[(dataHUM["secondaryAlignment"]==0)&(dataHUM["PCRdup"]==0)&(dataHUM["suppAl"]==0)&(dataHUM["noPassFilter"]==0)]
    dataHIV=dataHIV[(dataHIV["secondaryAlignment"]==0)&(dataHIV["PCRdup"]==0)&(dataHIV["suppAl"]==0)&(dataHIV["noPassFilter"]==0)]
    #now we can exclude paired-end reads which do not contain information from both references
    setHIV=set(dataHIV["QNAME"])
    setHUM=set(dataHUM["QNAME"])
    setQnames=setHIV.intersection(setHUM)
    dataHUM=dataHUM[dataHUM["QNAME"].isin(setQnames)]
    dataHIV=dataHIV[dataHIV["QNAME"].isin(setQnames)]
    if (len(setQnames)==0):
        return pd.DataFrame([]),pd.DataFrame([])
    return dataHUM, dataHIV

def createData(data,dataHUM,dataHIV):
    df=pd.DataFrame([])
    dataHUM=dataHUM[dataHUM['QNAME'].isin(set(data['QNAME']))]
    dataHUMR1=dataHUM[dataHUM['firstRead']==64][["QNAME",
                                                "Template_start",
                                                "Template_end",
                                                "RNAME",
                                                "Reference_start",
                                                "Reference_end",
                                                "READ_LEN",
                                                "SEQ",
                                                "QUAL",
                                                "MAPQ"]]
    dataHUMR2=dataHUM[dataHUM['lastRead']==128][["QNAME",
                                                "Template_start",
                                                "Template_end",
                                                "RNAME",
                                                "Reference_start",
                                                "Reference_end",
                                                "READ_LEN",
                                                "SEQ",
                                                "QUAL",
                                                "MAPQ"]]
    dataHIVR1=dataHIV[dataHIV['firstRead']==64][["QNAME",
                                                "Template_start",
                                                "Template_end",
                                                "RNAME",
                                                "Reference_start",
                                                "Reference_end",
                                                "READ_LEN",
                                                "SEQ",
                                                "QUAL",
                                                "MAPQ"]]
    dataHIVR2=dataHIV[dataHIV['lastRead']==128][["QNAME",
                                                "Template_start",
                                                "Template_end",
                                                "RNAME",
                                                "Reference_start",
                                                "Reference_end",
                                                "READ_LEN",
                                                "SEQ",
                                                "QUAL",
                                                "MAPQ"]]
    df[["QNAME",
        "R1HUM_TS",
        "R1HUM_TE",
        "R1HUM_ID",
        "R1HUM_RS",
        "R1HUM_RE",
        "READ_LEN_R1",
        "SEQ_R1",
        "QUAL_R1",
        "R1HUM_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHUMR1,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "READ_LEN",
                                                                                                        "SEQ",
                                                                                                        "QUAL",
                                                                                                        "MAPQ"]]
    df[["QNAME",
        "R2HUM_TS",
        "R2HUM_TE",
        "R2HUM_ID",
        "R2HUM_RS",
        "R2HUM_RE",
        "READ_LEN_R2",
        "SEQ_R2",
        "QUAL_R2",
        "R2HUM_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHUMR2,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "READ_LEN",
                                                                                                        "SEQ",
                                                                                                        "QUAL",
                                                                                                        "MAPQ"]]
    df[["QNAME",
        "R1HIV_TS",
        "R1HIV_TE",
        "R1HIV_ID",
        "R1HIV_RS",
        "R1HIV_RE",
        "R1HIV_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHIVR1,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "MAPQ"]]
    df[["QNAME",
        "R2HIV_TS",
        "R2HIV_TE",
        "R2HIV_ID",
        "R2HIV_RS",
        "R2HIV_RE",
        "R2HIV_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHIVR2,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "MAPQ"]]
    df[["R1HUM_ID",
        "R2HUM_ID",
        "R1HIV_ID",
        "R2HIV_ID",
        "R1HUM_MAPQ",
        "R2HUM_MAPQ",
        "R1HIV_MAPQ",
        "R2HIV_MAPQ"]].fillna('',inplace=True)
    df.fillna(0,inplace=True)
    return df

def createDataUnpaired(data,dataHUM,dataHIV):
    df=pd.DataFrame([])
    dataHUM=dataHUM[dataHUM['QNAME'].isin(set(data['QNAME']))]
    dataHUM=dataHUM[["QNAME",
                    "Template_start",
                    "Template_end",
                    "RNAME",
                    "Reference_start",
                    "Reference_end",
                    "READ_LEN",
                    "SEQ",
                    "QUAL",
                    "MAPQ"]]
    dataHIV=dataHIV[["QNAME",
                    "Template_start",
                    "Template_end",
                    "RNAME",
                    "Reference_start",
                    "Reference_end",
                    "READ_LEN",
                    "SEQ",
                    "QUAL",
                    "MAPQ"]]
    df[["QNAME",
        "HUM_TS",
        "HUM_TE",
        "HUM_ID",
        "HUM_RS",
        "HUM_RE",
        "READ_LEN",
        "SEQ",
        "QUAL",
        "HUM_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHUM,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "READ_LEN",
                                                                                                        "SEQ",
                                                                                                        "QUAL",
                                                                                                        "MAPQ"]]
    df[["QNAME",
        "HIV_TS",
        "HIV_TE",
        "HIV_ID",
        "HIV_RS",
        "HIV_RE",
        "HIV_MAPQ"]]=pd.DataFrame(pd.merge(data,dataHIV,how='left',on='QNAME')).reset_index(drop=True)[["QNAME",
                                                                                                        "Template_start",
                                                                                                        "Template_end",
                                                                                                        "RNAME",
                                                                                                        "Reference_start",
                                                                                                        "Reference_end",
                                                                                                        "MAPQ"]]
    df[["HUM_ID",
        "HUM_ID",
        "HIV_ID",
        "HIV_ID",
        "HUM_MAPQ",
        "HUM_MAPQ",
        "HIV_MAPQ",
        "HIV_MAPQ"]].fillna('',inplace=True)
    df.fillna(0,inplace=True)
    return df

# replace these two methods from stackOverflow
def e(s):
    prob=[float(s.count(c))/len(s) for c in dict.fromkeys(list(s))]
    e=-sum([p * math.log(p) / math.log(2.0) for p in prob])
    return e

def eI(l):
    p=1.0/l
    return -1.0*l*p*math.log(p)/math.log(2.0)

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

# How to compute the mimimum entropy for a string of a given length and characters

#consider using technique described in https://academic.oup.com/bioinformatics/article/27/8/1061/227307/Topological-entropy-of-DNA-sequences
#another paper: Look at figure 2 https://www.nature.com/articles/srep19788#f2
#more https://www.xycoon.com/normalized_entropy.htm

def processAligns(seq,qual,i1_hiv,i2_hiv,i1_hum,i2_hum):

    entropyScore_hiv=0
    entropyScore_hum=0
    meanQual_hiv=0
    meanQual_hum=0

    s_hiv=seq[i1_hiv:i2_hiv]
    if not len(s_hiv)==0:
        # ideal_hiv=eI(len(s_hiv))
        # if not ideal_hiv==0:
        #     entropyScore_hiv=abs(e(s_hiv))/ideal_hiv
        entropyScore_hiv=topologicalNormalizedEntropy(s_hiv)
        q_hiv=qual[i1_hiv:i2_hiv]
        meanQual_hiv=sum([ord(x)-33 for x in q_hiv])/len(q_hiv)

    s_hum=seq[i1_hum:i2_hum]
    if not len(s_hum)==0:
        # ideal_hum=eI(len(s_hum))
        # if not ideal_hum==0:
        #     entropyScore_hum=abs(e(s_hum))/ideal_hum
        entropyScore_hum=topologicalNormalizedEntropy(s_hum)
        q_hum=qual[i1_hum:i2_hum]
        meanQual_hum=sum([ord(x)-33 for x in q_hum])/len(q_hum)

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
def filterOverlapCombine(data,minLen):
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
              "HIV",
              "SEQ_R1",
              "SEQ_R2",
              "QUAL_R1",
              "QUAL_R2",
              "R1HUM_MAPQ",
              "R1HIV_MAPQ",
              "R2HUM_MAPQ",
              "R2HIV_MAPQ",
              "READ_LEN_R1",
              "READ_LEN_R2"]

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
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1HUM_ID+":"+dataR1Right.R1HIV_ID
    dataR1Right["orient"]="R1-hum:hiv"
    dataR1Right["overlap"]=dataR1Right["overlapR1"]
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
    dataR1Right["HIV_MAPQ"]=dataR1Right["R1HIV_MAPQ"].astype(int)
    dataR1Right["HUM_MAPQ"]=dataR1Right["R1HUM_MAPQ"].astype(int)
    dataR1Right["HIV_AL"]=dataR1Right["HIV_TE"]-dataR1Right["HIV_TS"]-dataR1Right["overlap"]
    dataR1Right["HUM_AL"]=dataR1Right["HUM_TE"]-dataR1Right["HUM_TS"]-dataR1Right["overlap"]
    dataR1Right=dataR1Right[(dataR1Right['HIV_AL']>minLen)&(dataR1Right['HUM_AL']>minLen)]
    if len(dataR1Right>0):
        dataR1Right[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR1Right.merge(dataR1Right.apply(lambda row: processAligns(row['SEQ_R1'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
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
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1HIV_ID+":"+dataR1Left.R1HUM_ID
    dataR1Left["orient"]="R1-hiv:hum"
    dataR1Left["overlap"]=dataR1Left["overlapR1"]
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
    dataR1Left["HIV_MAPQ"]=dataR1Left["R1HIV_MAPQ"].astype(int)
    dataR1Left["HUM_MAPQ"]=dataR1Left["R1HUM_MAPQ"].astype(int)
    dataR1Left["HIV_AL"]=dataR1Left["HIV_TE"]-dataR1Left["HIV_TS"]-dataR1Left["overlap"]
    dataR1Left["HUM_AL"]=dataR1Left["HUM_TE"]-dataR1Left["HUM_TS"]-dataR1Left["overlap"]
    dataR1Left=dataR1Left[(dataR1Left["HIV_AL"]>minLen)&(dataR1Left["HUM_AL"]>minLen)]
    if len(dataR1Left)>0:
        dataR1Left[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR1Left.merge(dataR1Left.apply(lambda row: processAligns(row['SEQ_R1'],
                                                                                                    row['QUAL_R1'],
                                                                                                    int(row['HIV_TS']),
                                                                                                    int(row['HIV_TE']-row['overlap']),
                                                                                                    int(row['HUM_TS']+row['overlap']),
                                                                                                    int(row['HUM_TE'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
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
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2HUM_ID+":"+dataR2Right.R2HIV_ID
    dataR2Right["orient"]="R2-hum:hiv"
    dataR2Right["overlap"]=dataR2Right["overlapR2"]
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
    dataR2Right["HIV_MAPQ"]=dataR2Right["R2HIV_MAPQ"].astype(int)
    dataR2Right["HUM_MAPQ"]=dataR2Right["R2HUM_MAPQ"].astype(int)
    dataR2Right["HIV_AL"]=dataR2Right["HIV_TE"]-dataR2Right["HIV_TS"]-dataR2Right["overlap"]
    dataR2Right["HUM_AL"]=dataR2Right["HUM_TE"]-dataR2Right["HUM_TS"]-dataR2Right["overlap"]
    dataR2Right=dataR2Right[(dataR2Right["HIV_AL"]>minLen)&(dataR2Right["HUM_AL"]>minLen)]
    if len(dataR2Right)>0:
        dataR2Right[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR2Right.merge(dataR2Right.apply(lambda row: processAligns(row['SEQ_R2'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
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
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2HIV_ID+":"+dataR2Left.R2HUM_ID
    dataR2Left["orient"]="R2-hiv:hum"
    dataR2Left["overlap"]=dataR2Left["overlapR2"]
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
    dataR2Left["HIV_MAPQ"]=dataR2Left["R2HIV_MAPQ"].astype(int)
    dataR2Left["HUM_MAPQ"]=dataR2Left["R2HUM_MAPQ"].astype(int)
    dataR2Left["HIV_AL"]=dataR2Left["HIV_TE"]-dataR2Left["HIV_TS"]-dataR2Left["overlap"]
    dataR2Left["HUM_AL"]=dataR2Left["HUM_TE"]-dataR2Left["HUM_TS"]-dataR2Left["overlap"]
    dataR2Left=dataR2Left[(dataR2Left["HIV_AL"]>minLen)&(dataR2Left["HUM_AL"]>minLen)]
    if len(dataR2Left)>0:
        dataR2Left[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataR2Left.merge(dataR2Left.apply(lambda row: processAligns(row['SEQ_R2'],
                                                                                                    row['QUAL_R2'],
                                                                                                    int(row['HIV_TS']),
                                                                                                    int(row['HIV_TE']-row['overlap']),
                                                                                                    int(row['HUM_TS']+row['overlap']),
                                                                                                    int(row['HUM_TE'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                   'meanQual_hum_y']]

    dataR2Left["READ_LEN"]=dataR2Left["READ_LEN_R2"]
    dataR2Left["R"]="R2"
    dataR2Left.drop(dropList,axis=1,inplace=True)

    frames=[dataR1Right,dataR1Left,dataR2Right,dataR2Left]
    df=pd.concat(frames).reset_index(drop=True)
    return df

def filterOverlapCombineUnpaired(data,minLen):
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
    dataRight["comb"]=dataRight.split+"@"+dataRight.HUM_ID+":"+dataRight.HIV_ID
    dataRight["orient"]="hum:hiv"
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
    dataRight["HIV_AL"]=dataRight["HIV_TE"]-dataRight["HIV_TS"]-dataRight["overlap"]
    dataRight["HUM_AL"]=dataRight["HUM_TE"]-dataRight["HUM_TS"]-dataRight["overlap"]
    dataRight=dataRight[(dataRight['HIV_AL']>minLen)&(dataRight['HUM_AL']>minLen)]
    if len(dataRight>0):
        dataRight[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataRight.merge(dataRight.apply(lambda row: processAligns(row['SEQ'],
                                                                                                    row['QUAL'],
                                                                                                    int(row['HIV_TS']+row['overlap']),
                                                                                                    int(row['HIV_TE']),
                                                                                                    int(row['HUM_TS']),
                                                                                                    int(row['HUM_TE']-row['overlap'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
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
    dataLeft["comb"]=dataLeft.split+"@"+dataLeft.HIV_ID+":"+dataLeft.HUM_ID
    dataLeft["orient"]="hiv:hum"
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
    dataLeft["HIV_AL"]=dataLeft["HIV_TE"]-dataLeft["HIV_TS"]-dataLeft["overlap"]
    dataLeft["HUM_AL"]=dataLeft["HUM_TE"]-dataLeft["HUM_TS"]-dataLeft["overlap"]
    dataLeft=dataLeft[(dataLeft["HIV_AL"]>minLen)&(dataLeft["HUM_AL"]>minLen)]
    if len(dataLeft)>0:
        dataLeft[['entropyScore_hiv',
                     'meanQual_hiv',
                     'entropyScore_hum',
                     'meanQual_hum']]=dataLeft.merge(dataLeft.apply(lambda row: processAligns(row['SEQ'],
                                                                                                    row['QUAL'],
                                                                                                    int(row['HIV_TS']),
                                                                                                    int(row['HIV_TE']-row['overlap']),
                                                                                                    int(row['HUM_TS']+row['overlap']),
                                                                                                    int(row['HUM_TE'])),axis=1),left_index=True,right_index=True)[['entropyScore_hiv_y',
                                                                                                                                                                   'meanQual_hiv_y',
                                                                                                                                                                   'entropyScore_hum_y',
                                                                                                                                                                   'meanQual_hum_y']]

    dataLeft["READ_LEN"]=dataLeft["READ_LEN"]
    dataLeft["R"]=0

    frames=[dataRight,dataLeft]
    df=pd.concat(frames).reset_index(drop=True)
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
def findSupport(dataOrig,minLen,args,unpaired):
    data=pd.DataFrame([])
    # if there are multiple minLen in the minLenList then check the the highest does yild more than one
    # if not, then discard and check the next minLen before constructing the dataPos

    if unpaired:
        data=filterOverlapCombineUnpaired(dataOrig,minLen)
    else:
        data=filterOverlapCombine(dataOrig,minLen)
    
    aggregations={
        'QNAME':{
            'count':'count',
            'reads': lambda x: set(x)
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
            'HIV_AL':'mean'
        },
        'HUM_AL':{
            'HUM_AL':'mean'
        },
        'READ_LEN':{
            'READ_LEN':'mean'
        },
        'entropyScore_hiv':{
            'entropyScore_hiv':'mean'
        },
        'entropyScore_hum':{
            'entropyScore_hum':'mean'
        },
        'meanQual_hiv':{
            'meanQual_hiv':'mean'
        },
        'meanQual_hum':{
            'meanQual_hum':'mean'
        },
        'HIV_MAPQ':{
            'HIV_MAPQ':'mean'
        },
        'HUM_MAPQ':{
            'HUM_MAPQ':'mean'
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
                                                                                    "HUM_MAPQ"]].agg(aggregations)).reset_index()
    dataPos.rename(columns={'HUM_ID':'chr'}, inplace=True)
    dataPos["READ_LEN"]=dataPos["READ_LEN"].astype(int)
    dataPos["HIV_AL_score"]=1-((dataPos["HIV_AL"]-dataPos["READ_LEN"]/2).abs()/(dataPos["READ_LEN"]/2))
    dataPos["HUM_AL_score"]=1-((dataPos["HUM_AL"]-dataPos["READ_LEN"]/2).abs()/(dataPos["READ_LEN"]/2))
    dataPos["count_score"]=(dataPos["count"].astype(float)/float(args.minCount))/(((1.0/float(args.minCount))+(dataPos["count"].astype(float)/float(args.minCount))**2)**0.5) # algebraic sigmoid function
    dataPos["HIV_MAPQ_score"]=1-(10**(-(dataPos["HIV_MAPQ"]/10)))
    dataPos["HUM_MAPQ_score"]=1-(10**(-(dataPos["HUM_MAPQ"]/10)))
    # dataPos["score"]=(dataPos['HIV_AL_score'] \
    #                 *dataPos['entropyScore_hiv'] \
    #                 *dataPos['HIV_MAPQ_score'] \
    #                 *dataPos['HUM_AL_score'] \
    #                 *dataPos['entropyScore_hum'] \
    #                 *dataPos['HUM_MAPQ_score']) \
    #                 *dataPos['count_score']
    dataPos["score"]=dataPos['HIV_AL_score'] \
                    *dataPos['entropyScore_hiv'] \
                    *dataPos['HUM_AL_score'] \
                    *dataPos['entropyScore_hum'] \
                    *dataPos['count_score']
    dataPos=dataPos.round({'score': 3})
            
    return dataPos

def approxCloseness(split1,split2):
    res=abs(int(split1.split(":")[0])-int(split2.split(":")[0]))+abs(int(split1.split(":")[2])-int(split2.split(":")[2]))
    try:
        res2=math.log(res)
    except:
        return 0
    return res2

def getSet(row):
    tmpLocalQNAME=list(row["reads"])
    return tmpLocalQNAME

def writeReadNames(outDir,row,fileNameR1,fileNameR2,baseName,dirPath):
    outD=os.path.abspath(outDir)
    tempF=outD+'/tmp/fq'
    stringReads=row['reads'].replace(';','|')

    outPath1_All="'"+outD+"/Positions/"+baseName+"/"+str(row['hum_nearest_SS'].strip('\n')+"@"+row["R"])+"_R1.all.fa'"
    outPath2_All="'"+outD+"/Positions/"+baseName+"/"+str(row['hum_nearest_SS'].strip('\n')+"@"+row["R"])+"_R2.all.fa'"
    cmdR1="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR1+" | seqtk seq -a - > "+outPath1_All
    cmdR2="egrep -A 3 '"+stringReads+"' "+tempF+"/"+fileNameR2+" | seqtk seq -a - > "+outPath2_All
    os.system(cmdR1)
    os.system(cmdR2)

# def writeReadNames(dataPos,dataHIV,outDir,fileName):
#     def helper(row,dataHIV,outDir,fileName):
#         pathR1=
#         pathR2=
#         reads=row['reads'].split(";")
#         dataHIV[(dataHIV['QNAME'].isin(reads))&(dataHIV["R"]=="R1")][["QNAME","SEQ"]].to_csv(outDir+"/")
#     dataPos.apply(lambda row: helper(row,dataHIV,outDir,fileName),axis=1)

def getStats(data,baseName,outDir):
    numSpliceJunctions=0
    with open(outDir+"/hisat/"+baseName+".junctions") as f:
        for i, l in enumerate(f):
            pass
    numSplits=len(data)
    numReads=data["count"].sum()
    numSpliceJunctions=numSpliceJunctions+1
    return pd.DataFrame([[baseName,numSplits,numReads,numSpliceJunctions]],columns=["name",
                                                                                    "numSplits",
                                                                                    "numReads",
                                                                                    "numSpliceHIV"])
    
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

def addSpanOld(row,dataSep):
    dataSpan=pd.DataFrame([])
    humR=row['R']
    hivR='R1'
    if humR=='R1':
        hivR='R2'
    if 'hiv:hum' in row['orient']:
        dataSpan=dataSep[dataSep['HIV'].str.contains('sep'+humR)]  #check that the hum alignment is on the same side
        dataSpan=dataSpan[~(dataSpan[hivR+'HIV_ID']==0)] # check that the hiv is on the same side
        dataSpan=dataSpan[dataSpan[humR+'HUM_ID']==row['chr']] # check that the chromosomes match
        dataSpan=dataSpan[(dataSpan[hivR+'HIV_RS']-row['HIV_RS']>-500)&(dataSpan[hivR+'HIV_RS']-row['HIV_RS']<0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataSpan=dataSpan[(dataSpan[humR+'HUM_RE']-row['HUM_RE']<500)&(dataSpan[humR+'HUM_RE']-row['HUM_RE']>0)] # check that the end of the hum is after the end of the hum in the grouped dataframe
        if len(dataSpan)>0:
            return [set(dataSpan["QNAME"]),len(dataSpan)] # return set of reads and count of reads
    if 'hum:hiv' in row['orient']:
        dataSpan=dataSep[dataSep['HIV'].str.contains('sep'+humR)]  #check that the hum alignment is on the same side
        dataSpan=dataSpan[~(dataSpan[hivR+'HIV_ID']==0)] # check that the hiv is on the same side
        dataSpan=dataSpan[dataSpan[humR+'HUM_ID']==row['chr']] # check that the chromosomes match
        dataSpan=dataSpan[(dataSpan[hivR+'HIV_RE']-row['HIV_RE']<500)&(dataSpan[hivR+'HIV_RE']-row['HIV_RE']>0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataSpan=dataSpan[(dataSpan[humR+'HUM_RS']-row['HUM_RS']>-500)&(dataSpan[humR+'HUM_RS']-row['HUM_RS']<0)] # check that the end of the hum is after the end of the hum in the grouped dataframe
        if len(dataSpan)>0:
            return [set(dataSpan["QNAME"]),len(dataSpan)] # return set of reads and count of reads

    return [{''},0]

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
        'score':{
            'score':'sum'
        }
    }
    
    dfg=pd.DataFrame(data.groupby(by=["hum_nearest_SS",
                                      "chr",
                                      "R"])[["comb",
                                            "reads",
                                            "count",
                                            "spanCount",
                                            "spanR1-R2",
                                            "score"]].agg(aggregations)).reset_index()

    return dfg.sort_values(by='score',ascending=False).reset_index(drop=True)

def groupBySpliceSitesUnpaired(data):
    colNames=["orient"]
        
    data.drop(colNames,axis=1,inplace=True)
    aggregations={
        'reads':{
            'groupsCount':'count',
            'reads': lambda x: ';'.join(set(x))
        },
        'count':{
            'count':'sum'
        },
        'comb':{
            'comb': lambda x: ';'.join(set(x))
        },
        'score':{
            'score':'sum'
        }
    }
    
    dfg=pd.DataFrame(data.groupby(by=["hum_nearest_SS",
                                      "chr",
                                      "R"])[["comb",
                                            "reads",
                                            "count",
                                            "score"]].agg(aggregations)).reset_index()

    return dfg.sort_values(by='score',ascending=False).reset_index(drop=True)

def addAnnotation(row,baseName,outDir):
    firstRead=row['reads'].split(';')[0]
    hum_nearest_SS='-'
    if row['R']=='R1':
        SS5_CMD="grep '"+firstRead+"' "+outDir+"/"+baseName+".txt"+" | awk -F '\t' '{print $4}'"
        hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
        SS3_CMD="grep '"+firstRead+"' "+outDir+"/"+baseName+".txt"+" | awk -F '\t' '{print $5}'"
        hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)
    else:
        SS5_CMD="grep '"+firstRead+"' "+outDir+"/"+baseName+".txt"+" | awk -F '\t' '{print $9}'"
        hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
        SS3_CMD="grep '"+firstRead+"' "+outDir+"/"+baseName+".txt"+" | awk -F '\t' '{print $10}'"
        hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)

    return pd.Series({'hum_nearest_5SS': hum_nearest_5SS, 'hum_nearest_3SS': hum_nearest_3SS})

def annotate(dataBed,annPath,data):
    data.reset_index(drop=True,inplace=True)

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

    # annotation should return - for those which do not intersect a genomic region

# the following function will produce plots for the data
# most importantly it will generate a snapshot plots for genomic regions as follows:
# x-axis is the genomic region such as chromosome
# y-axis is the coverage of the position
def plotSnapshots(data):
    x = np.linspace(0,10,100)
    x[75:] = np.linspace(40,42.5,25)

    y = np.sin(x)

    f,(ax,ax2) = plt.subplots(1,2,sharey=True, facecolor='w')

    # plot the same data on both axes
    ax.plot(x, y)
    ax2.plot(x, y)

    ax.set_xlim(0,7.5)
    ax2.set_xlim(40,42.5)

    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labelright='off')
    ax2.yaxis.tick_right()

    # This looks pretty good, and was fairly painless, but you can get that
    # cut-out diagonal lines look with just a bit more work. The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0,0), (0,1),
    # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.

    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d), (-d,+d), **kwargs)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)

    # What's cool about this is that now if we vary the distance between
    # ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
    # the diagonal lines will move accordingly, and stay right at the tips
    # of the spines they are 'breaking'

    plt.show()

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
        print("Alignment is empty")
        return
    
    outDirPOS=outDir+"/Positions/"
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDir+"/Positions/"+baseName)):
        os.mkdir(os.path.abspath(outDir+"/Positions/"+baseName))

    data=pd.DataFrame([],columns=['comb',
                                    'split',
                                    'reads',
                                    'count',
                                    'orient'])

    if len(dataHIV)>0 and len(dataHUM)>0:
        # extract flag information
        extractFlagBits(dataHIV)
        extractFlagBits(dataHUM)
        unpaired=False
        if len(dataHIV[dataHIV["paired"]>0])==0 and len(dataHUM[dataHUM["paired"]>0])==0:
            unpaired=True
        # extract start and end for both template and reference
        print("1124")
        dataHIV=extractStartEnd(dataHIV)
        dataHUM=extractStartEnd(dataHUM)
        print("1127")
        dataHUM,dataHIV=filterReads(dataHUM,dataHIV)
        if len(dataHUM)==0:
            return

        data=pd.DataFrame(dataHIV["QNAME"]).reset_index(drop=True)
        dataHUM=dataHUM.reset_index(drop=True)
        dataHIV=dataHIV.reset_index(drop=True)
        print("1135")
        if unpaired:
            data=createDataUnpaired(data,dataHUM,dataHIV)
        else:
            data=createData(data,dataHUM,dataHIV)

        data.replace('', np.nan,inplace=True)
        data.fillna(0,inplace=True)
        print("1140")
        if unpaired:
            data=leftRightUnpaired(data,minLen)
            print("1142")
            data=overlapUnpaired(data)
        else:
            data=leftRight(data,minLen)
            print("1142")
            data=overlap(data)
        # data.to_csv(os.path.abspath(args.out)+".chim"+end+".csv",index=False)
        # drop duplicated reads - preserve first occurence
        data.drop_duplicates(inplace=True)
        print('1147')
        dataPos=findSupport(data,minLen,args,unpaired)

        if len(dataPos)>0:
            print('1151')
            if not unpaired:
                dataPos=addSpan(data,dataPos)
                dataPos.loc[dataPos['spanCount'].isnull(),['spanCount']]=dataPos.loc[dataPos['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
                dataPos.loc[dataPos['spanR1-R2'].isnull(),['spanR1-R2']]=dataPos.loc[dataPos['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())

            dataBed=dataPos[['chr','HUM_RS','HUM_RE']].drop_duplicates()
            # dataBed.to_csv(os.path.abspath(args.out)+".bed",sep='\t',header=False,index=False)
            # annotate(os.path.abspath(args.out)+".bed",os.path.abspath(args.annotation),dataPos)
            print('1159')
            dataPos=annotate(dataBed,os.path.abspath(args.annotation),dataPos)
            # bedCMD='bedtools intersect -a '+outDir+'/beds/'+baseName+'.bed'+' -b '+args.annotation+' -wo > '+outDir+'/beds/'+baseName+'.bed.out'
            # os.system(bedCMD)
            print("1163")
            dataPos=dataPos.drop(["split"],axis=1).sort_values(by='score',ascending=False).reset_index(drop=True)
            dataPos["reads"]=dataPos["reads"].str.join(";")
            if not unpaired:
                dataPos["spanR1-R2"]=dataPos["spanR1-R2"].str.join(";")
                dataPos["spanR1-R2"].replace("","-",inplace=True)

            #Now add the nearest human splice site for each chimeric position
            # dataPos[['hum_nearest_5SS',
            #         'hum_nearest_3SS']]=dataPos.apply(lambda row: addAnnotation(row,baseName,outDir),axis=1)
            # dataPos.to_csv(os.path.abspath(args.out)+"_Test"+end+".csv",index=False)

            dataPosMultiAlignment=dataPos[(dataPos['HIV_MAPQ']<args.minQual)|(dataPos['HUM_MAPQ']<args.minQual)]
            dataPosMultiAlignment.to_csv(os.path.abspath(args.out)+"_Pos_lowMAPQ"+end+".csv",index=False)
            dataPos=dataPos[(dataPos['HIV_MAPQ']>=args.minQual)&(dataPos['HUM_MAPQ']>=args.minQual)]

            dataPos.to_csv(os.path.abspath(args.out)+"_Pos"+end+".csv",index=False)

            dataPos=dataPos[~(dataPos['hum_nearest_SS']=="-") \
                            &(dataPos['entropyScore_hiv']>args.minEntropy) \
                            &(dataPos['entropyScore_hum']>args.minEntropy) \
                            &(dataPos['HIV_MAPQ']>args.minQual) \
                            &(dataPos['HUM_MAPQ']>args.minQual) \
                            &(dataPos['score']>0.3)]
            dataPos.to_csv(os.path.abspath(args.out)+"_Pos"+end+".clean.csv",index=False)

            if unpaired:
                data=groupBySpliceSitesUnpaired(dataPos)
                data.to_csv(os.path.abspath(args.out)+"_Group"+end+".csv",index=False)
                data=data[~(data['hum_nearest_SS']=='-')].reset_index(drop=True)
                data.to_csv(os.path.abspath(args.out)+"_Group"+end+".clean.csv")
            else:
                data=groupBySpliceSites(dataPos)
                data.to_csv(os.path.abspath(args.out)+"_Group"+end+".csv",index=False)
                data=data[~(data['hum_nearest_SS']=='-')].reset_index(drop=True)
                data.to_csv(os.path.abspath(args.out)+"_Group"+end+".clean.csv")

            #completely forgot that we can write reads from the sam file
            #should make it much faster
            #for each comb in the Pos.csv file
            #find corresponding positions in the original human file
            #save to csv with \n delimeter

            # childPIDS=[]
            # fileNameR1=baseName+"_R1.fastq"
            # fileNameR2=baseName+"_R2.fastq"
            # data.apply(lambda row: writeReadNames(os.path.abspath(outDir),row,fileNameR1,fileNameR2,baseName,dirPath),axis=1)
            # os.remove(os.path.abspath(outDir)+"/tempF/"+baseName+"_R1"+baseEnd+".fastq")
            # os.remove(os.path.abspath(outDir)+"/tempF/"+baseName+"_R2"+baseEnd+".fastq")
    
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

# another point to implement is to implement duplicate removal.
# Simple duplicates removal (should be sufficient to reduce the number of reads in the df can be as follows):
# 1. in the sam dataframe if multiple reads which are on the same R1/R2 have identical SEQs then remove all but one with highest quality
# 2. Implement quality calculation for both human and hiv alignments

#Think carefully about this:
# perhaps makes more sense to just make step3 as extraction step. Give original fastq file and extract fasta from there
# Reasons:
# 1. Sam does not preserve barcode information in the read name
# 2. R1/R2 pair is not guaranteed to be present

# Scoring. COnsider the following:
# 1. Hum ALignment Length
# 2. HIV alignment length
# 3. Hum alignment quality score
# 4. HIV alignment quality score
# 5. number of reads in groups
# 6. number of reads in splice groups
# 7. entropy of the read complexity

# for scoring the following needs to be done:
# all individual scores need to be normalized to [0,1]
# for counts the following should be done
# the minimum number is defined
# the count is divided (normazlized) by the minimum
# such results in deflation of the overal score when the count is less than minimum
# and inflation of the count is greater than minimum
# in order to achive this the product of all scores can be 

# also need to run old annotation into a separate column to see if different from the new implementation

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# for each annotated region annotation seems to also return a duplicate annotated as -
# needs fixing
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# apply soft thresholding to count similar to the signum function and not a sigmoid

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# try to make extractStartEnd more efficient
# perhaps look at the regular expressions used by Ella in her scripts
# and implement a regular-expression based extraction
# could be faster
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# calculate additional score for Pos without counts
# recalculate when grouping by human sites
# when grouping by sites - so not filter by score prior to grouping. Only after

# figure out how to deal with mapq of 0
# one approach is to produce separate reports for likely uniquely mapped reads - if any chimeras are identified
# then we can show the location on the chromosome
# and build a report for no unique mappings which will tell which reads are chimeric, but without clear human mapping

# Need to separate groupings also by confidence of mapping
# for instance, split data by MAPQ
# for those where MAPQ == 0, we can not make a judgement about the nearest human splice site
# however they may still be very good reads, hence in that section do include MAPQ in the score calculation

#==DONE
# investigate leftRight:
# currently R1 has R1:r/l as opposed to R2:right/left
# need to check which one is being checked later in the software and change to uniform representation

# Again, try to do approximation of closeness of different chimeras by position
# grouping by splicesite does sometimes group rather distant positions
# if a specific closeness paramter is implemented such results could be better trusted

# If use the reference from the ucsc, then must also use the corresponding annotation
# Need to doenload

# Check the reset_index() for the dfs

#===========================================
#=================NOTES=====================
#===========================================
# 1. have a logger ready for the following tasks:
#       - log which minLen values were discarded
#       - log general stats
# 2. #it seems the add.sh is pretty slow as well
#       - perhaps should timeit the execution of the custom python method against the add.sh
#       - would be best not to apply any external methods or as little as possible
#once the add.sh is replaced by the actual function - do one run with both included in order to
#compare the results and make sure the new method is working properly
#verification should be done as follows:
#1. check if anything that is identified in the add.sh column is also present in the same way in the custom column
#2. if any differences are observed - check why - likely a mistake
#3. then check for any which are identified by the custom method but not by the add.sh

#splicing information is currently missing
#this is due to the fact that we can not build a reliable consensus sequence from alignment to the multi-fasta reference
#what we need to do is write the following script without including in the main pipeline
#1. align to a regular selected HIV reference
#2. run hisat2 on that
#3. perform step2 of the analysis pipeline on that data
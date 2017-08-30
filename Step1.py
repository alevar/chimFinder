import pandas as pd
import numpy as np
import numba
import os
import subprocess
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

@numba.guvectorize([(numba.int32[:], numba.uint32[:],numba.uint32[:],numba.int32[:])], '(n),(),()->()')
def calcAlignmentTemplateStart(cigar,reversedCurr,lenSeq,res):
    start=0
    end=0
    cigarListRAW=[]
    for i in cigar:
        cigarListRAW.append(chr(i))
    cigarTxt=''.join(cigarListRAW)
    it=itertools.groupby(cigarTxt, key=str.isdigit)
    cigarList=[]
    for _,x in it:
        tmp="".join(x)
        cigarList.append(tmp.rstrip('\x00'))
    if 'M' in cigarList:
        idxfM=cigarList.index("M") #index of the first M
        idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
        start=start+sum(list(map(int,cigarList[:idxfM-1][::2])))
        end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
        if (reversedCurr==16):
            start=lenSeq-1-end
    res[0]=start

@numba.guvectorize([(numba.int32[:], numba.uint32[:],numba.uint32[:],numba.int32[:])], '(n),(),()->()')
def calcAlignmentTemplateEnd(cigar,reversedCurr,lenSeq,res):
    start=0
    end=0
    cigarListRAW=[]
    for i in cigar:
        cigarListRAW.append(chr(i))
    cigarTxt=''.join(cigarListRAW)
    it=itertools.groupby(cigarTxt, key=str.isdigit)
    cigarList=[]
    for _,x in it:
        tmp="".join(x)
        cigarList.append(tmp.rstrip('\x00'))
    if 'M' in cigarList:
        idxfM=cigarList.index("M") #index of the first M
        idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
        start=start+sum(list(map(int,cigarList[:idxfM-1][::2])))
        end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
        if (reversedCurr==16):
            end=lenSeq-1-start
    res[0]=end

@numba.guvectorize([(numba.int32[:], numba.uint32[:],numba.uint32[:],numba.uint32[:],numba.int32[:])], '(n),(),(),()->()')
def calcAlignmentReferenceStart(cigar,reversedCurr,lenSeq,pos,res):
    start=pos
    end=0
    cigarListRAW=[]
    for i in cigar:
        cigarListRAW.append(chr(i))
    cigarTxt=''.join(cigarListRAW)
    it=itertools.groupby(cigarTxt, key=str.isdigit)
    cigarList=[]
    for _,x in it:
        tmp="".join(x)
        cigarList.append(tmp.rstrip('\x00'))
    if 'M' in cigarList:
        idxfM=cigarList.index("M") #index of the first M
        idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
        end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
    res[0]=start

@numba.guvectorize([(numba.int32[:], numba.uint32[:],numba.uint32[:],numba.uint32[:],numba.int32[:])], '(n),(),(),()->()')
def calcAlignmentReferenceEnd(cigar,reversedCurr,lenSeq,pos,res):
    start=pos
    end=0
    cigarListRAW=[]
    for i in cigar:
        cigarListRAW.append(chr(i))
    cigarTxt=''.join(cigarListRAW)
    it=itertools.groupby(cigarTxt, key=str.isdigit)
    cigarList=[]
    for _,x in it:
        tmp="".join(x)
        cigarList.append(tmp.rstrip('\x00'))
    if 'M' in cigarList:
        idxfM=cigarList.index("M") #index of the first M
        idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
        end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
    res[0]=end

def alMap(row,data1,data2):
    dataTMP_1=data1[data1['QNAME']==row["QNAME"]]
    dataTMP_2=data2[data2['QNAME']==row["QNAME"]]
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
def leftRight(row):
    t=[]
    if not row["R1HUM_ID"]==0:
        if not row["R1HIV_ID"]==0:
            if row["R1HIV_TS"]>row["R1HUM_TS"]:
                t.append("R1:right")
            else:
                t.append("R1:left")
        else:
            t.append("sepR1:2")


    if not row["R2HUM_ID"]==0:
        if not row["R2HIV_ID"]==0:
            if row['R2HIV_TS']>row["R2HUM_TS"]:
                t.append("R2:right")
            else:
                t.append("R2:left")
        else:
            t.append("sepR2:1")

    return "-".join(t)

def overlapR1(row):
    return len(set(range(int(row["R1HUM_TS"]),int(row["R1HUM_TE"]))).intersection(set(range(int(row["R1HIV_TS"]),int(row["R1HIV_TE"])))))

def overlapR2(row):
    return len(set(range(int(row["R2HUM_TS"]),int(row["R2HUM_TE"]))).intersection(set(range(int(row["R2HIV_TS"]),int(row["R2HIV_TE"])))))

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

    data["SEQ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$").replace(np.nan,0).astype(int)
    data["END"]=data.SEQ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[S]").replace(np.nan,0).astype(int)

    data16=data[data["reversedCurr"]==16]
    data0=data[data["reversedCurr"]==0]
    data16["Template_start"]=data16.SEQ_LEN-data16.END
    data16["Template_end"]=data16.SEQ_LEN-data16.CIGAR_PRE-1
    data0["Template_start"]=data0.CIGAR_PRE
    data0["Template_end"]=data0.END

    data16["Reference_start"]=data16.SEQ_LEN-data16.END+data16.POS-data16.Template_start
    data16["Reference_end"]=data16.SEQ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE

    data=pd.concat([data16,data0]).reset_index(drop=True)
    data.drop(["SEQ_LEN","CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
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
    dataHUMR1=dataHUM[dataHUM['firstRead']==64][["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    dataHUMR2=dataHUM[dataHUM['lastRead']==128][["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    dataHIVR1=dataHIV[dataHIV['firstRead']==64][["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    dataHIVR2=dataHIV[dataHIV['lastRead']==128][["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    df[["QNAME","R1HUM_TS","R1HUM_TE","R1HUM_ID","R1HUM_RS","R1HUM_RE"]]=pd.DataFrame(pd.merge(data,dataHUMR1,how='left',on='QNAME')).reset_index(drop=True)[["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    df[["QNAME","R2HUM_TS","R2HUM_TE","R2HUM_ID","R2HUM_RS","R2HUM_RE"]]=pd.DataFrame(pd.merge(data,dataHUMR2,how='left',on='QNAME')).reset_index(drop=True)[["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    df[["QNAME","R1HIV_TS","R1HIV_TE","R1HIV_ID","R1HIV_RS","R1HIV_RE"]]=pd.DataFrame(pd.merge(data,dataHIVR1,how='left',on='QNAME')).reset_index(drop=True)[["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    df[["QNAME","R2HIV_TS","R2HIV_TE","R2HIV_ID","R2HIV_RS","R2HIV_RE"]]=pd.DataFrame(pd.merge(data,dataHIVR2,how='left',on='QNAME')).reset_index(drop=True)[["QNAME","Template_start","Template_end","RNAME","Reference_start","Reference_end"]]
    df[["R1HUM_ID","R2HUM_ID","R1HIV_ID","R2HIV_ID"]].fillna('',inplace=True)
    df.fillna(0,inplace=True)
    return df

def filterOverlapCombine(data,minLen):
    # this function filters by overlap and flanking alignment length
    # further it removes unnecessary data and combines into a single full dataframe with unified naming avoiding R1/R2 conventions
    # output can be saved as .full.csv and then grouped by split position all at once

    # this function allows identifying best reads
    # However, since greater minLen values for HIV and HUM alignments will yield fewer reads but at higher confidence
    # support reads could ignore the min len requirement as long as the split position is identical
    # or perhaps the minLen requirement for the support reads should be lower

    # so the pipeline here should be as follows:
    # 1. Run this function with lower minLenHUM and minLenHIV parameters in order to identfy all possible support reads
    # 2. Run this function with greater minLenHUM and minLenHIV parameters on the dataframe generated in step 1
    #    in order to identify split positions
    # 3. For each split position from df with higher minLen threshold - find all support reads from both dfs
    # 4. Perhaps find high confidence and low confidence reads separately

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
              "HIV"]

    # R1 right
    dataR1Right=data[data["HIV"].str.contains("R1:right")]
    dataR1Right=dataR1Right[~((dataR1Right["R1HUM_TS"]<dataR1Right["R1HIV_TS"])&(dataR1Right["R1HUM_TE"]>dataR1Right["R1HIV_TE"]))]
    dataR1Right=dataR1Right[~((dataR1Right["R1HIV_TS"]<dataR1Right["R1HUM_TS"])&(dataR1Right["R1HIV_TE"]>dataR1Right["R1HUM_TE"]))]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(str)+":"+dataR1Right['R1HIV_RS'].astype(str)
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1HUM_ID+":"+dataR1Right.R1HIV_ID
    dataR1Right["orient"]="R1-hum:hiv"
    dataR1Right["overlap"]=dataR1Right["overlapR1"]
    dataR1Right["HUM_TS"]=dataR1Right["R1HUM_TS"]
    dataR1Right["HUM_TE"]=dataR1Right["R1HUM_TE"]
    dataR1Right["HUM_RS"]=dataR1Right["R1HUM_RS"]
    dataR1Right["HUM_RE"]=dataR1Right["R1HUM_RE"]
    dataR1Right["HIV_TS"]=dataR1Right["R1HIV_TS"]
    dataR1Right["HIV_TE"]=dataR1Right["R1HIV_TE"]
    dataR1Right["HIV_RS"]=dataR1Right["R1HIV_RS"]
    dataR1Right["HIV_RE"]=dataR1Right["R1HIV_RE"]
    dataR1Right["HUM_ID"]=dataR1Right["R1HUM_ID"]
    dataR1Right["HIV_ID"]=dataR1Right["R1HIV_ID"]
    dataR1Right=dataR1Right[(dataR1Right["HIV_TE"]-dataR1Right["HIV_TS"]-dataR1Right["overlap"])>minLen]
    dataR1Right=dataR1Right[(dataR1Right["HUM_TE"]-dataR1Right["HUM_TS"]-dataR1Right["overlap"])>minLen]
    dataR1Right["R"]="R1"
    dataR1Right.drop(dropList,axis=1,inplace=True)

    # R1 left
    dataR1Left=data[data["HIV"].str.contains("R1:left")]
    dataR1Left=dataR1Left[~((dataR1Left["R1HUM_TS"]<dataR1Left["R1HIV_TS"])&(dataR1Left["R1HUM_TE"]>dataR1Left["R1HIV_TE"]))]
    dataR1Left=dataR1Left[~((dataR1Left["R1HIV_TS"]<dataR1Left["R1HUM_TS"])&(dataR1Left["R1HIV_TE"]>dataR1Left["R1HUM_TE"]))]
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(str)+":"+dataR1Left['R1HUM_RS'].astype(str)
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1HIV_ID+":"+dataR1Left.R1HUM_ID
    dataR1Left["orient"]="R1-hiv:hum"
    dataR1Left["overlap"]=dataR1Left["overlapR1"]
    dataR1Left["HUM_TS"]=dataR1Left["R1HUM_TS"]
    dataR1Left["HUM_TE"]=dataR1Left["R1HUM_TE"]
    dataR1Left["HUM_RS"]=dataR1Left["R1HUM_RS"]
    dataR1Left["HUM_RE"]=dataR1Left["R1HUM_RE"]
    dataR1Left["HIV_TS"]=dataR1Left["R1HIV_TS"]
    dataR1Left["HIV_TE"]=dataR1Left["R1HIV_TE"]
    dataR1Left["HIV_RS"]=dataR1Left["R1HIV_RS"]
    dataR1Left["HIV_RE"]=dataR1Left["R1HIV_RE"]
    dataR1Left["HUM_ID"]=dataR1Left["R1HUM_ID"]
    dataR1Left["HIV_ID"]=dataR1Left["R1HIV_ID"]
    dataR1Left=dataR1Left[(dataR1Left["HIV_TE"]-dataR1Left["HIV_TS"]-dataR1Left["overlap"])>minLen]
    dataR1Left=dataR1Left[(dataR1Left["HUM_TE"]-dataR1Left["HUM_TS"]-dataR1Left["overlap"])>minLen]
    dataR1Left["R"]="R1"
    dataR1Left.drop(dropList,axis=1,inplace=True)

    # R2 right
    dataR2Right=data[data["HIV"].str.contains("R2:right")]
    dataR2Right=dataR2Right[~((dataR2Right["R2HUM_TS"]<dataR2Right["R2HIV_TS"])&(dataR2Right["R2HUM_TE"]>dataR2Right["R2HIV_TE"]))]
    dataR2Right=dataR2Right[~((dataR2Right["R2HIV_TS"]<dataR2Right["R2HUM_TS"])&(dataR2Right["R2HIV_TE"]>dataR2Right["R2HUM_TE"]))]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(str)+":"+dataR2Right['R2HIV_RS'].astype(str)
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2HUM_ID+":"+dataR2Right.R2HIV_ID
    dataR2Right["orient"]="R2-hum:hiv"
    dataR2Right["overlap"]=dataR2Right["overlapR2"]
    dataR2Right["HUM_TS"]=dataR2Right["R2HUM_TS"]
    dataR2Right["HUM_TE"]=dataR2Right["R2HUM_TE"]
    dataR2Right["HUM_RS"]=dataR2Right["R2HUM_RS"]
    dataR2Right["HUM_RE"]=dataR2Right["R2HUM_RE"]
    dataR2Right["HIV_TS"]=dataR2Right["R2HIV_TS"]
    dataR2Right["HIV_TE"]=dataR2Right["R2HIV_TE"]
    dataR2Right["HIV_RS"]=dataR2Right["R2HIV_RS"]
    dataR2Right["HIV_RE"]=dataR2Right["R2HIV_RE"]
    dataR2Right["HUM_ID"]=dataR2Right["R2HUM_ID"]
    dataR2Right["HIV_ID"]=dataR2Right["R2HIV_ID"]
    dataR2Right=dataR2Right[(dataR2Right["HIV_TE"]-dataR2Right["HIV_TS"]-dataR2Right["overlap"])>minLen]
    dataR2Right=dataR2Right[(dataR2Right["HUM_TE"]-dataR2Right["HUM_TS"]-dataR2Right["overlap"])>minLen]
    dataR2Right["R"]="R2"
    dataR2Right.drop(dropList,axis=1,inplace=True)

    # R2 left
    dataR2Left=data[data["HIV"].str.contains("R2:left")]
    dataR2Left=dataR2Left[~((dataR2Left["R2HUM_TS"]<dataR2Left["R2HIV_TS"])&(dataR2Left["R2HUM_TE"]>dataR2Left["R2HIV_TE"]))]
    dataR2Left=dataR2Left[~((dataR2Left["R2HIV_TS"]<dataR2Left["R2HUM_TS"])&(dataR2Left["R2HIV_TE"]>dataR2Left["R2HUM_TE"]))]
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(str)+":"+dataR2Left['R2HUM_RS'].astype(str)
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2HIV_ID+":"+dataR2Left.R2HUM_ID
    dataR2Left["orient"]="R2-hiv:hum"
    dataR2Left["overlap"]=dataR2Left["overlapR2"]
    dataR2Left["HUM_TS"]=dataR2Left["R2HUM_TS"]
    dataR2Left["HUM_TE"]=dataR2Left["R2HUM_TE"]
    dataR2Left["HUM_RS"]=dataR2Left["R2HUM_RS"]
    dataR2Left["HUM_RE"]=dataR2Left["R2HUM_RE"]
    dataR2Left["HIV_TS"]=dataR2Left["R2HIV_TS"]
    dataR2Left["HIV_TE"]=dataR2Left["R2HIV_TE"]
    dataR2Left["HIV_RS"]=dataR2Left["R2HIV_RS"]
    dataR2Left["HIV_RE"]=dataR2Left["R2HIV_RE"]
    dataR2Left["HUM_ID"]=dataR2Left["R2HUM_ID"]
    dataR2Left["HIV_ID"]=dataR2Left["R2HIV_ID"]
    dataR2Left=dataR2Left[(dataR2Left["HIV_TE"]-dataR2Left["HIV_TS"]-dataR2Left["overlap"])>minLen]
    dataR2Left=dataR2Left[(dataR2Left["HUM_TE"]-dataR2Left["HUM_TS"]-dataR2Left["overlap"])>minLen]
    dataR2Left["R"]="R2"
    dataR2Left.drop(dropList,axis=1,inplace=True)

    frames=[dataR1Right,dataR1Left,dataR2Right,dataR2Left]
    df=pd.concat(frames).reset_index(drop=True)
    # data["overlap"]=data.apply(lambda row: 0-row["ins"] if row["overlap"]==0 else row["overlap"],axis=1)
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
        dataPosR1[['spanR1-R2','spanCount']]=pd.DataFrame([x for x in dataPosR1.apply(lambda row: testR1(row,dataHIVR1),axis=1)])

    dataHIVR2=data[data['HIV'].str.contains('sepR1:2')]
    dataHIVR2=dataHIVR2[~(dataHIVR2['R2HIV_ID']==0)] # check that the hiv is on the same side
    dataPosR2=dataPos[dataPos['orient'].str.contains('hum:hiv')]
    if len(dataPosR2)>0:
        dataPosR2[['spanR1-R2','spanCount']]=pd.DataFrame([x for x in dataPosR2.apply(lambda row: testR2(row,dataHIVR2),axis=1)])

    frames=[dataPosR1,dataPosR2]
    df=pd.concat(frames).reset_index(drop=True)
    return df

def findSupport(dataOrig, minLenList):
    minLenList_tmp=list(minLenList)
    # this function should do the following:
    # 1. group data by split position
    # 2. add a column of high confidence support reads from the first parameter DF
    # 3. add a column of low confidence support reads from the second parameter DF

    # First group by split and extract high confidence support reads
    minLenList_tmp.sort(reverse=True)
    data=pd.DataFrame([])
    # if there are multiple minLen in the minLenList then check the the highest does yild more than one
    # if not, then discard and check the next minLen before constructing the dataPos

    while len(data)==0 and minLenList_tmp:
        minLen=minLenList_tmp.pop(0)
        data=filterOverlapCombine(dataOrig,minLen)
    dataPos=pd.DataFrame([])
    
    aggregations={
        'QNAME':{
            'count_'+str(minLen):'count',
            'reads_'+str(minLen): lambda x: set(x)
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
        }
    }

    dataPos=pd.DataFrame(data.groupby(by=["comb","split","HUM_ID","R","orient"])[["QNAME","HUM_RS","HUM_RE","HIV_RS","HIV_RE"]].agg(aggregations)).reset_index()
    dataPos.rename(columns={'HUM_ID':'chr'}, inplace=True)
    if not len(dataPos)==0:
        dataPos['allReads']=dataPos['reads_'+str(minLen)]
        
        pastMinLen=minLen
        pastMinLens=[minLen]
        
        for minLen in minLenList_tmp:
            aggregations={
                'QNAME':{
                    'count_'+str(minLen):'count',
                    'reads_'+str(minLen): lambda x: set(x)
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
                }
            }
            
            
            data=filterOverlapCombine(dataOrig,minLen)
            dataPosT=pd.DataFrame([])
            dataPosT=pd.DataFrame(data.groupby(by=["comb","split","HUM_ID","R","orient"])[["QNAME","HUM_RS","HUM_RE","HIV_RS","HIV_RE"]].agg(aggregations)).reset_index()
            dataPosT.rename(columns={'HUM_ID':'chr'}, inplace=True)
            dataPosT['allReads']=dataPosT['reads_'+str(minLen)]
            dataPos=pd.merge(left=dataPos,right=dataPosT, how='right', on=["comb",
                                                                          "split",
                                                                          "chr",
                                                                          "R",
                                                                          "orient",
                                                                          "HUM_RE",
                                                                          "HIV_RS",
                                                                          "HUM_RS",
                                                                          "HIV_RE"],copy=True)
            
            dataPos.loc[dataPos['allReads_x'].isnull(),['allReads_x']]=dataPos.loc[dataPos['allReads_x'].isnull(),'allReads_x'].apply(lambda x: set())
            dataPos['allReads']=pd.Series([set1.union(set2) for set1, set2 in zip(dataPos['allReads_x'], dataPos['allReads_y'])])
            dataPos['tmp_reads_'+str(minLen)]=dataPos['reads_'+str(minLen)]-dataPos['allReads_x']
            dataPos.drop('reads_'+str(minLen),axis=1,inplace=True)
            dataPos.rename(columns={'tmp_reads_'+str(minLen):'reads_'+str(minLen)}, inplace=True)
            dataPos.loc[dataPos['reads_'+str(pastMinLen)].isnull(),['reads_'+str(pastMinLen)]]=dataPos.loc[dataPos['reads_'+str(pastMinLen)].isnull(),'reads_'+str(pastMinLen)].apply(lambda x: set())
            dataPos.loc[dataPos['count_'+str(pastMinLen)].isnull(),['count_'+str(pastMinLen)]]=dataPos.loc[dataPos['count_'+str(pastMinLen)].isnull(),'count_'+str(pastMinLen)].apply(lambda x: 0)
            dataPos.drop(['allReads_x','allReads_y'],axis=1,inplace=True)
            
            pastMinLen=minLen
            pastMinLens.append(minLen)
            
        for pastMinLen in pastMinLens:
            dataPos.loc[dataPos['reads_'+str(pastMinLen)].isnull(),['reads_'+str(pastMinLen)]]=dataPos.loc[dataPos['reads_'+str(pastMinLen)].isnull(),'reads_'+str(pastMinLen)].apply(lambda x: set())
            dataPos.loc[dataPos['count_'+str(pastMinLen)].isnull(),['count_'+str(pastMinLen)]]=dataPos.loc[dataPos['count_'+str(pastMinLen)].isnull(),'count_'+str(pastMinLen)].apply(lambda x: 0)
            
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

def writeReadNames(path,dataPos,dirPath,fileName,outPath1,outPath2):
    readsFile=open(path,'w+')
    readsList=[]
    if ";" in dataPos:
        for x in dataPos.split(";"):
            readsList.append(x)
    else:
        readsList.append(dataPos)

    for QNAME in readsList:
        readsFile.write(QNAME+"\n")
    readsFile.close()

    scriptCMD="./extractReadsPerSplit.sh "+dirPath+" "+fileName+" "+path+" "+outPath1+" "+outPath2
    os.system(scriptCMD)

def getStats(data,baseName,outDir):
    numSpliceJunctions=0
    with open(outDir+"/hisat/"+baseName+".junctions") as f:
        for i, l in enumerate(f):
            pass
    numSplits=len(data)
    numReads=data["count"].sum()
    numSpliceJunctions=numSpliceJunctions+1
    return pd.DataFrame([[baseName,numSplits,numReads,numSpliceJunctions]],columns=["name","numSplits","numReads","numSpliceHIV"])

def allSamples(out,paths,end):
    patientCodes=set([x.split("/")[-1].split('.')[0].split("_")[0] for x in paths])

    data=pd.DataFrame([])

    for patient in patientCodes:
        for sampleFile in glob.glob(os.path.abspath(out)+"/"+patient+"*Pos"+end+".csv"):
            dfT=pd.read_csv(sampleFile)
            dfT["patientName"]=patient
            dfT["sampleName"]=sampleFile.split("/")[-1].split(".")[0]
            data=pd.concat([data,dfT])
    dropList=[x for x in list(data) if 'reads' in x]
    data=data.reset_index(drop=True).drop(dropList,axis=1)
    data.replace(np.nan,"",inplace=True)
    dropList=[x for x in list(data) if 'count' in x]
    data["count"]=data["totalCount"]
    data=data.drop(dropList,axis=1)

    aggregations = {
        'count': {
            'numSamples': 'count',
            'numReads': 'sum',
            'median': 'median',
            'mean': 'mean'
        }
    }

    df=pd.DataFrame(data.groupby(by=["comb","orient"])["count"].agg(aggregations)).reset_index()
    df.columns=df.columns.droplevel(0)
    columns=list(df)
    columns[0]="split"
    columns[1]="orient"
    df.columns=columns
    df[["mean","median","numReads","numSamples"]]=df[["mean","median","numReads","numSamples"]].astype(int)
    df=df.sort_values(by="numSamples",ascending=False).reset_index(drop=True)
    df["sampleNames"]=df.apply(lambda row: ";".join(list(set(data[data["comb"]==row["split"]]["sampleName"]))),axis=1)
    df["patientNames"]=df.apply(lambda row: ";".join(list(set(data[data["comb"]==row["split"]]["patientName"]))),axis=1)
    df["hum_nearest_3SS"]=df.apply(lambda row: ";".join(list(set(data[data["comb"]==row["split"]]["hum_nearest_3SS"]))),axis=1)
    df["hum_nearest_5SS"]=df.apply(lambda row: ";".join(list(set(data[data["comb"]==row["split"]]["hum_nearest_5SS"]))),axis=1)
    def countNegPos(row):
        l=list(data[data["comb"]==row["split"]]["sampleName"])
        countNeg=0
        countPos=0
        for el in l:
            if "neg" in el:
                countNeg=countNeg+1
            else:
                countPos=countPos+1
        return [countNeg,countPos]
    df[["numNeg","numPos"]]=pd.DataFrame([x for x in df.apply(lambda row: countNegPos(row),axis=1)])
    df.to_csv(os.path.abspath(out)+"/allSamples"+end+".csv",index=False)

# the following function is designed to combine the two dataframes together
def combineLocalFull(dataPosLocal,dataPosFull,minLens):
    def generateColumnsList(minLens):
        colList=['comb']
        for minLen in minLens:
            colList=colList+["count_"+str(minLen),"reads_"+str(minLen)]
        colList=colList+['orient','prim','chr','split','R','spanR1-R2','allReads','spanCount']
        return colList
    
    #the first thing to do is to add any missing reads and counts column and populate them with blank {} and 0 respectively
    def addMissing(dataPosLocal,dataPosFull):
        setColumnsLocal=set(list(dataPosLocal))
        setColumnsFull=set(list(dataPosFull))
        missingLocal=setColumnsFull-setColumnsLocal
        missingFull=setColumnsLocal-setColumnsFull
        for col in missingLocal:
            if 'reads' in col:
                dataPosLocal[col]=''
            if 'count' in col:
                dataPosLocal[col]=0
        for col in missingFull:
            if 'reads' in col:
                dataPosFull[col]=''
            if 'count' in col:
                dataPosFull[col]=0
        return dataPosLocal,dataPosFull
    
    dataPosLocal,dataPosFull=addMissing(dataPosLocal,dataPosFull)
    
    colList=generateColumnsList(minLens)

    data=pd.DataFrame([],columns=colList)

    setLocalPos=set(dataPosLocal["comb"])
    setFullPos=set(dataPosFull["comb"])
    intersect=setFullPos.intersection(setLocalPos)
    diff=setFullPos.symmetric_difference(setLocalPos)

    if len(intersect)>0:
        for el in intersect:
            df2=pd.DataFrame(0, index=np.arange(1), columns=colList)
            l1=set(list(dataPosLocal[dataPosLocal['comb']==el]['spanR1-R2'])[0])
            l2=set(list(dataPosFull[dataPosFull['comb']==el]['spanR1-R2'])[0])
            setL=l1.union(l2)
            l1a=set(list(dataPosLocal[dataPosLocal['comb']==el]['allReads'])[0])
            l2a=set(list(dataPosFull[dataPosFull['comb']==el]['allReads'])[0])
            setLa=l1a.union(l2a)
            df2[['comb','orient','prim','chr','split','R','spanR1-R2','allReads','spanCount']]=[el,
                                                            dataPosLocal[dataPosLocal["comb"]==el]["orient"].iloc[0],
                                                            10,
                                                            dataPosLocal[dataPosLocal["comb"]==el]["chr"].iloc[0],
                                                            dataPosLocal[dataPosLocal["comb"]==el]["split"].iloc[0],
                                                            dataPosLocal[dataPosLocal["comb"]==el]["R"].iloc[0],
                                                            ";".join(list(setL)),
                                                            ";".join(list(setLa)),
                                                            len(setL)]

            for minLen in minLens:
                l1=set(list(dataPosLocal[dataPosLocal['comb']==el]['reads_'+str(minLen)])[0])
                l2=set(list(dataPosFull[dataPosFull['comb']==el]['reads_'+str(minLen)])[0])
                setL=l1.union(l2)
                df2[['count_'+str(minLen),'reads_'+str(minLen)]]=[len(setL),";".join(list(setL))]

            data=data.append(df2)
            
        data=data.reset_index(drop=True)

    dataLocalPosDiff=dataPosLocal[(dataPosLocal['comb'].isin(diff))]
    if len(dataLocalPosDiff)>0:
        for minLen in minLens:
            dataLocalPosDiff["reads_"+str(minLen)]=dataLocalPosDiff.apply(lambda row: ";".join(list(row['reads_'+str(minLen)])),axis=1)
            dataLocalPosDiff['count_'+str(minLen)]=dataLocalPosDiff['reads_'+str(minLen)].str.split(";").str.len()
            dataLocalPosDiff.loc[(dataLocalPosDiff['reads_'+str(minLen)]==""),"count_"+str(minLen)] = 0
        dataLocalPosDiff['spanR1-R2']=dataLocalPosDiff.apply(lambda row: ";".join(list(row['spanR1-R2'])),axis=1)
        dataLocalPosDiff['allReads']=dataLocalPosDiff.apply(lambda row: ";".join(list(row['allReads'])),axis=1)
    dataLocalPosDiff["prim"]=1
    dataFullPosDiff=dataPosFull[(dataPosFull['comb'].isin(diff))]
    if len(dataFullPosDiff)>0:
        for minLen in minLens:
            dataFullPosDiff["reads_"+str(minLen)]=dataFullPosDiff.apply(lambda row: ";".join(list(row['reads_'+str(minLen)])),axis=1)
            dataFullPosDiff['count_'+str(minLen)]=dataFullPosDiff['reads_'+str(minLen)].str.split(";").str.len()
            dataFullPosDiff.loc[(dataFullPosDiff['reads_'+str(minLen)]==""),"count_"+str(minLen)] = 0
        dataFullPosDiff['spanR1-R2']=dataFullPosDiff.apply(lambda row: ";".join(list(row['spanR1-R2'])),axis=1)
        dataFullPosDiff['allReads']=dataFullPosDiff.apply(lambda row: ";".join(list(row['allReads'])),axis=1)
    dataFullPosDiff["prim"]=0
    data=pd.concat([data,dataLocalPosDiff,dataFullPosDiff])
    data.reset_index(drop=True)
    
    colList=[]
    for minLen in minLens:
        colList=colList+["count_"+str(minLen)]
    data['totalCount']=data[colList].sum(axis=1)

    data.replace("","-",inplace=True)
    # for minLen in minLens:
    #     data.replace({"reads_"+str(minLen): {set([]):"-"},
    #                   "count_"+str(minLen): {set([]):0},
    #                   "splitR1-R2": {set([]):'-'}} ,inplace=True)
    return data

# after this step is done
# need to incorporate changes to the merging algorithm (when merging dataPosFull and dataPosLocal)
# also verify that counts after merging are correct
# perhaps would be reasonable to do the following:
# 1. only return sets of reads when building separate dataPosFull and dataPosLocal
# 2. count number of reads after merging the dataPosFull and dataPosLocal dataFrames together

#Another aspect to consider:
#when using multiple hiv references it might complicate things to guarantee that the hiv id of the spanning region is the same as that of the identified chimera
#thus we are making the assumption that the references are roughly the same in terms of the positioning of genomic regions
#however it would be best to add another note somewhere as to which particular gi number in the multiple reference the alignment belongs

def addSpanOld(row,dataSep):
    dataSpan=pd.DataFrame([])
    humR=row['R']
    hivR='R1'
    if humR=='R1':
        hivR='R2'
    if 'hiv:hum' in row['orient']:
        #&(dataSep[r+'HIV_RE']>row['HUM_RE'])&(dataSep[r+'HUM_RE']>row['HUM_RE'])
        dataSpan=dataSep[dataSep['HIV'].str.contains('sep'+humR)]  #check that the hum alignment is on the same side
        dataSpan=dataSpan[~(dataSpan[hivR+'HIV_ID']==0)] # check that the hiv is on the same side
        dataSpan=dataSpan[dataSpan[humR+'HUM_ID']==row['chr']] # check that the chromosomes match
        dataSpan=dataSpan[(dataSpan[hivR+'HIV_RS']-row['HIV_RS']>-500)&(dataSpan[hivR+'HIV_RS']-row['HIV_RS']<0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataSpan=dataSpan[(dataSpan[humR+'HUM_RE']-row['HUM_RE']<500)&(dataSpan[humR+'HUM_RE']-row['HUM_RE']>0)] # check that the end of the hum is after the end of the hum in the grouped dataframe
#         dataSpan=dataSpan[] # check that the distance between the beginning of the hiv and the end of the hum is reasonable. this step can likely be left out
        if len(dataSpan)>0:
            return [set(dataSpan["QNAME"]),len(dataSpan)] # return set of reads and count of reads
    if 'hum:hiv' in row['orient']:
        dataSpan=dataSep[dataSep['HIV'].str.contains('sep'+humR)]  #check that the hum alignment is on the same side
        dataSpan=dataSpan[~(dataSpan[hivR+'HIV_ID']==0)] # check that the hiv is on the same side
        dataSpan=dataSpan[dataSpan[humR+'HUM_ID']==row['chr']] # check that the chromosomes match
        dataSpan=dataSpan[(dataSpan[hivR+'HIV_RE']-row['HIV_RE']<500)&(dataSpan[hivR+'HIV_RE']-row['HIV_RE']>0)] # check that the start of the hiv is before the start of the hiv in the grouped dataframe
        dataSpan=dataSpan[(dataSpan[humR+'HUM_RS']-row['HUM_RS']>-500)&(dataSpan[humR+'HUM_RS']-row['HUM_RS']<0)] # check that the end of the hum is after the end of the hum in the grouped dataframe
#         dataSpan=dataSpan[] # check that the distance between the beginning of the hiv and the end of the hum is reasonable. this step can likely be left out
        if len(dataSpan)>0:
            return [set(dataSpan["QNAME"]),len(dataSpan)] # return set of reads and count of reads

    return [{''},0]

# this function will produce a dataframe with information grouped by the SpliceSites
# those reads that do not contain a valid spliceSite shall be discarded
def groupBySpliceSites(df):
    df.drop(['HIV_RE',
             'HIV_RS',
             'HUM_RE',
             'HUM_RS',
             'count_20',
             'count_30',
             'count_40',
             'orient',
             'prim',
             'reads_20',
             'reads_30',
             'reads_40',
             'spanCount',
             'spanR1-R2'],axis=1,inplace=True)
    aggregations={
        'allReads':{
            'groupsCount':'count',
            'allReads': lambda x: set(x)
        },
        'totalCount':{
            'allReads_count':'sum'
        }
    }
    dfg=pd.DataFrame(df.groupby(by=["hum_nearest_5SS","hum_nearest_3SS"])[["comb","R","allReads","totalCount"]].agg(aggregations)).reset_index()
    dfg.sort_values(by='groupsCount',ascending=False)

# this function is to replace the add.sh script
def addAnnotation(row,baseName,outDir):
    #need to add joint reads column to each dataframe which will contain all reads from reads_40 ..30..etc and spanReads
    #then can use this column below to find the first read (first read must always be with the greatest available minLen)

    firstRead=row['allReads'].split(';')[0]
    hum_nearest_5SS='-'
    hum_nearest_3SS='-'
    if row['prim']==0 or row['prim']==10:
        if row['R']=='R1':
            SS5_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".full.txt"+" | awk -F '\t' '{print $4}'"
            hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
            SS3_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".full.txt"+" | awk -F '\t' '{print $5}'"
            hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)
        else:
            SS5_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".full.txt"+" | awk -F '\t' '{print $9}'"
            hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
            SS3_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".full.txt"+" | awk -F '\t' '{print $10}'"
            hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)
    else:
        if row['R']=='R1':
            SS5_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".chim.txt"+" | awk -F '\t' '{print $4}'"
            hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
            SS3_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".chim.txt"+" | awk -F '\t' '{print $5}'"
            hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)
        else:
            SS5_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".chim.txt"+" | awk -F '\t' '{print $9}'"
            hum_nearest_5SS=subprocess.check_output(SS5_CMD,shell=True)
            SS3_CMD="grep '"+firstRead+"' "+outDir+"/splices/"+baseName+".chim.txt"+" | awk -F '\t' '{print $10}'"
            hum_nearest_3SS=subprocess.check_output(SS3_CMD,shell=True)

    return [hum_nearest_5SS,hum_nearest_3SS]

def wrapper(outDir,baseName,dirPath,fileName,minLenList,end,args):
    print(">>>>>>>>>>>>>   Begin analyses")
    # load data from local alignments
    dataLocalHIV = pd.read_csv(outDir+"/localAlignments/"+baseName+".chim.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    dataLocalHUM = pd.read_csv(outDir+"/localAlignments/"+baseName+".chim.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    # load data from full alignments
    dataFullHIV = pd.read_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    dataFullHUM = pd.read_csv(outDir+"/fullAlignments/"+baseName+".full.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    print("<<<<<<<<<<<<<   Done reading in") 
    outDirPOS=outDir+"/Positions/"+baseName
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDirPOS)):
        os.mkdir(os.path.abspath(outDirPOS))

    if (len(dataLocalHIV)==0 or len(dataLocalHUM)==0) and (len(dataFullHIV)==0):
        return

    data=pd.DataFrame([],columns=['comb','split','reads','count','orient','prim'])

    if len(dataLocalHIV)>0 and len(dataFullHIV)>0 and len(dataLocalHUM)>0:
        print('826')
        # extract flag information
        extractFlagBits(dataLocalHIV)
        extractFlagBits(dataLocalHUM)
        extractFlagBits(dataFullHIV)
        extractFlagBits(dataFullHUM)
        # extract start and end for both template and reference
        dataLocalHIV=extractStartEnd(dataLocalHIV)
        dataLocalHUM=extractStartEnd(dataLocalHUM)
        dataFullHIV=extractStartEnd(dataFullHIV)
        dataFullHUM=extractStartEnd(dataFullHUM)
        dataLocalHUM,dataLocalHIV=filterReads(dataLocalHUM,dataLocalHIV)
        if len(dataLocalHUM)==0:
            return

        # dataLocalHIV["lenAlign"]=dataLocalHIV.apply(lambda row: len(row["SEQ"]),axis=1)
        # dataLocalHUM["lenAlign"]=dataLocalHUM.apply(lambda row: len(row["SEQ"]),axis=1)
        #calculate the percent aligned (num bp aligned/total read length bp)
        # dataLocalHIV["percentAlign"]=(dataLocalHIV["Template_end"]-dataLocalHIV["Template_start"])/dataLocalHIV["lenAlign"]
        # dataLocalHUM["percentAlign"]=(dataLocalHUM["Template_end"]-dataLocalHUM["Template_start"])/dataLocalHUM["lenAlign"]
        dataLocal=pd.DataFrame(dataLocalHIV["QNAME"]).reset_index(drop=True)
        dataLocalHUM=dataLocalHUM.reset_index(drop=True)
        dataLocalHIV=dataLocalHIV.reset_index(drop=True)
        dataLocal=createData(dataLocal,dataLocalHUM,dataLocalHIV)
        dataLocalHUM.to_csv(outDir+"/localAlignments/"+baseName+".chim.hum"+end+".csv",index=False)
        dataLocalHIV.to_csv(outDir+"/localAlignments/"+baseName+".chim.hiv"+end+".csv",index=False)

        dataLocal.replace('', np.nan,inplace=True)
        dataLocal.fillna(0,inplace=True)
        dataLocal["overlapR1"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR1(row),axis=1))
        dataLocal["overlapR2"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR2(row),axis=1))
        dataLocal["HIV"]=dataLocal.apply(lambda row: leftRight(row),axis=1)
        dataLocal.to_csv(outDir+"/"+baseName+".chim"+end+".csv",index=False)
        # drop duplicated reads - preserve first occurence
        dataLocal.drop_duplicates(inplace=True)
        dataPosLocal=findSupport(dataLocal,minLenList)
        if len(dataPosLocal)>0:
            dataPosLocal=addSpan(dataLocal,dataPosLocal)
            dataPosLocal.loc[dataPosLocal['spanCount'].isnull(),['spanCount']]=dataPosLocal.loc[dataPosLocal['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
            dataPosLocal.loc[dataPosLocal['spanR1-R2'].isnull(),['spanR1-R2']]=dataPosLocal.loc[dataPosLocal['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())

        dataFullHUM,dataFullHIV=filterReads(dataFullHUM,dataFullHIV)
        if not len(dataFullHUM)==0:
            print('867')
            dataFull=pd.DataFrame(dataFullHIV["QNAME"]).reset_index(drop=True)
            dataFullHUM=dataFullHUM.reset_index(drop=True)
            dataFullHIV=dataFullHIV.reset_index(drop=True)
            dataFull=createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum"+end+".csv",index=False)
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv"+end+".csv",index=False)

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full"+end+".csv",index=False)
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataPosFull=findSupport(dataFull,minLenList)
            if len(dataPosFull)>0:
                dataPosFull=addSpan(dataFull,dataPosFull)
                # dataPosFull[['spanR1-R2','spanCount']]=pd.DataFrame([x for x in dataPosFull.apply(lambda row: addSpan(row,dataFull[dataFull['HIV'].str.contains('sep')]),axis=1)])
                dataPosFull.loc[dataPosFull['spanCount'].isnull(),['spanCount']]=dataPosFull.loc[dataPosFull['spanCount'].isnull(),'spanCount'].apply(lambda x: 0)
                dataPosFull.loc[dataPosFull['spanR1-R2'].isnull(),['spanR1-R2']]=dataPosFull.loc[dataPosFull['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())
            
            data=combineLocalFull(dataPosLocal,dataPosFull,minLenList)

            # colList=[]
            # for minLen in minLenList:
            #     colList=colList+["count_"+str(minLen)]
            # data['totalCount']=data[colList].sum(axis=1).astype(int)

            # data['totalCount']=data['totalCount']+data['spanCount']
            
            # #lets now also add a column for all the reads that belong to a achimeric position including the reads which span r1-r2
            # def joinReads(row,colList):
            #     colList=['reads_'+str(x) for x in colList]+['spanR1-R2']
            #     reads=[]
            #     for col in colList:
            #         if not pd.isnull(row[col]):
            #             reads.append(row[col])
            #     return ';'.join(reads)

            # colList=[int(x.split('_')[-1]) for x in list(data) if 'reads_' in x]
            # colList.sort(reverse=True)
            # data['allReads']=data['allReads'].apply(lambda x: ';'.join(x))
            # colList=['reads_'+str(x) for x in colList]+['spanR1-R2']
            # reads=[]
            # for col in colList:
            #     print(col)
            #     data[col]=data[col].apply(lambda x: ';'.join(x) if len(x)>0 or pd.isnull(x) else '-')

            #Now add the nearest human splice site for each chimeric position
            # data[['hum_nearest_5SS','hum_nearest_3SS']]=pd.DataFrame([x for x in data.apply(lambda row: addAnnotation(row,baseName,outDir),axis=1)])
            data.drop(["split","chr","HIV_RE","HIV_RS","HUM_RE","HUM_RS"],axis=1).sort_values(by='totalCount',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos"+end+".csv",index=False)

            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            childPIDS=[]
            #  cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            #  cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            #  os.system(cmdR1)
            #  os.system(cmdR2)
            #  print("write read names")
            #  data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            #  os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            #  os.remove(dirPath+"/"+baseName+"_R2_001.fastq")

        else:
            print('936')
            if len(dataPosLocal)>0:
                for minLen in minLenList:
                    dataPosLocal["reads_"+str(minLen)]=dataPosLocal.apply(lambda row: ";".join(list(row['reads_'+str(minLen)])),axis=1)
                    dataPosLocal['count_'+str(minLen)]=dataPosLocal['reads_'+str(minLen)].str.split(";").str.len()
                    dataPosLocal.loc[(dataPosLocal['reads_'+str(minLen)]==""),"count_"+str(minLen)] = 0

            dataPosLocal["prim"]=1
            data=pd.concat([data,dataPosLocal])
            data.reset_index(drop=True)

            colList=[]
            for minLen in minLenList:
                colList=colList+["count_"+str(minLen)]
            data['totalCount']=data[colList].sum(axis=1)

            data['allReads']=data['allReads'].apply(lambda x: ';'.join(x))
            data.loc[data['spanCount'].isnull(),['spanCount']]=0
            data.loc[data['spanR1-R2'].isnull(),['spanR1-R2']]=data.loc[data['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())
            data['spanR1-R2']=data['spanR1-R2'].apply(lambda x: ';'.join(x) if len(x)>0 or pd.isnull(x) else '-')

            #Now add the nearest human splice site for each chimeric position
            # data[['hum_nearest_5SS','hum_nearest_3SS']]=pd.DataFrame([x for x in data.apply(lambda row: addAnnotation(row,baseName,outDir),axis=1)])
            data.drop(["split","chr","HIV_RE","HIV_RS","HUM_RE","HUM_RS"],axis=1).sort_values(by='totalCount',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos"+end+".csv",index=False)
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            childPIDS=[]
            #  cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            #  cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            #  os.system(cmdR1)
            #  os.system(cmdR2)
            #  print("write read names")
            #  data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            #  os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            #  os.remove(dirPath+"/"+baseName+"_R2_001.fastq")

    if len(dataLocalHIV)==0 and len(dataFullHIV)>0:
        print('987')
        # extract flag information
        extractFlagBits(dataFullHIV)
        extractFlagBits(dataFullHUM)
        # extract start and end for both template and reference
        dataFullHIV=extractStartEnd(dataFullHIV)
        dataFullHUM=extractStartEnd(dataFullHUM)

        dataFullHUM,dataFullHIV=filterReads(dataFullHUM,dataFullHIV)
        if not len(dataFullHUM)==0:
            # dataFullHIV["lenAlign"]=dataFullHIV.apply(lambda row: len(row["SEQ"]),axis=1)
            # dataFullHUM["lenAlign"]=dataFullHUM.apply(lambda row: len(row["SEQ"]),axis=1)
            #calculate the percent aligned (num bp aligned/total read length bp)
            # dataFullHIV["percentAlign"]=(dataFullHIV["Template_end"]-dataFullHIV["Template_start"])/dataFullHIV["lenAlign"]
            # dataFullHUM["percentAlign"]=(dataFullHUM["Template_end"]-dataFullHUM["Template_start"])/dataFullHUM["lenAlign"]
            dataFull=pd.DataFrame(dataFullHUM["QNAME"]).reset_index(drop=True)
            dataFullHUM=dataFullHUM.reset_index(drop=True)
            dataFullHIV=dataFullHIV.reset_index(drop=True)
            dataFull=createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum"+end+".csv",index=False)
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv"+end+".csv",index=False)

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full"+end+".csv",index=False)
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataPosFull=findSupport(dataFull,minLenList)
            data=pd.DataFrame([])
            if len(dataPosFull)>0:
                data=addSpan(dataFull,dataPosFull)

            data['allReads']=data['allReads'].apply(lambda x: ';'.join(x))
            data.loc[data['spanCount'].isnull(),['spanCount']]=0
            data.loc[data['spanR1-R2'].isnull(),['spanR1-R2']]=data.loc[data['spanR1-R2'].isnull(),'spanR1-R2'].apply(lambda x: set())
            data['spanR1-R2']=data['spanR1-R2'].apply(lambda x: ';'.join(x) if len(x)>0 or pd.isnull(x) else '-')
            
            for minLen in minLenList:
                data["reads_"+str(minLen)]=data.apply(lambda row: ";".join(list(row['reads_'+str(minLen)])),axis=1)
                data['count_'+str(minLen)]=data['reads_'+str(minLen)].str.split(";").str.len()
                data.loc[(data['reads_'+str(minLen)]==""),"count_"+str(minLen)] = 0
                
            data.replace("","-",inplace=True)
            data.reset_index(drop=True)
            data["prim"]=0

            colList=[]
            for minLen in minLenList:
                colList=colList+["count_"+str(minLen)]
            data['totalCount']=data[colList].sum(axis=1)
            #save the bed file of the positions
            # data[['chr','HUM_RS','HUM_RE']].to_csv(outDir+"/beds/"+baseName+".bed",sep='\t',header=False,index=False)
            # bedCMD='bedtools intersect -a '+outDir+'/beds/'+baseName+'.bed'+' -b '+args.annotation+' -wo > '+outDir+'/beds/'+baseName+'.bed.out'
            # os.system(bedCMD)

            #Now add the nearest human splice site for each chimeric position
            # data[['hum_nearest_5SS','hum_nearest_3SS']]=pd.DataFrame([x for x in data.apply(lambda row: addAnnotation(row,baseName,outDir),axis=1)])
            data.drop(["split","chr","HIV_RE","HIV_RS","HUM_RE","HUM_RS"],axis=1).sort_values(by='totalCount',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos"+end+".csv",index=False)
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            # childPIDS=[]
            #  cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            #  cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            #  os.system(cmdR1)
            #  os.system(cmdR2)
            #  print("write read names")
            #  data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+end+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            #  os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            #  os.remove(dirPath+"/"+baseName+"_R2_001.fastq")
    
    return 1

#this function is responsible for identifying the relevant annotation line from the bedtoold intersect output for each position in the dataPos
def findAnnotation(row,intersection):
    prefOrder=["exon",
               "CDS",
               "mRNA",
               "gene",
               "ncRNA",
               "transcript",
               "primary_transcript",
               "rRNA",
               "tRNA",
               "locus"]


# 1. add HIV strain information along the chromosome to the comb column of the final dataframe
# 2. group by the splice site
# 3. incorporate splice site information gathering into the main script instead of running perl separately
# 4. currently novel sites with lower minLen of alignment are discarded/not included as part of the final dataframe
#       need to be included and higher minLen where absent needs to be approprietaly marked as such
# 5. have a logger ready for the following tasks:
#       - log which minLen values were discarded
#       - log general stats
# 6. add R1-R2 split information. Could be achieved as follows:
#       - somewhere R1-R2/R2-R1 should be added
#       - get the subset of the dataFrame from the dataLocal/dataFull which has the appropriate R1-R2/R2-R1 label
#       - create a new column in the dataPos for R1-R2 support
#       - add read names from the subset DataFrame to this new column in case the split position happens to be somewhere within the region spanned by the R1-R2 split

def testT(outDir,baseName,dirPath,fileName,minLenList,end,args):
    data=pd.read_csv(outDir+"/"+baseName+"_Pos"+end+".csv")
    #lets now also add a column for all the reads that belong to a achimeric position including the reads which span r1-r2
    def joinReads(row,colList):
        colList=['reads_'+str(x) for x in colList]+['spanR1-R2']
        reads=[]
        for col in colList:
            if not pd.isnull(row[col]):
                reads.append(row[col])
        return ';'.join(reads)

    colList=[int(x.split('_')[-1]) for x in list(data) if 'reads_' in x]
    colList.sort(reverse=True)
    data['allReads']=data.apply(lambda row: joinReads(row,colList),axis=1)

    #Now add the nearest human splice site for each chimeric position
    data[['hum_nearest_5SS','hum_nearest_3SS']]=pd.DataFrame([x for x in data.apply(lambda row: addAnnotation(row,baseName,outDir),axis=1)])
    data.to_csv(outDir+"/"+baseName+"_Pos"+end+".csv",index=False)

#yet another idea:
#do rum the remove duplication scripts on the alignments and perform the analysis on the cleaned data

#it seems the add.sh is pretty slow as well
#perhaps should timeit the execution of the custom python method against the add.sh
#would be best not to apply any external methods or as little as possible

#Theres something wrong with the total counts - spancounts seem to be included but some count_s are not

def main(args):
    print(args.minLen)
    end=""
    if args.end==True:
        end='.no_dup'

    for file in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
        fullPath=os.path.abspath(file)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])
        baseName="_R1".join(fileName.split("_R1")[:-1])
        scriptCMD="./kraken.sh "+dirPath+" "+fileName+" "+args.out+" "+args.krakenDB+" "+args.hivDB+" "+args.humDB+" "+args.annotation+" "+str(args.threads)

        setInter=['154_pos_12_S4',
                 '154_pos_17_S5',
                 '154_pos_18_S6',
                 '154_pos_21_S7',
                 '154_pos_6_S1',
                 '154_pos_7_S2',
                 '154_pos_8_S3']

        # setDone=['154_pos_7_S2',
        #          '154_pos_21_S7',
        #          '154_pos_17_S5']

        setDone=['218_pos_2_S11',
                 '348_pos_3_S15']
        
        if (baseName in setInter) and (baseName not in setDone):
            print("Analyzing: ",baseName)
            if args.shell==True:
                print("Running main shell script")
                os.system(scriptCMD)
            # if os.path.exists(os.path.abspath(os.path.abspath(args.out)+"/"+baseName+"_Pos"+end+".csv")):
                # testT(os.path.abspath(args.out),baseName,dirPath,fileName,args.minLen,end,args)
            resultsRow=wrapper(os.path.abspath(args.out),baseName,dirPath,fileName,args.minLen,end,args)

    #once the add.sh is replaced by the actual function - do one run with both included in order to
    #compare the results and make sure the new method is working properly
    #verification should be done as follows:
    #1. check if anything that is identified in the add.sh column is also present in the same way in the custom column
    #2. if any differences are observed - check why - likely a mistake
    #3. then check for any which are identified by the custom method but not by the add.sh
    
    os.system("./add.sh "+os.path.abspath(args.out))
    # paths=glob.glob(os.path.abspath(args.out)+"/*Pos"+end+".csv")
    # allSamples(args.out,paths,end)
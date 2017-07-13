#!/usr/bin/env python
#./post2.py -o ./out -i ./data/YHo060517 -k ./customDB -v ./refs/hiv89.6/hiv89.6 -g ./refs/bowtie2/hg38
import pandas as pd
import numpy as np
import numba
import os
import shutil
import sys
import glob
import time
import itertools
import argparse
import scipy
pd.set_option('display.max_columns', None)
import matplotlib as mpl
mpl.use('Agg')
# from mpl_toolkits.axes_grid1 import host_subplot
# import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
# %matplotlib inline
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
#     cigarList=["".join(x) for _, x in itertools.groupby(cigarTxt, key=str.isdigit)]
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
#     cigarList=["".join(x) for _, x in itertools.groupby(cigarTxt, key=str.isdigit)]
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
#     cigarList=["".join(x) for _, x in itertools.groupby(cigarTxt, key=str.isdigit)]
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
#     cigarList=["".join(x) for _, x in itertools.groupby(cigarTxt, key=str.isdigit)]
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
            t.append("sep:2")
            
        
    if not row["R2HUM_ID"]==0:            
        if not row["R2HIV_ID"]==0:
            if row['R2HIV_TS']>row["R2HUM_TS"]:
                t.append("R2:right")
            else:
                t.append("R2:left")
        else:
            t.append("sep:1")
            
    return "-".join(t)

def overlapR1(row):
    return len(set(range(row["R1HUM_TS"],row["R1HUM_TE"])).intersection(set(range(row["R1HIV_TS"],row["R1HIV_TE"]))))

def overlapR2(row):
    return len(set(range(row["R2HUM_TS"],row["R2HUM_TE"])).intersection(set(range(row["R2HIV_TS"],row["R2HIV_TE"]))))

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
    data0["Reference_end"]=data0.END+data0.POS

    data=pd.concat([data16,data0]).reset_index()
    data.drop(["index","SEQ_LEN","CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
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
    data[["R1HUM_TS","R1HUM_TE","R1HUM_ID","R1HUM_RS","R1HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)])>0 else [0,0,"",0,0],axis=1)])
    data[["R2HUM_TS","R2HUM_TE","R2HUM_ID","R2HUM_RS","R2HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)])>0 else [0,0,"",0,0],axis=1)])
    data[["R1HIV_TS","R1HIV_TE","R1HIV_ID","R1HIV_RS","R1HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)])>0 else [0,0,"",0,0],axis=1)])
    data[["R2HIV_TS","R2HIV_TE","R2HIV_ID","R2HIV_RS","R2HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)])>0 else [0,0,"",0,0],axis=1)])
    data["HUM_AL_MAP"]=data.apply(lambda row: alMap(row,dataHUM,dataHIV),axis=1)
    data["HIV_AL_MAP"]=data.apply(lambda row: alMap(row,dataHIV,dataHUM),axis=1)

def findSupportOld(data): # deprecated
    # get all local chimeric reads where both hiv and hum are on the R1 strand the R2 is empty for HIV
    dataR1=data[~(data["R1HUM_TE"]==0)&~(data["R1HIV_TE"]==0)&(data["R2HIV_TE"]==0)]
    dataR2=data[~(data["R2HUM_TE"]==0)&~(data["R2HIV_TE"]==0)&(data["R1HIV_TE"]==0)]
    dataSepHivR1=data[~(data["R2HUM_TE"]==0)&~(data["R1HIV_TE"]==0)&(data["R2HIV_TE"]==0)&(data["R1HUM_TE"]==0)]
    dataSepHivR2=data[~(data["R1HUM_TE"]==0)&~(data["R2HIV_TE"]==0)&(data["R1HIV_TE"]==0)&(data["R2HUM_TE"]==0)]
    allSep=list(dataR1["QNAME"])+list(dataR2["QNAME"])+list(dataSepHivR1["QNAME"])+list(dataSepHivR2["QNAME"])
    dataOther=data[data["QNAME"].isin(set(data["QNAME"]).difference(set(allSep)))]

    # find support for R1 right
    dataR1Right=dataR1[dataR1["HIV"].str.contains("R1:right")]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(str)+":"+dataR1Right["ins"].astype(str)+":"+dataR1Right['R1HIV_RS'].astype(str)
    dataPosR1Right=pd.DataFrame([])
    dataPosR1Right[["split","count"]]=pd.DataFrame(dataR1Right.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Right)==0:
        dataPosR1Right["readsLocal"]=dataPosR1Right.apply(lambda row: list(dataR1Right[dataR1Right["split"]==row["split"]]["QNAME"]),axis=1)
        dataPosR1Right["Read:orient"]="R1:right"

    # find support for R1 left
    dataR1Left=dataR1[dataR1["HIV"].str.contains("R1:left")]
    dataR1Left
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(str)+":"+dataR1Left["ins"].astype(str)+":"+dataR1Left['R1HUM_RS'].astype(str)
    dataPosR1Left=pd.DataFrame([])
    dataPosR1Left[["split","count"]]=pd.DataFrame(dataR1Left.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Left)==0:
        dataPosR1Left["readsLocal"]=dataPosR1Left.apply(lambda row: list(dataR1Left[dataR1Left["split"]==row["split"]]["QNAME"]),axis=1)
        dataPosR1Left["Read:orient"]="R1:left"

    # find support for R2 right
    dataR2Right=dataR2[dataR2["HIV"].str.contains("R2:right")]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(str)+":"+dataR2Right["ins"].astype(str)+":"+dataR2Right['R2HIV_RS'].astype(str)
    dataPosR2Right=pd.DataFrame([])
    dataPosR2Right[["split","count"]]=pd.DataFrame(dataR2Right.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Right)==0:
        dataPosR2Right["readsLocal"]=dataPosR2Right.apply(lambda row: list(dataR2Right[dataR2Right["split"]==row["split"]]["QNAME"]),axis=1)
        dataPosR2Right["Read:orient"]="R2:right"

    # find support for R2 left
    dataR2Left=dataR2[dataR2["HIV"].str.contains("R2:left")]
    dataR2Left
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(str)+":"+dataR2Left["ins"].astype(str)+":"+dataR2Left['R2HUM_RS'].astype(str)
    dataPosR2Left=pd.DataFrame([])
    dataPosR2Left[["split","count"]]=pd.DataFrame(dataR2Left.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Left)==0:
        dataPosR2Left["readsLocal"]=dataPosR2Left.apply(lambda row: list(dataR2Left[dataR2Left["split"]==row["split"]]["QNAME"]),axis=1)
        dataPosR2Left["Read:orient"]="R2:left"

    frames=[dataPosR1Right,dataPosR1Left,dataPosR2Right,dataPosR2Left]
    return pd.concat(frames).reset_index()

def findSupportOverlap(data): # current method

    # find support for R1 right
    dataR1Right=data[data["HIV"].str.contains("R1:right")]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(str)+":"+dataR1Right["ins"].astype(str)+":"+dataR1Right['R1HIV_RS'].astype(str)
    dataPosR1Right=pd.DataFrame([])
    dataPosR1Right[["split","count"]]=pd.DataFrame(dataR1Right.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Right)==0:
        dataPosR1Right["readsLocal"]=dataPosR1Right.apply(lambda row: set(list(dataR1Right[dataR1Right["split"]==row["split"]]["QNAME"])),axis=1)
        dataPosR1Right["Read:orient"]="R1:right"

    # find support for R1 left
    dataR1Left=data[data["HIV"].str.contains("R1:left")]
    dataR1Left
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(str)+":"+dataR1Left["ins"].astype(str)+":"+dataR1Left['R1HUM_RS'].astype(str)
    dataPosR1Left=pd.DataFrame([])
    dataPosR1Left[["split","count"]]=pd.DataFrame(dataR1Left.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Left)==0:
        dataPosR1Left["readsLocal"]=dataPosR1Left.apply(lambda row: set(list(dataR1Left[dataR1Left["split"]==row["split"]]["QNAME"])),axis=1)
        dataPosR1Left["Read:orient"]="R1:left"

    # find support for R2 right
    dataR2Right=data[data["HIV"].str.contains("R2:right")]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(str)+":"+dataR2Right["ins"].astype(str)+":"+dataR2Right['R2HIV_RS'].astype(str)
    dataPosR2Right=pd.DataFrame([])
    dataPosR2Right[["split","count"]]=pd.DataFrame(dataR2Right.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Right)==0:
        dataPosR2Right["readsLocal"]=dataPosR2Right.apply(lambda row: set(list(dataR2Right[dataR2Right["split"]==row["split"]]["QNAME"])),axis=1)
        dataPosR2Right["Read:orient"]="R2:right"

    # find support for R2 left
    dataR2Left=data[data["HIV"].str.contains("R2:left")]
    dataR2Left
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(str)+":"+dataR2Left["ins"].astype(str)+":"+dataR2Left['R2HUM_RS'].astype(str)
    dataPosR2Left=pd.DataFrame([])
    dataPosR2Left[["split","count"]]=pd.DataFrame(dataR2Left.groupby(by=["split"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Left)==0:
        dataPosR2Left["readsLocal"]=dataPosR2Left.apply(lambda row: set(list(dataR2Left[dataR2Left["split"]==row["split"]]["QNAME"])),axis=1)
        dataPosR2Left["Read:orient"]="R2:left"

    frames=[dataPosR1Right,dataPosR1Left,dataPosR2Right,dataPosR2Left]
    return pd.concat(frames).reset_index()

def findSupport(data): # current method
    # find support for R1 right
    dataR1Right=data[data["HIV"].str.contains("R1:right")]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(str)+":"+dataR1Right['R1HIV_RS'].astype(str)
    dataR1Right["comb"]=dataR1Right.split+"_"+dataR1Right.R1HUM_ID
    dataPosR1Right=pd.DataFrame([])
    dataPosR1Right[["comb","split","chr","count"]]=pd.DataFrame(dataR1Right.groupby(by=["comb","split","R1HUM_ID"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Right)==0:
        dataPosR1Right["readsLocal"]=dataPosR1Right.apply(lambda row: set(list(dataR1Right[(dataR1Right["split"]==row["split"])&(dataR1Right["R1HUM_ID"]==row["chr"])]["QNAME"])),axis=1)
        dataPosR1Right["Read:orient"]="R1:right"

    # find support for R1 left
    dataR1Left=data[data["HIV"].str.contains("R1:left")]
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(str)+":"+dataR1Left['R1HUM_RS'].astype(str)
    dataR1Left["comb"]=dataR1Left.split+"_"+dataR1Left.R1HUM_ID
    dataPosR1Left=pd.DataFrame([])
    dataPosR1Left[["comb","split","chr","count"]]=pd.DataFrame(dataR1Left.groupby(by=["comb","split","R1HUM_ID"])["QNAME"].count()).reset_index()
    if not len(dataPosR1Left)==0:
        dataPosR1Left["readsLocal"]=dataPosR1Left.apply(lambda row: set(list(dataR1Left[(dataR1Left["split"]==row["split"])&(dataR1Left["R1HUM_ID"]==row["chr"])]["QNAME"])),axis=1)
        dataPosR1Left["Read:orient"]="R1:left"

    # find support for R2 right
    dataR2Right=data[data["HIV"].str.contains("R2:right")]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(str)+":"+dataR2Right['R2HIV_RS'].astype(str)
    dataR2Right["comb"]=dataR2Right.split+"_"+dataR2Right.R2HUM_ID
    dataPosR2Right=pd.DataFrame([])
    dataPosR2Right[["comb","split","chr","count"]]=pd.DataFrame(dataR2Right.groupby(by=["comb","split","R2HUM_ID"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Right)==0:
        dataPosR2Right["readsLocal"]=dataPosR2Right.apply(lambda row: set(list(dataR2Right[(dataR2Right["split"]==row["split"])&(dataR2Right["R2HUM_ID"]==row["chr"])]["QNAME"])),axis=1)
        dataPosR2Right["Read:orient"]="R2:right"

    # find support for R2 left
    dataR2Left=data[data["HIV"].str.contains("R2:left")]
    dataR2Left
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(str)+":"+dataR2Left['R2HUM_RS'].astype(str)
    dataR2Left["comb"]=dataR2Left.split+"_"+dataR2Left.R2HUM_ID
    dataPosR2Left=pd.DataFrame([])
    dataPosR2Left[["comb","split","chr","count"]]=pd.DataFrame(dataR2Left.groupby(by=["comb","split","R2HUM_ID"])["QNAME"].count()).reset_index()
    if not len(dataPosR2Left)==0:
        dataPosR2Left["readsLocal"]=dataPosR2Left.apply(lambda row: set(list(dataR2Left[(dataR2Left["split"]==row["split"])&(dataR2Left["R2HUM_ID"]==row["chr"])]["QNAME"])),axis=1)
        dataPosR2Left["Read:orient"]="R2:left"

    frames=[dataPosR1Right,dataPosR1Left,dataPosR2Right,dataPosR2Left]
    return pd.concat(frames).reset_index()

def approxCloseness(split1,split2):
    res=abs(int(split1.split(":")[0])-int(split2.split(":")[0]))+abs(int(split1.split(":")[2])-int(split2.split(":")[2]))
    try:
        res2=math.log(res)
    except:
        return 0
    return res2

def getSet(row):
    tmpLocalQNAME=list(row["readsLocal"])
    return tmpLocalQNAME

def writeReadNames(path,dataPos,dirPath,fileName,outPath):
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

    scriptCMD="./extractReadsPerSplit.sh "+dirPath+" "+fileName+" "+path+" "+outPath
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

def wrapper(outDir,baseName,dirPath,fileName):

    # load data from local alignments
    dataLocalHIV = pd.read_csv(outDir+"/localAlignments/"+baseName+".chim.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    dataLocalHUM = pd.read_csv(outDir+"/localAlignments/"+baseName+".chim.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    # load data from full alignments
    dataFullHIV = pd.read_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    dataFullHUM = pd.read_csv(outDir+"/fullAlignments/"+baseName+".full.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    
    outDirPOS=outDir+"/Positions/"+baseName
    if not os.path.exists(os.path.abspath(outDir+"/Positions/")):
        os.mkdir(os.path.abspath(outDir+"/Positions/"))
    if not os.path.exists(os.path.abspath(outDirPOS)):
        os.mkdir(os.path.abspath(outDirPOS))

    if (len(dataLocalHIV)==0 or len(dataLocalHUM)==0) and (len(dataFullHIV)==0):
        return

    data=pd.DataFrame([],columns=['comb','split','readsLocal','count','Read:orient','prim'])

    if len(dataLocalHIV)>0 and len(dataFullHIV)>0 and len(dataLocalHUM)>0:

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

        dataLocalHIV["lenAlign"]=dataLocalHIV.apply(lambda row: len(row["SEQ"]),axis=1)
        dataLocalHUM["lenAlign"]=dataLocalHUM.apply(lambda row: len(row["SEQ"]),axis=1)
        #calculate the percent aligned (num bp aligned/total read length bp)
        dataLocalHIV["percentAlign"]=(dataLocalHIV["Template_end"]-dataLocalHIV["Template_start"])/dataLocalHIV["lenAlign"]
        dataLocalHUM["percentAlign"]=(dataLocalHUM["Template_end"]-dataLocalHUM["Template_start"])/dataLocalHUM["lenAlign"]
        dataLocal=pd.DataFrame(dataLocalHUM["QNAME"]).reset_index().drop("index",axis=1)
        dataLocalHUM=dataLocalHUM.reset_index().drop("index",axis=1)
        dataLocalHIV=dataLocalHIV.reset_index().drop("index",axis=1)
        createData(dataLocal,dataLocalHUM,dataLocalHIV)
        dataLocalHUM.to_csv(outDir+"/localAlignments/"+baseName+".chim.hum.csv")
        dataLocalHIV.to_csv(outDir+"/localAlignments/"+baseName+".chim.hiv.csv")

        dataLocal.replace('', np.nan,inplace=True)
        dataLocal.fillna(0,inplace=True)
        dataLocal["overlapR1"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR1(row),axis=1))
        dataLocal["overlapR2"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR2(row),axis=1))
        dataLocal["HIV"]=dataLocal.apply(lambda row: leftRight(row),axis=1)
        dataLocal.to_csv(outDir+"/"+baseName+".chim.csv")
        # drop duplicated reads - preserve first occurence
        dataLocal.drop_duplicates(inplace=True)
        dataPosLocal=findSupport(dataLocal)
        dataPosLocal=dataPosLocal.drop("index",axis=1)

        dataFullHUM,dataFullHIV=filterReads(dataFullHUM,dataFullHIV)
        if not len(dataFullHUM)==0:
            # dataFullHIV["lenAlign"]=dataFullHIV.apply(lambda row: len(row["SEQ"]),axis=1)
            # dataFullHUM["lenAlign"]=dataFullHUM.apply(lambda row: len(row["SEQ"]),axis=1)
            #calculate the percent aligned (num bp aligned/total read length bp)
            # dataFullHIV["percentAlign"]=(dataFullHIV["Template_end"]-dataFullHIV["Template_start"])/dataFullHIV["lenAlign"]
            # dataFullHUM["percentAlign"]=(dataFullHUM["Template_end"]-dataFullHUM["Template_start"])/dataFullHUM["lenAlign"]
            dataFull=pd.DataFrame(dataFullHUM["QNAME"]).reset_index().drop("index",axis=1)
            dataFullHUM=dataFullHUM.reset_index().drop("index",axis=1)
            dataFullHIV=dataFullHIV.reset_index().drop("index",axis=1)
            createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum.csv")
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.csv")

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full.csv")
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataPosFull=findSupport(dataFull)
            dataPosFull=dataPosFull.drop("index",axis=1)
            # dataPosLocal["set"]=np.nan
            # if len(dataPosLocal)>0:
            #     dataPosLocal["set"]=dataPosLocal.apply(lambda row: getSet(row),axis=1)
            # dataPosFull["set"]=np.nan
            # if len(dataPosFull)>0:
            #     dataPosFull["set"]=dataPosFull.apply(lambda row: getSet(row),axis=1)
            setLocalPos=set(dataPosLocal["comb"])
            setFullPos=set(dataPosFull["comb"])
            intersect=setFullPos.intersection(setLocalPos)
            if len(intersect)>0:
                for el in intersect:
                    l1=set(list(dataPosLocal[dataPosLocal["comb"]==el]["readsLocal"])[0])
                    l2=set(list(dataPosFull[dataPosFull["comb"]==el]["readsLocal"])[0])
                    newSet=l1.union(l2)
                    df2 = pd.DataFrame([[el,len(newSet),";".join(list(newSet)),dataPosLocal[dataPosLocal["comb"]==el]["Read:orient"].iloc[0],10]],columns=['comb','split','count','readsLocal','Read:orient','prim'])
                    data=data.append(df2)

                data=data.reset_index().drop("index",axis=1)
            dataLocalPosDiff=dataPosLocal[~(dataPosLocal['split'].isin(intersect))]
            if len(dataLocalPosDiff)>0:
                dataLocalPosDiff["readsLocal"]=dataLocalPosDiff.apply(lambda row: ";".join(list(row['readsLocal'])),axis=1)
            # dataLocalPosDiff=dataLocalPosDiff.drop(["set"],axis=1)
            dataLocalPosDiff["prim"]=1
            dataFullPosDiff=dataPosFull[~(dataPosFull['split'].isin(intersect))]
            if len(dataFullPosDiff)>0:
                dataFullPosDiff["readsLocal"]=dataFullPosDiff.apply(lambda row: ";".join(list(row['readsLocal'])),axis=1)
            # dataFullPosDiff=dataFullPosDiff.drop(["set"],axis=1)
            dataFullPosDiff["prim"]=0
            data=pd.concat([data,dataLocalPosDiff,dataFullPosDiff])
            data.to_csv(outDir+"/"+baseName+"_Pos.csv")

            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["count"]))+".txt",row["readsLocal"],dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["count"]))+".fq"),axis=1)

        else:
            if len(dataPosLocal)>0:
                dataPosLocal["readsLocal"]=dataPosLocal.apply(lambda row: ";".join(list(row['readsLocal'])),axis=1)
            # dataPosLocal=dataPosLocal.drop(["set"],axis=1)
            dataPosLocal["prim"]=1
            data=pd.concat([data,dataPosLocal])
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["count"]))+".txt",row["readsLocal"],dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["count"]))+".fq"),axis=1)
            data.to_csv(outDir+"/"+baseName+"_Pos.csv")

    if len(dataLocalHIV)==0 and len(dataFullHIV)>0:
        # extract flag information
        extractFlagBits(dataFullHIV)
        extractFlagBits(dataFullHUM)
        # extract start and end for both template and reference
        dataFullHIV=extractStartEnd(dataFullHIV)
        dataFullHUM=extractStartEnd(dataFullHUM)

        dataFullHUM,dataFullHIV=filterReads(dataFullHUM,dataFullHIV)
        if not len(dataFullHUM)==0:
            dataFullHIV["lenAlign"]=dataFullHIV.apply(lambda row: len(row["SEQ"]),axis=1)
            dataFullHUM["lenAlign"]=dataFullHUM.apply(lambda row: len(row["SEQ"]),axis=1)
            #calculate the percent aligned (num bp aligned/total read length bp)
            dataFullHIV["percentAlign"]=(dataFullHIV["Template_end"]-dataFullHIV["Template_start"])/dataFullHIV["lenAlign"]
            dataFullHUM["percentAlign"]=(dataFullHUM["Template_end"]-dataFullHUM["Template_start"])/dataFullHUM["lenAlign"]
            dataFull=pd.DataFrame(dataFullHUM["QNAME"]).reset_index().drop("index",axis=1)
            dataFullHUM=dataFullHUM.reset_index().drop("index",axis=1)
            dataFullHIV=dataFullHIV.reset_index().drop("index",axis=1)
            createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum.csv")
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.csv")

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full.csv")
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataPosFull=findSupport(dataFull)

            # dataPosFull["set"]=np.nan
            if len(dataPosFull)>0:
                # dataPosFull["set"]=dataPosFull.apply(lambda row: getSet(row),axis=1)
                dataPosFull["readsLocal"]=dataPosFull.apply(lambda row: ";".join(list(row['readsLocal'])),axis=1)
                # dataPosFull=dataPosFull.drop(["set"],axis=1)
            dataPosFull["prim"]=0
            data=pd.concat([data,dataPosFull])
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["count"]))+".txt",row["readsLocal"],dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["count"]))+".fq"),axis=1)
            data.to_csv(outDir+"/"+baseName+"_Pos.csv")

    # return getStats(data,baseName,outDir)
    return 1

def mainRun(args):

    # finalStatsDF=pd.DataFrame([],columns=["name","numSplits","numReads","numSpliceHIV"])
    for file in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
        fullPath=os.path.abspath(file)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])

        baseName="_R1".join(fileName.split("_R1")[:-1])
        scriptCMD="./kraken.sh "+dirPath+" "+fileName+" "+args.out+" "+args.krakenDB+" "+args.hivDB+" "+args.humDB
        if not baseName in ["Y430_pos_12_S51",
                            "PH029_pos_3_S39",
                            "Y430_pos_3_S42",
                            "Y430_pos_10_S49",
                            "Y430_neg_14_S17",
                            "Y111_pos_3_S62",
                            "Y354_neg_1_S34",
                            "PH029_neg_1_S1",
                            "Y159_neg_2_S21",
                            "Y430_neg_5_S8",
                            "Y159_pos_2_S57",
                            "Y430_pos_4_S43",
                            "Y111_pos_6_S65",
                            "PH029_neg_3_S3",
                            "Y430_neg_13_S16",
                            "Y111_pos_4_S63",
                            "Y241_pos_2_S68",
                            "Y430_pos_6_S45",
                            "Y111_pos_5_S64",
                            "Y430_pos_5_S44",
                            "Y430_neg_8_S11",
                            "Y241_neg_2_S32",
                            "Y354_neg_2_S35",
                            "Y159_neg_1_S20",
                            "Y111_neg_4_S27",
                            "Y159_pos_1_S56",
                            "Y241_pos_1_S67",
                            "Y430_pos_11_S50",
                            "Y430_neg_2_S5",
                            "Y111_neg_3_S26",
                            "Y430_pos_8_S47",
                            "Y241_neg_1_S31",
                            "Y430_neg_11_S14",
                            "PH029_pos_2_S38",
                            "Y430_pos_15_S54",
                            "Y354_pos_1_S69",
                            "Y159_pos_4_S59",
                            "Y430_neg_15_S18",
                            "Y354_neg_3_S36",
                            "Y430_neg_9_S12",
                            "Y430_neg_6_S9",
                            "Y111_pos_1_S60",
                            "Y111_neg_6_S29",
                            "Y430_pos_1_S40",
                            "Y159_pos_3_S58",
                            "Y111_neg_2_S25",
                            "Y241_neg_3_S33",
                            "Y159_neg_4_S23",
                            "Y111_neg_1_S24",
                            "Y430_neg_7_S10",
                            "Y111_pos_2_S61",
                            "Y111_neg_5_S28",
                            "Y430_pos_2_S41",
                            "Y111_neg_7_S30",
                            "Y430_pos_13_S52",
                            "PH029_pos_1_S37"]:
            # resultsRow=pd.DataFrame([])
            # if baseName in ["Y159_pos_2_S57"]:
            print(baseName)
            # os.system(scriptCMD)
            resultsRow=wrapper(os.path.abspath(args.out),baseName,dirPath,fileName)
            # finalStatsDF=pd.concat([finalStatsDF,resultsRow])

    # finalStatsDF=finalStatsDF.reset_index().drop("index",axis=1)
    # finalStatsDF.to_csv(os.path.abspath(args.out)+"/results.csv")

def main(argv):

    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('-i',
                                '--input',
                                required=True,
                                type=str,
                                help="directory which contains fastq.gz files")
    parser.add_argument('-o',
                                '--out',
                                required=False,
                                type=str,
                                default="./out",
                                help="output directory")
    parser.add_argument('-k',
                                '--krakenDB',
                                required=True,
                                type=str,
                                help="path to kraken database")
    parser.add_argument('-v',
                                '--hivDB',
                                required=True,
                                type=str,
                                help="path to the hiv reference")
    parser.add_argument('-g',
                                '--humDB',
                                required=True,
                                type=str,
                                help="path to the hg38 reference")

    parser.set_defaults(func=mainRun)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])
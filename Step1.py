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
    data[["R1HUM_TS","R1HUM_TE","R1HUM_ID","R1HUM_RS","R1HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)])>0 else [0,0,"",0,0],axis=1)])
    data[["R2HUM_TS","R2HUM_TE","R2HUM_ID","R2HUM_RS","R2HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)])>0 else [0,0,"",0,0],axis=1)])
    data[["R1HIV_TS","R1HIV_TE","R1HIV_ID","R1HIV_RS","R1HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)])>0 else [0,0,"",0,0],axis=1)])
    data[["R2HIV_TS","R2HIV_TE","R2HIV_ID","R2HIV_RS","R2HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)][["Template_start","Template_end","RNAME","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)])>0 else [0,0,"",0,0],axis=1)])
    data["HUM_AL_MAP"]=data.apply(lambda row: alMap(row,dataHUM,dataHIV),axis=1)
    data["HIV_AL_MAP"]=data.apply(lambda row: alMap(row,dataHIV,dataHUM),axis=1)

def filterOverlapCombine(data,minLenHUM,minLenHIV):
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
              "HUM_AL_MAP",
              "HIV_AL_MAP",
              "overlapR1",
              "overlapR2",
              "HIV"]

    # R1 right
    dataR1Right=data[data["HIV"].str.contains("R1:right")]
    dataR1Right=dataR1Right[~((dataR1Right["R1HUM_TS"]<dataR1Right["R1HIV_TS"])&(dataR1Right["R1HUM_TE"]>dataR1Right["R1HIV_TE"]))]
    dataR1Right=dataR1Right[~((dataR1Right["R1HIV_TS"]<dataR1Right["R1HUM_TS"])&(dataR1Right["R1HIV_TE"]>dataR1Right["R1HUM_TE"]))]
    dataR1Right["ins"]=dataR1Right["R1HIV_TS"]-dataR1Right["R1HUM_TE"]
    dataR1Right["split"]=dataR1Right['R1HUM_RE'].astype(str)+":"+dataR1Right['R1HIV_RS'].astype(str)
    dataR1Right["comb"]=dataR1Right.split+"@"+dataR1Right.R1HUM_ID
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
    dataR1Right=dataR1Right[(dataR1Right["HIV_TE"]-dataR1Right["HIV_TS"]-dataR1Right["overlap"])>minLenHIV]
    dataR1Right=dataR1Right[(dataR1Right["HUM_TE"]-dataR1Right["HUM_TS"]-dataR1Right["overlap"])>minLenHUM]
    dataR1Right["R"]="R1"
    dataR1Right.drop(dropList,axis=1,inplace=True)

    # R1 left
    dataR1Left=data[data["HIV"].str.contains("R1:left")]
    dataR1Left=dataR1Left[~((dataR1Left["R1HUM_TS"]<dataR1Left["R1HIV_TS"])&(dataR1Left["R1HUM_TE"]>dataR1Left["R1HIV_TE"]))]
    dataR1Left=dataR1Left[~((dataR1Left["R1HIV_TS"]<dataR1Left["R1HUM_TS"])&(dataR1Left["R1HIV_TE"]>dataR1Left["R1HUM_TE"]))]
    dataR1Left["ins"]=dataR1Left["R1HUM_TS"]-dataR1Left["R1HIV_TE"]
    dataR1Left["split"]=dataR1Left['R1HIV_RE'].astype(str)+":"+dataR1Left['R1HUM_RS'].astype(str)
    dataR1Left["comb"]=dataR1Left.split+"@"+dataR1Left.R1HUM_ID
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
    dataR1Left=dataR1Left[(dataR1Left["HIV_TE"]-dataR1Left["HIV_TS"]-dataR1Left["overlap"])>minLenHIV]
    dataR1Left=dataR1Left[(dataR1Left["HUM_TE"]-dataR1Left["HUM_TS"]-dataR1Left["overlap"])>minLenHUM]
    dataR1Left["R"]="R1"
    dataR1Left.drop(dropList,axis=1,inplace=True)

    # R2 right
    dataR2Right=data[data["HIV"].str.contains("R2:right")]
    dataR2Right=dataR2Right[~((dataR2Right["R2HUM_TS"]<dataR2Right["R2HIV_TS"])&(dataR2Right["R2HUM_TE"]>dataR2Right["R2HIV_TE"]))]
    dataR2Right=dataR2Right[~((dataR2Right["R2HIV_TS"]<dataR2Right["R2HUM_TS"])&(dataR2Right["R2HIV_TE"]>dataR2Right["R2HUM_TE"]))]
    dataR2Right["ins"]=dataR2Right["R2HIV_TS"]-dataR2Right["R2HUM_TE"]
    dataR2Right["split"]=dataR2Right['R2HUM_RE'].astype(str)+":"+dataR2Right['R2HIV_RS'].astype(str)
    dataR2Right["comb"]=dataR2Right.split+"@"+dataR2Right.R2HUM_ID
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
    dataR2Right=dataR2Right[(dataR2Right["HIV_TE"]-dataR2Right["HIV_TS"]-dataR2Right["overlap"])>minLenHIV]
    dataR2Right=dataR2Right[(dataR2Right["HUM_TE"]-dataR2Right["HUM_TS"]-dataR2Right["overlap"])>minLenHUM]
    dataR2Right["R"]="R2"
    dataR2Right.drop(dropList,axis=1,inplace=True)

    # R2 left
    dataR2Left=data[data["HIV"].str.contains("R2:left")]
    dataR2Left=dataR2Left[~((dataR2Left["R2HUM_TS"]<dataR2Left["R2HIV_TS"])&(dataR2Left["R2HUM_TE"]>dataR2Left["R2HIV_TE"]))]
    dataR2Left=dataR2Left[~((dataR2Left["R2HIV_TS"]<dataR2Left["R2HUM_TS"])&(dataR2Left["R2HIV_TE"]>dataR2Left["R2HUM_TE"]))]
    dataR2Left["ins"]=dataR2Left["R2HUM_TS"]-dataR2Left["R2HIV_TE"]
    dataR2Left["split"]=dataR2Left['R2HIV_RE'].astype(str)+":"+dataR2Left['R2HUM_RS'].astype(str)
    dataR2Left["comb"]=dataR2Left.split+"@"+dataR2Left.R2HUM_ID
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
    dataR2Left=dataR2Left[(dataR2Left["HIV_TE"]-dataR2Left["HIV_TS"]-dataR2Left["overlap"])>minLenHIV]
    dataR2Left=dataR2Left[(dataR2Left["HUM_TE"]-dataR2Left["HUM_TS"]-dataR2Left["overlap"])>minLenHUM]
    dataR2Left["R"]="R2"
    dataR2Left.drop(dropList,axis=1,inplace=True)

    frames=[dataR1Right,dataR1Left,dataR2Right,dataR2Left]
    df=pd.concat(frames).reset_index(drop=True)
    # data["overlap"]=data.apply(lambda row: 0-row["ins"] if row["overlap"]==0 else row["overlap"],axis=1)
    return df

def findSupport(dataHC, dataLC):
    # this function should do the following:
    # 1. group data by split position
    # 2. add a column of high confidence support reads from the first parameter DF
    # 3. add a column of low confidence support reads from the second parameter DF

    # First group by split and extract high confidence support reads
    dataPos=pd.DataFrame([])
    dataPos[["comb","split","chr","R","orient","countHC"]]=pd.DataFrame(dataHC.groupby(by=["comb","split","HUM_ID","R","orient"])["QNAME"].count()).reset_index()
    if not len(dataPos)==0:
        dataPos["readsHC"]=dataPos.apply(lambda row: set(list(dataHC[(dataHC["comb"]==row["comb"])&(dataHC["HUM_ID"]==row["chr"])]["QNAME"])),axis=1)
        dataPos["readsLC"]=dataPos.apply(lambda row: set(list(dataLC[(dataLC["comb"]==row["comb"])&(dataLC["HUM_ID"]==row["chr"])]["QNAME"])).symmetric_difference(row["readsHC"]),axis=1)
        dataPos["countLC"]=dataPos.readsLC.str.len()
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

def allSamples(out):
    paths=glob.glob(os.path.abspath(out)+"/*Pos.csv")
    patientCodes=set([x.split("/")[-1].split('.')[0].split("_")[0] for x in paths])

    data=pd.DataFrame([])

    for patient in patientCodes:
        for sampleFile in glob.glob(os.path.abspath(out)+"/"+patient+"*Pos.csv"):
            dfT=pd.read_csv(sampleFile)
            dfT["patientName"]=patient
            dfT["sampleName"]=sampleFile.split("/")[-1].split(".")[0]
            data=pd.concat([data,dfT])
    data=data.reset_index(drop=True).drop(["readsHC","readsLC"],axis=1)
    data.replace(np.nan,"",inplace=True)
    data["count"]=data["countLC"]+data["countHC"]
    data=data.drop(["countHC","countLC"],axis=1)

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
    df.to_csv(os.path.abspath(out)+"/allSamples.csv",index=False)

def wrapper(outDir,baseName,dirPath,fileName,minLenHC_HUM,minLenLC_HUM,minLenHC_HIV,minLenLC_HIV):

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

    data=pd.DataFrame([],columns=['comb','split','reads','count','orient','prim'])

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

        # dataLocalHIV["lenAlign"]=dataLocalHIV.apply(lambda row: len(row["SEQ"]),axis=1)
        # dataLocalHUM["lenAlign"]=dataLocalHUM.apply(lambda row: len(row["SEQ"]),axis=1)
        #calculate the percent aligned (num bp aligned/total read length bp)
        # dataLocalHIV["percentAlign"]=(dataLocalHIV["Template_end"]-dataLocalHIV["Template_start"])/dataLocalHIV["lenAlign"]
        # dataLocalHUM["percentAlign"]=(dataLocalHUM["Template_end"]-dataLocalHUM["Template_start"])/dataLocalHUM["lenAlign"]
        dataLocal=pd.DataFrame(dataLocalHUM["QNAME"]).reset_index(drop=True)
        dataLocalHUM=dataLocalHUM.reset_index(drop=True)
        dataLocalHIV=dataLocalHIV.reset_index(drop=True)
        createData(dataLocal,dataLocalHUM,dataLocalHIV)
        dataLocalHUM.to_csv(outDir+"/localAlignments/"+baseName+".chim.hum.csv",index=False)
        dataLocalHIV.to_csv(outDir+"/localAlignments/"+baseName+".chim.hiv.csv",index=False)

        dataLocal.replace('', np.nan,inplace=True)
        dataLocal.fillna(0,inplace=True)
        dataLocal["overlapR1"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR1(row),axis=1))
        dataLocal["overlapR2"]=pd.DataFrame(dataLocal.apply(lambda row: overlapR2(row),axis=1))
        dataLocal["HIV"]=dataLocal.apply(lambda row: leftRight(row),axis=1)
        dataLocal.to_csv(outDir+"/"+baseName+".chim.csv",index=False)
        # drop duplicated reads - preserve first occurence
        dataLocal.drop_duplicates(inplace=True)
        dataLocalHC=filterOverlapCombine(dataLocal,minLenHC_HUM,minLenHC_HIV)
        dataLocalLC=filterOverlapCombine(dataLocal,minLenLC_HUM,minLenLC_HIV)
        dataPosLocal=findSupport(dataLocalHC,dataLocalLC)

        dataFullHUM,dataFullHIV=filterReads(dataFullHUM,dataFullHIV)
        if not len(dataFullHUM)==0:
            dataFull=pd.DataFrame(dataFullHUM["QNAME"]).reset_index(drop=True)
            dataFullHUM=dataFullHUM.reset_index(drop=True)
            dataFullHIV=dataFullHIV.reset_index(drop=True)
            createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum.csv",index=False)
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.csv",index=False)

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full.csv",index=False)
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataFullHC=filterOverlapCombine(dataFull,minLenHC_HUM,minLenHC_HIV)
            dataFullLC=filterOverlapCombine(dataFull,minLenLC_HUM,minLenLC_HIV)
            dataPosFull=findSupport(dataFullHC,dataFullLC)
            setLocalPos=set(dataPosLocal["comb"])
            setFullPos=set(dataPosFull["comb"])
            intersect=setFullPos.intersection(setLocalPos)
            diff=setFullPos.symmetric_difference(setLocalPos)
            if len(intersect)>0:
                for el in intersect:
                    # new set of high confidence/support reads
                    l1HC=set(list(dataPosLocal[dataPosLocal["comb"]==el]["readsHC"])[0])
                    l2HC=set(list(dataPosFull[dataPosFull["comb"]==el]["readsHC"])[0])
                    setHC=l1HC.union(l2HC)
                    #new set of low confidence/support reads
                    l1LC=set(list(dataPosLocal[dataPosLocal["comb"]==el]["readsLC"])[0])
                    l2LC=set(list(dataPosFull[dataPosFull["comb"]==el]["readsLC"])[0])
                    setLC=l1LC.union(l2LC)
                    df2 = pd.DataFrame([[el,
                                        len(setHC),
                                        ";".join(list(setHC)),
                                        len(setLC),
                                        ";".join(list(setLC)),
                                        dataPosLocal[dataPosLocal["comb"]==el]["orient"].iloc[0],
                                        10,
                                        dataPosLocal[dataPosLocal["comb"]==el]["chr"].iloc[0],
                                        dataPosLocal[dataPosLocal["comb"]==el]["split"].iloc[0],
                                        dataPosLocal[dataPosLocal["comb"]==el]["R"].iloc[0]]],columns=['comb','countHC','readsHC','countLC','readsLC','orient','prim','chr','split','R'])
                    data=data.append(df2)

                data=data.reset_index(drop=True)
            dataLocalPosDiff=dataPosLocal[(dataPosLocal['comb'].isin(diff))]
            if len(dataLocalPosDiff)>0:
                dataLocalPosDiff["readsHC"]=dataLocalPosDiff.apply(lambda row: ";".join(list(row['readsHC'])),axis=1)
                dataLocalPosDiff["readsLC"]=dataLocalPosDiff.apply(lambda row: ";".join(list(row['readsLC'])),axis=1)
            dataLocalPosDiff["prim"]=1
            dataFullPosDiff=dataPosFull[(dataPosFull['comb'].isin(diff))]
            if len(dataFullPosDiff)>0:
                dataFullPosDiff["readsHC"]=dataFullPosDiff.apply(lambda row: ";".join(list(row['readsHC'])),axis=1)
                dataFullPosDiff["readsLC"]=dataFullPosDiff.apply(lambda row: ";".join(list(row['readsLC'])),axis=1)
            dataFullPosDiff["prim"]=0
            data=pd.concat([data,dataLocalPosDiff,dataFullPosDiff])
            data.reset_index(drop=True)
            data.drop(["count","reads","split","chr"],axis=1).sort_values(by='countHC',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos.csv",index=False)

            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            childPIDS=[]
            cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            os.system(cmdR1)
            os.system(cmdR2)
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            os.remove(dirPath+"/"+baseName+"_R2_001.fastq")

        else:
            if len(dataPosLocal)>0:
                dataPosLocal["readsHC"]=dataPosLocal.apply(lambda row: ";".join(list(row['readsHC'])),axis=1)
                dataPosLocal["readsLC"]=dataPosLocal.apply(lambda row: ";".join(list(row['readsLC'])),axis=1)
            dataPosLocal["prim"]=1
            data=pd.concat([data,dataPosLocal])
            data.reset_index(drop=True)
            data.drop(["count","reads","split","chr"],axis=1).sort_values(by='countHC',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos.csv",index=False)
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            childPIDS=[]
            cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            os.system(cmdR1)
            os.system(cmdR2)
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            os.remove(dirPath+"/"+baseName+"_R2_001.fastq")

    if len(dataLocalHIV)==0 and len(dataFullHIV)>0:
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
            createData(dataFull,dataFullHUM,dataFullHIV)
            dataFullHUM.to_csv(outDir+"/fullAlignments/"+baseName+".full.hum.csv",index=False)
            dataFullHIV.to_csv(outDir+"/fullAlignments/"+baseName+".full.hiv.csv",index=False)

            dataFull.replace('', np.nan,inplace=True)
            dataFull.fillna(0,inplace=True)
            dataFull["overlapR1"]=pd.DataFrame(dataFull.apply(lambda row: overlapR1(row),axis=1))
            dataFull["overlapR2"]=pd.DataFrame(dataFull.apply(lambda row: overlapR2(row),axis=1))
            dataFull["HIV"]=dataFull.apply(lambda row: leftRight(row),axis=1)
            dataFull.to_csv(outDir+"/"+baseName+".full.csv",index=False)
            # drop duplicated reads - preserve first occurence
            dataFull.drop_duplicates(inplace=True)
            dataFullHC=filterOverlapCombine(dataFull,minLenHC_HUM,minLenHC_HIV)
            dataFullLC=filterOverlapCombine(dataFull,minLenLC_HUM,minLenLC_HIV)
            dataPosFull=findSupport(dataFullHC,dataFullLC)

            if len(dataPosFull)>0:
                dataPosFull["readsHC"]=dataPosFull.apply(lambda row: ";".join(list(row['readsHC'])),axis=1)
                dataPosFull["readsLC"]=dataPosFull.apply(lambda row: ";".join(list(row['readsLC'])),axis=1)
            dataPosFull["prim"]=0
            data=pd.concat([data,dataPosFull])
            data.reset_index(drop=True)
            data.drop(["count","reads","split","chr"],axis=1).sort_values(by='countHC',ascending=False).reset_index(drop=True).to_csv(outDir+"/"+baseName+"_Pos.csv",index=False)
            if not os.path.exists(os.path.abspath(outDirPOS+"/fq/")):
                os.mkdir(os.path.abspath(outDirPOS+"/fq/"))

            childPIDS=[]
            cmdR1="zcat "+dirPath+"/"+baseName+"_R1_001.fastq.gz > "+dirPath+"/"+baseName+"_R1_001.fastq"
            cmdR2="zcat "+dirPath+"/"+baseName+"_R2_001.fastq.gz > "+dirPath+"/"+baseName+"_R2_001.fastq"
            os.system(cmdR1)
            os.system(cmdR2)
            data.apply(lambda row: writeReadNames(outDirPOS+"/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+".txt",";".join([row["readsHC"],row["readsLC"]]),dirPath,fileName,outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R1.fq",outDirPOS+"/fq/"+str(row["comb"])+"_"+str(int(row["countLC"]+int(row["countHC"])))+"_R2.fq"),axis=1)
            os.remove(dirPath+"/"+baseName+"_R1_001.fastq")
            os.remove(dirPath+"/"+baseName+"_R2_001.fastq")
    return 1

def main(args):
    minLenHC_HUM=int(args.minLen1.split(":")[0])
    minLenLC_HUM=int(args.minLen1.split(":")[1])
    minLenHC_HIV=int(args.minLen2.split(":")[0])
    minLenLC_HIV=int(args.minLen2.split(":")[1])

    fileList=["154_neg_2_S9",
    "154_neg_3_S10",
    "154_pos_18_S6",
    "154_pos_21_S7",
    "154_pos_7_S2",
    "218_neg_1_S17",
    "218_pos_1_S13",
    "218_pos_3_S15",
    "256_pos_1_S11",
    "348_pos_1_S20"]

    for file in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
        fullPath=os.path.abspath(file)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])

        baseName="_R1".join(fileName.split("_R1")[:-1])

        if not baseName in fileList:
            scriptCMD="./kraken.sh "+dirPath+" "+fileName+" "+args.out+" "+args.krakenDB+" "+args.hivDB+" "+args.humDB+" "+args.annotation+" "+str(args.threads)

            print("Analyzing: ",baseName)
            # os.system(scriptCMD)
            # resultsRow=wrapper(os.path.abspath(args.out),baseName,dirPath,fileName,minLenHC_HUM,minLenLC_HUM,minLenHC_HIV,minLenLC_HIV)

    # os.system("./add.sh "+os.path.abspath(args.out))
    allSamples(args.out)

    scriptCMD="./results.sh "+os.path.abspath(args.out)+" "+os.path.abspath(args.input)
    os.system(scriptCMD)

#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
import sys
import glob
import time
import itertools
import scipy

# The function below returns a dataframe which consists of paired-end reads
# where each read is confirmed to contain a segment from both viral and human genomes
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
def calcAlignmentStartEnd(cigar,forward,start):
    # first, break cigar into parseable entities
    cigarList=["".join(x) for _, x in itertools.groupby(cigar, key=str.isdigit)]
    idxfM=cigarList.index("M") #index of the first M
    idxlM=len(cigarList)-1-cigarList[::-1].index("M") #index of the last M
    start=start+sum(list(map(int,cigarList[:idxfM-1][::2])))
    end=start+sum(list(map(int,cigarList[idxfM-1:idxlM][::2])))-1
    return str(start)+":"+str(end)

def main():
    dataHIV = pd.read_csv("./108PI_pos_S2.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21'])
    dataHUM = pd.read_csv("./108PI_pos_S2.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21'])
    dataHIV["paired"]=dataHIV["FLAG"]               &1 #template having multiple segments in sequencing
    dataHIV["aligned2Mates"]=dataHIV["FLAG"]        &2 #each segment properly aligned according to the aligner
    dataHIV["unmappedCurr"]=dataHIV["FLAG"]         &4 #segment unmapped
    dataHIV["unmappedMate"]=dataHIV["FLAG"]         &8 #next segment in the template unmapped
    dataHIV["reversedCurr"]=dataHIV["FLAG"]         &16 #SEQ being reverse complemented
    dataHIV["reversedMate"]=dataHIV["FLAG"]         &32 #SEQ of the next segment in the template being reverse complemented
    dataHIV["firstRead"]=dataHIV["FLAG"]            &64 #the first segment in the template
    dataHIV["lastRead"]=dataHIV["FLAG"]             &128 #the last segment in the template
    dataHIV["secondaryAlignment"]=dataHIV["FLAG"]   &256 # secondary alignment
    dataHIV["noPassFilter"]=dataHIV["FLAG"]         &512 #not passing filters, such as platform/vendor quality controls
    dataHIV["PCRdup"]=dataHIV["FLAG"]               &1024 #PCR or optical duplicate
    dataHIV["suppAl"]=dataHIV["FLAG"]               &2048 #supplementary alignment
    dataHIV["Template_start:end"]=dataHIV.apply(lambda row: calcAlignmentStartEnd(row['CIGAR'],True,0),axis=1)
    dataHIV[["Template_start","Template_end"]]=dataHIV["Template_start:end"].str.split(':', expand=True).astype(int)
    dataHIV["Reference_start:end"]=dataHIV.apply(lambda row: calcAlignmentStartEnd(row['CIGAR'],True,row['POS']),axis=1)
    dataHIV[["Reference_start","Reference_end"]]=dataHIV["Reference_start:end"].str.split(':', expand=True).astype(int)
    
    dataHUM["paired"]=dataHUM["FLAG"]               &1 #template having multiple segments in sequencing
    dataHUM["aligned2Mates"]=dataHUM["FLAG"]        &2 #each segment properly aligned according to the aligner
    dataHUM["unmappedCurr"]=dataHUM["FLAG"]         &4 #segment unmapped
    dataHUM["unmappedMate"]=dataHUM["FLAG"]         &8 #next segment in the template unmapped
    dataHUM["reversedCurr"]=dataHUM["FLAG"]         &16 #SEQ being reverse complemented
    dataHUM["reversedMate"]=dataHUM["FLAG"]         &32 #SEQ of the next segment in the template being reverse complemented
    dataHUM["firstRead"]=dataHUM["FLAG"]            &64 #the first segment in the template
    dataHUM["lastRead"]=dataHUM["FLAG"]             &128 #the last segment in the template
    dataHUM["secondaryAlignment"]=dataHUM["FLAG"]   &256 # secondary alignment
    dataHUM["noPassFilter"]=dataHUM["FLAG"]         &512 #not passing filters, such as platform/vendor quality controls
    dataHUM["PCRdup"]=dataHUM["FLAG"]               &1024 #PCR or optical duplicate
    dataHUM["suppAl"]=dataHUM["FLAG"]               &2048 #supplementary alignment
    dataHUM["Template_start:end"]=dataHUM.apply(lambda row: calcAlignmentStartEnd(row['CIGAR'],True,0),axis=1)
    dataHUM[["Template_start","Template_end"]]=dataHUM["Template_start:end"].str.split(':', expand=True).astype(int)
    dataHUM["Reference_start:end"]=dataHUM.apply(lambda row: calcAlignmentStartEnd(row['CIGAR'],True,row['POS']),axis=1)
    dataHUM[["Reference_start","Reference_end"]]=dataHUM["Reference_start:end"].str.split(':', expand=True).astype(int)
    
    outCols=list(dataHUM)
    outCols.remove("Template_start:end")
    outCols.remove("Reference_start:end")
    dataHUM[outCols].to_csv("./dataHUM.csv")
    dataHIV[outCols].to_csv("./dataHIV.csv")

if __name__=="__main__":
    main()
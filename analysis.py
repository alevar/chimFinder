#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
import sys
import glob
import time
import itertools
import scipy
import argparse

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

def wrapper(outDir,baseName):
    dataHIV = pd.read_csv(outDir+"/localAlignments/"+baseName+".hiv.chim.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21'])
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
   
    dataHUM = pd.read_csv(outDir+"/localAlignments/"+baseName+".hum.chim.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21']) 
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


    #filtering the reads based on the flags:
    #remove all reads that belong to secondary or supplementary alignments and did not have PCR duplicates
    dataHUM=dataHUM[(dataHUM["secondaryAlignment"]==0)&(dataHUM["PCRdup"]==0)&(dataHUM["suppAl"]==0)&(dataHUM["noPassFilter"]==0)]
    dataHIV=dataHIV[(dataHIV["secondaryAlignment"]==0)&(dataHIV["PCRdup"]==0)&(dataHIV["suppAl"]==0)&(dataHIV["noPassFilter"]==0)]
    #now we can exclude paired-end reads which do not contain information from both references
    setHIV=set(dataHIV["QNAME"])
    setHUM=set(dataHUM["QNAME"])
    setQnames=setHIV.intersection(setHUM)
    dataHUM=dataHUM[dataHUM["QNAME"].isin(setQnames)]
    dataHIV=dataHIV[dataHIV["QNAME"].isin(setQnames)]
    #calculate full length of the read
    dataHIV["lenAlign"]=dataHIV.apply(lambda row: len(row["SEQ"]),axis=1)
    dataHUM["lenAlign"]=dataHUM.apply(lambda row: len(row["SEQ"]),axis=1)

    #calculate the percent aligned (num bp aligned/total read length bp)
    dataHIV["percentAlign"]=(dataHIV["Template_end"]-dataHIV["Template_start"])/dataHIV["lenAlign"]
    dataHUM["percentAlign"]=(dataHUM["Template_end"]-dataHUM["Template_start"])/dataHUM["lenAlign"]

    #deduce whether a read is aligned on the 5' or 3'
    #this should be done once the reads from two alignments are put together
    #that way we can look at the relative position of the alignments within reads
    #if the read aligns to both references we can then deduce the prime end of the read to which each alignment belongs
    #1  - 5'
    #-1 - 3'
    #0  - unavailable. Marked if the read only aligns to a single reference and not the other

    #==========================================
    #Let's now focus attention on the reads which align to both HIV and HUM reference within a single mate
    setHUMf=set(dataHUM[dataHUM["firstRead"]==64]["QNAME"]) #get all forward reads
    setHUMl=set(dataHUM[dataHUM["lastRead"]==128]["QNAME"]) #get all last reads
    setHIVf=set(dataHIV[dataHIV["firstRead"]==64]["QNAME"]) #get all forward reads
    setHIVl=set(dataHIV[dataHIV["lastRead"]==128]["QNAME"]) #get all last reads

    intersectF=setHUMf.intersection(setHIVf)
    intersectL=setHUMl.intersection(setHIVl)

    #for the reads identified run the following check:
    #if the hiv alignment is on the inside of the pair - unlikely
    #if on the outer should look for the neighboring reads

    #we may express the overlap between the hiv and human alignments within a read as follows:
    #+n - there is n nucleotides separating the two alignments
    #-n - the two alignments overlap by n nucleotides

    #==========================================
    #Now let's look at those reads where forward mate aligned to one reference and the reverse mate to the other
    cleanHUM=dataHUM[(dataHUM["paired"]==1)&(dataHUM["aligned2Mates"]==0)]
    cleanHUM_First=cleanHUM[(cleanHUM["firstRead"]==64)&(cleanHUM["lastRead"]==0)]
    cleanHUM_Last=cleanHUM[(cleanHUM["firstRead"]==0)&(cleanHUM["lastRead"]==128)]

    cleanHIV=dataHIV[(dataHIV["paired"]==1)&(dataHIV["aligned2Mates"]==0)]
    cleanHIV_First=cleanHIV[(cleanHIV["firstRead"]==64)&(cleanHIV["lastRead"]==0)]
    cleanHIV_Last=cleanHIV[(cleanHIV["firstRead"]==0)&(cleanHIV["lastRead"]==128)]
    #find mates within paired-end reads which aligned to both hiv and hum
    setHIV_First=set(cleanHIV_First["QNAME"])
    setHUM_First=set(cleanHUM_First["QNAME"])
    setHIV_Last=set(cleanHIV_Last["QNAME"])
    setHUM_Last=set(cleanHUM_Last["QNAME"])
    splitMatesF=setHIV_First.intersection(setHUM_Last)
    splitMatesF=splitMatesF.difference(intersectF)
    splitMatesF=splitMatesF.difference(intersectL)
    splitMatesL=setHIV_Last.intersection(setHUM_First)
    splitMatesL=splitMatesL.difference(intersectF)
    splitMatesL=splitMatesL.difference(intersectL)

    #==========================================
    #Once we identify the reads where each mate aligns to a single but different reference
    #(eg. forward mate aligned only to the human reference and the reverse mate aligned to the hiv reference)
    #we can look for previously identified split reads which are contained
    #within the region of the reference spanned by the pair
    #need to take into account the possibility that the split read is not contained in the region spanned
    #by the alignment but one of the paired-end reads themselves is also a split read.
    #Just need to look for the intersection between the split reads and current sets.
    
    outCols=list(dataHUM)
    outCols.remove("Template_start:end")
    outCols.remove("Reference_start:end")
    dataHUM[outCols].to_csv(outDir+"/localAlignments/"+baseName+".hum.chim.csv")
    dataHIV[outCols].to_csv(outDir+"/localAlignments/"+baseName+".hiv.chim.csv")

    data=pd.DataFrame(dataHUM["QNAME"]).reset_index().drop("index",axis=1)
    dataHUM.reset_index().drop("index",axis=1)
    dataHIV.reset_index().drop("index",axis=1)
    data[["R1HUM_TS","R1HUM_TE","R1HUM_RS","R1HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)][["Template_start","Template_end","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["firstRead"]==64)])>0 else [0,0,0,0],axis=1)])
    data[["R2HUM_TS","R2HUM_TE","R2HUM_RS","R2HUM_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)][["Template_start","Template_end","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHUM[(dataHUM["QNAME"]==row["QNAME"])&(dataHUM["lastRead"]==128)])>0 else [0,0,0,0],axis=1)])
    data[["R1HIV_TS","R1HIV_TE","R1HIV_RS","R1HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)][["Template_start","Template_end","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["firstRead"]==64)])>0 else [0,0,0,0],axis=1)])
    data[["R2HIV_TS","R2HIV_TE","R2HIV_RS","R2HIV_RE"]] = pd.DataFrame([x for x in data.apply(lambda row: dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)][["Template_start","Template_end","Reference_start","Reference_end"]].values.tolist()[0] if len(dataHIV[(dataHIV["QNAME"]==row["QNAME"])&(dataHIV["lastRead"]==128)])>0 else [0,0,0,0],axis=1)])
    data.to_csv(outDir+"/localAlignments/"+baseName+".chim.csv")

def mainRun(args):
    for file in glob.glob(os.path.abspath(args.input)+"/*R1_001.fastq.gz"):
    # for file in glob.glob(os.path.abspath(args.input)+"/*ht2"):
        fullPath=os.path.abspath(file)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])
        
        scriptCMD="./kraken.sh "+dirPath+" "+fileName+" "+args.out+" "+args.krakenDB+" "+args.hivDB+" "+args.humDB
        os.system(scriptCMD)

        baseName="_R1".join(fileName.split("_R1")[:-1])
        wrapper(os.path.abspath(args.out),baseName)


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
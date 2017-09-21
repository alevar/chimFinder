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

# get positions where n's occur
def getContigPos(consensus):
    consensus=consensus.replace("\n","")
    if "n" in consensus:
        #enumeration=[i for i,x in enumerate(consensus) if ((x=="n") and not i==len(consensus)-1 and not consensus[i+1]=="n") or (x=="n" and not consensus[i-1]=="n") and not i==0]
        #if not consensus[len(consensus)-1]=="n":
        #    enumeration.append(len(consensus)-1)
        #pairs=list(zip(enumeration[0::2],enumeration[1::2]))
        #return pairs
        return [(m.start(0)-1, m.end(0)) for m in re.finditer(re.compile(r'([^n]+)'), consensus)]
    else:
        return []

# filter out short distances.
# for instance if the threshold is set to 70 - do not count any contigPairs which such limited length
def filterDistances(contigPairs,threshold):
    for contigPair in contigPairs:
        if contigPair[1][0]-contigPair[0][1]<threshold:
            contigPairs.remove(contigPair)
    return contigPairs

# filter contigs smaller in length than set by a parameter
def filterLength(contigs,length):
    return [contig for contig in contigs if contig[1]-contig[0]>length]

# build all possible forward pairs of positions which could be potential links between splices
def buildPairs(contigs):
    return list(itertools.combinations(contigs, 2))

# trim positions to specified length on each side
def trimContigPos(contigsPaired,length,totalLen):
    newContigList=[]
    for contigPair in contigsPaired:
        ex1=contigPair[0][1]-contigPair[0][0]
        ex2=contigPair[1][1]-contigPair[1][0]
        if contigPair[0][1]>length and ex1>length and totalLen-contigPair[1][0]>length and ex2>length:
            newcontigPair=((contigPair[0][1]-length,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]<=length and totalLen-contigPair[1][0]>length and ex2>length:
            newcontigPair=((contigPair[0][0]+1,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][0]+length))
            newContigList.append(newcontigPair)

        if contigPair[0][1]>length and ex1>length and totalLen-contigPair[1][0]<=length:
            newcontigPair=((contigPair[0][1]-length,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][1]))
            newContigList.append(newcontigPair)

        if contigPair[0][1]<=length and totalLen-contigPair[1][0]<=length:
            newcontigPair=((contigPair[0][0]+1,contigPair[0][1]),(contigPair[1][0]+1,contigPair[1][1]))
            newContigList.append(newcontigPair)

    return newContigList

# extract string for each contig
def getStrings(coordinates,sequence):
    return [(sequence[cord[0][0]:cord[0][1]],sequence[cord[1][0]:cord[1][1]],cord[0][1],cord[1][0]) for cord in coordinates]

# build a possible spliced read
def getUnifiedRef(contigs):
    return [("".join([contigPair[0],contigPair[1]]),contigPair[2],contigPair[3],len(contigPair[0])) for contigPair in contigs]

# write splicesitedata to fasta files to be used as reference in alignment
def writeRefs(finalOut,finalOutUni,header,outDir,baseName):
    outDir=os.path.abspath(outDir)+"/spliceRefs/"
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    outDir=outDir+baseName
    # create a directory for the sample
    # in the directory save each possible splice fasta as endPosPre:startPosPost
    if not os.path.exists(os.path.abspath(outDir)):
        os.mkdir(os.path.abspath(outDir))

    # write unified splice sites first
    for spliceSite in finalOutUni:
        spliceSiteUniPath=os.path.abspath(outDir)+"/u-"+str(spliceSite[1])+":"+str(spliceSite[2])+"_"+str(spliceSite[3])+".fa"
        readsFile=open(spliceSiteUniPath,'w+')
        readsFile.write(header)
        readsFile.write(spliceSite[0])
        readsFile.close()

    for spliceSite in finalOut:
        spliceSitePath1=os.path.abspath(outDir)+"/i-1_"+str(spliceSite[2])+":"+str(spliceSite[3])+".fa"
        spliceSitePath2=os.path.abspath(outDir)+"/i-2_"+str(spliceSite[2])+":"+str(spliceSite[3])+".fa"
        readsFile1=open(spliceSitePath1,'w+')
        readsFile1.write(header)
        readsFile1.write(spliceSite[0])
        readsFile1.close()
        readsFile2=open(spliceSitePath2,'w+')
        readsFile2.write(header)
        readsFile2.write(spliceSite[1])
        readsFile2.close()

# run the shell script to build indeces and perform alignments
def indexAlign(args,finalOut,baseName,posName,posS):
    os.system("./spliceDetect.sh "+os.path.abspath(args.out)+" "+baseName+" "+args.input+" "+posName+" "+posS)

def confirmSplice(baseName,filePath,args,baseEnd):
    header=""
    finalOut=[]
    finalOutUni=[]
    count=0
    for line in open(filePath,"r"):
        if count==0:
            header=line
            # write first line to the new file
            count=count+1
        else:
            contigs=getContigPos(line)
            contigs=filterLength(contigs,args.cLen)
            if len(contigs)>0:
                contigsPaired=buildPairs(contigs)
                contigsPaired=filterDistances(contigsPaired,args.minLen)
                contigsPairedTrimmed=trimContigPos(contigsPaired,args.len,len(line))
                contigsStrs=getStrings(contigsPairedTrimmed,line)
                contigStrsUni=getUnifiedRef(contigsStrs)
                finalOut=contigsStrs
                finalOutUni=contigStrsUni

    writeRefs(finalOut,finalOutUni,header,args.out,baseName+baseEnd)

    os.system("./prepIndexAlign.sh "+os.path.abspath(args.out)+" "+baseName+" "+args.input+" "+baseEnd)
    for pos in finalOutUni:
        indexAlign(args,finalOut,baseName+baseEnd,str(pos[1])+":"+str(pos[2]),str(pos[3]))
    return [(x[1],x[2]) for x in finalOutUni]

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

def overlapLR(sample,base,outDir,pos,s):
    dataL=pd.read_csv(os.path.abspath(outDir)+"/splicing/"+sample+"/i-1_"+base+".sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    extractFlagBits(dataL)
    dataL=extractStartEnd(dataL)
    dataR=pd.read_csv(os.path.abspath(outDir)+"/splicing/"+sample+"/i-2_"+base+".sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    extractFlagBits(dataR)
    dataR=extractStartEnd(dataR)
    sL=set(list(dataL[pos-dataL["Reference_end"]<s]["QNAME"]))
    sR=set(list(dataR[(dataR["Reference_start"]-1)<s]["QNAME"]))
    sI=sL.intersection(sR)
    dataL=dataL[dataL["QNAME"].isin(sI)]
    dataR=dataR[dataR["QNAME"].isin(sI)]
    data=pd.merge(dataL,dataR,on="QNAME")
    return len(set(list(data[abs(data["Template_end_x"]-data["Template_start_y"])<5]["QNAME"])))

def overlapUF(sample,base,outDir,pos):
    dataU=pd.read_csv(os.path.abspath(outDir)+"/splicing/"+sample+"/u-"+base+"_"+str(pos)+".f.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    extractFlagBits(dataU)
    dataU=extractStartEnd(dataU)
    return len(dataU[(dataU["Reference_start"]<pos)&(pos-dataU["Reference_start"]>20)&(dataU["Reference_end"]>pos)&(dataU["Reference_end"]-pos>20)])

def overlapUL(sample,base,outDir,pos):
    dataU=pd.read_csv(os.path.abspath(outDir)+"/splicing/"+sample+"/u-"+base+"_"+str(pos)+".sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL'])
    extractFlagBits(dataU)
    dataU=extractStartEnd(dataU)
    return len(dataU[(dataU["Reference_start"]<pos)&(pos-dataU["Reference_start"]>20)&(dataU["Reference_end"]>pos)&(dataU["Reference_end"]-pos>20)])


# read .sam files and compare results
def confirm(args):
    # begin analyzing the sam alignment output
    outDir=args.out
    df=pd.DataFrame([])
    for sample in os.listdir(os.path.abspath(outDir)+"/splicing"):
        print("analyzing:", sample)
        globList=[[x.split("/")[-1].split(".")[0].split("-")[1].split("_")[0],x.split("/")[-1].split(".")[0].split("-")[1].split("_")[1]] for x in list(glob.glob(os.path.abspath(outDir)+"/splicing/"+sample+"/u*f.sam"))]
        data=pd.DataFrame({"spliceSite":[x[0] for x in globList],"pos":[int(x[1]) for x in globList]})
        data["sample"]=sample
        if len(data)>0:
            data["supportLR"]=data.apply(lambda row: overlapLR(sample,row["spliceSite"],outDir,row["pos"],2),axis=1)
            data["supportUF"]=data.apply(lambda row: overlapUF(sample,row["spliceSite"],outDir,row["pos"]),axis=1)
            data["supportUL"]=data.apply(lambda row: overlapUL(sample,row["spliceSite"],outDir,row["pos"]),axis=1)
            df=pd.concat([df,data])

    # df.dropna(axis=0,inplace=True)
    df=df[(df["supportLR"]>0)|(df["supportUF"]>0)|(df["supportUL"]>0)]
    df.to_csv(os.path.abspath(outDir)+"/splicing/splicePos.csv")

def main(args):
    for fileN in glob.glob(os.path.abspath(args.input)+"/*R1*.fastq.gz"):
        fullPath=os.path.abspath(fileN)
        fileName=fullPath.split('/')[-1]
        dirPath="/".join(fullPath.split('/')[:-1])
        baseName="_R1".join(fileName.split("_R1")[:-1])
        baseEnd=fileName.split(baseName)[-1].split(".fastq.gz")[0].split('R1')[-1]
        print("Working on: "+baseName)
        consensusPath=os.path.abspath(args.out)+"/consensusHIV/"+baseName+baseEnd+".fa"
        if os.path.exists(consensusPath):
            poss=confirmSplice(baseName,consensusPath,args,baseEnd)

    confirm(args)

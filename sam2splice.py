#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import sys
import subprocess
from pybedtools import BedTool
import re
import argparse

#first need to extract splice junctions
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

def se(row):
    cigar=row["CIGAR"]
    chars=re.findall(r"[\D']+", cigar)
    ints=[int(x) for x in re.findall(r"[\d']+",cigar)]
    readLen=0
    pre=0
    post=0
    n=0
    m_pre_tem=0
    m_pre_ref=0
    m_post_tem=0
    m_post_ref=0
    indexN=0
    di=0
    blockCount=1
    blockSizes=[0]
    tStarts=[0]
    if "N" in chars:
        indexN=chars.index("N")
        blockCount=len(cigar.split("N"))
    for i in range(len(chars)):
        if i==0 and chars[i] in "SH":
            pre=ints[i]
            readLen=readLen+ints[i]
#             tStarts[0]=tStarts[0]+ints[i]
        if i==len(chars)-1 and chars[i] in "SH":
            post=ints[i]
            readLen=readLen+ints[i]
        if chars[i]=="N":
            n=n+ints[i]
            tStarts.append(tStarts[-1]+blockSizes[-1]+ints[i])
            blockSizes.append(0)
        if i<indexN and chars[i]=="M":
            m_pre_tem=m_pre_tem+ints[i]
            m_pre_ref=m_pre_ref+ints[i]
            readLen=readLen+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
        if i>=indexN and chars[i]=="M":
            m_post_tem=m_post_tem+ints[i]
            m_post_ref=m_post_ref+ints[i]
            readLen=readLen+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
        if i<indexN and chars[i]=="D":
            m_pre_ref=m_pre_ref+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
            di=di+ints[i]
        if i>=indexN and chars[i]=="D":
            m_post_ref=m_post_ref+ints[i]
            blockSizes[-1]=blockSizes[-1]+ints[i]
            di=di+ints[i]
        if i<indexN and chars[i]=="I":
            readLen=readLen+ints[i]
            m_pre_tem=m_pre_tem+ints[i]
            di=di+ints[i]
        if i>=indexN and chars[i]=="I":
            readLen=readLen+ints[i]
            m_post_tem=m_post_tem+ints[i]
            di=di+ints[i]
    return pd.Series([pre,post,m_pre_ref,m_pre_tem,m_post_ref,m_post_tem,n,readLen,di,blockCount,",".join([str(x) for x in blockSizes]),",".join([str(x) for x in tStarts])])

def parseCIGAR(data):
    data["CIGAR"].replace("*",np.nan,inplace=True)
    data.dropna(axis=0,inplace=True)
    data.reset_index(drop=True,inplace=True)

#     data["READ_LEN"]=data.SEQ.str.len()
    data["CIGAR_POST"]=data.CIGAR.str.extract("[M]([0-9]+)[A-Z]$",expand=False).replace(np.nan,0).astype(int)
    data["END"]=data.READ_LEN-data.CIGAR_POST
    data["CIGAR_PRE"]=data.CIGAR.str.extract("^([0-9]+)[SH]",expand=False).replace(np.nan,0).astype(int)

    data16=data[data["reversedCurr"]==16].reset_index(drop=True)
    data0=data[data["reversedCurr"]==0].reset_index(drop=True)
    data16["Template_start"]=data16.READ_LEN-data16.END
    data16["Template_end"]=data16.READ_LEN-data16.CIGAR_PRE
    data0["Template_start"]=data0.CIGAR_PRE
    data0["Template_end"]=data0.END

    data16["Reference_start"]=data16.READ_LEN-data16.END+data16.POS-data16.Template_start
    data16["Reference_end"]=data16.READ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start+data16.N
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE+data0.N 
    
    data=pd.concat([data16,data0]).reset_index(drop=True)
    data.drop(["CIGAR_POST","CIGAR_PRE"],axis=1,inplace=True)
    return data

def tStarts(row):
    re=row["Reference_start"]
    tStarts=[]
    qStarts=row.qStarts.split(",")
    qs2=[int(x)-int(qStarts[0]) for x in qStarts]
    blockSizes=row.blockSizes.split(",")
    for i in range(row["blockCount"]):
        tStarts.append(re+int(qStarts[i]))
    return pd.Series([",".join([str(x) for x in tStarts])])

def getSS(row):
    ns=row["N"]
    bs=[int(x) for x in row["blockSizes"].split(",")]
    ts=[int(x) for x in row["tStarts"].split(",")]
    res=[]
    for i in range(len(bs)-1):
        res.append(str(ts[i]+bs[i]-1)+"-"+str(ts[i+1]-1))
    return ",".join(res)

def findBlocks(row): # returns the minimum block size for a given junction
    sjs=row["fullSS"].split(",") # list of spliceJunctions
    cursj=row["ss"]
    fbs=row["blockSizes"].split(",") # list of blocks
    fbIDX=sjs.index(cursj) # index of the first junction
    curBlocks=[int(x) for x in fbs[fbIDX:fbIDX+2]]
    return min(curBlocks)

# main part of the tool
def sam2splice(args):
    df=pd.read_csv(os.path.abspath(args.input),sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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
    df=df[df["CIGAR"].str.contains("N")].reset_index(drop=True)
    # df.to_csv(args.out+".spliced.sam",index=False,header=False,sep="\t")
    df[["QNAME","SEQ"]].to_csv(args.out+".spliced.fastq",sep="\n",index=False,header=False)
    if len(df)==0:
        return
    extractFlagBits(df)

    df["PRE"]=np.nan
    df["POST"]=np.nan
    df["MPRER"]=np.nan
    df["MPRET"]=np.nan
    df["MPOSTR"]=np.nan
    df["MPOSTT"]=np.nan
    df["N"]=np.nan
    df["blockCount"]=0
    df["blockSizes"]=""
    df["qStarts"]=""
    df["tStarts"]=""
    df[["PRE","POST","MPRER","MPRET","MPOSTR","MPOSTT","N","READ_LEN","DI","blockCount","blockSizes","qStarts"]]=pd.DataFrame(df.apply(lambda row: se(row),axis=1))
    df=parseCIGAR(df)
    df["tStarts"]=df.apply(lambda row: tStarts(row),axis=1)
    df.drop(["FLAG","QUAL","paired","aligned2Mates","unmappedCurr","unmappedMate","reversedMate","firstRead","lastRead","secondaryAlignment","noPassFilter","PCRdup","suppAl","qStarts"],axis=1,inplace=True)

    df["ss"]=df.apply(lambda row: getSS(row),axis=1)
    df.reset_index(drop=True,inplace=True)

    df["dup"]=df.duplicated("QNAME").astype(int)
    df["QNAME"]=df["QNAME"]+"_"+df["dup"].astype(str)

    tdf=pd.concat([pd.Series(row['QNAME'], row['ss'].split(',')) for _, row in df.iterrows()]).reset_index()
    tdf.columns=["ss","QNAME"]

    tdf=tdf.merge(df[["QNAME","MAPQ","SEQ","CIGAR","blockSizes","ss"]],on="QNAME",how="left")
    tdf.columns=["ss","QNAME","MAPQ","SEQ","CIGAR","blockSizes","fullSS"]

    seqdf=tdf[["QNAME","SEQ","ss"]]
    seqdf["QNAME"]=seqdf["QNAME"].str.split("_")[0]
    for sj in list(seqdf["ss"]):
        seqdf[seqdf["ss"]==sj][["QNAME","SEQ"]].to_csv(args.out+"_"+sj+".fasta",delimiter="\n",index=False,header=False)

    # create a copy and compute best reads by alignment quality
    seqDF=tdf.copy(deep=True)

    # first get all with highest mapq
    gdf=pd.DataFrame(seqDF.groupby(by="ss")["MAPQ"].max()).reset_index()
    gdf.columns=["ss","maxMAPQ"]
    tdf=tdf.merge(gdf,on="ss",how="left")
    tdf=tdf[tdf["MAPQ"]==tdf["maxMAPQ"]].reset_index(drop=True)

    tdf["minBlock"]=tdf.apply(lambda row: findBlocks(row),axis=1)
    tdf.sort_values(by="QNAME")

    maxdf=pd.DataFrame(tdf.groupby(by=["ss"])["minBlock"].max()).reset_index()
    maxdf.columns=["ss","mb"]
    maxdf["mb"]=maxdf["mb"].astype(str)
    tdf["minBlock"]=tdf["minBlock"].astype(str)
    tdf=tdf.merge(maxdf,on="ss",how="left")
    fdf=tdf[tdf["minBlock"]==tdf["mb"]].reset_index(drop=True)
    
    fdf.drop_duplicates(subset=["ss"],keep="first",inplace=True)
    fdf.reset_index(inplace=True,drop=True)
    
    gdf=tdf.groupby(by=["ss"])["QNAME"].count().reset_index()
    gdf=gdf.merge(fdf[["ss","SEQ"]],on="ss",how="left")

    donors={726:"D1b",
            743:"D1",
            747:"D1c",
            4721:"D1a",
            4962:"D2",
            5463:"D3",
            6045:"D4",
            6723:"D5",
            8433:"D6"}

    acceptors={4542:"A1a",
                4912:"A1",
                5389:"A2",
                5777:"A3",
                5936:"A4c",
                5954:"A4a",
                5960:"A4b",
                5976:"A5",
                5980:"A5a",
                5983:"A5b",
                6607:"A6",
                6649:"A6a",
                8354:"A7a",
                8350:"A7b",
                8354:"A7c",
                8372:"A7d",
                8378:"A7",
                8391:"A7e",
                8495:"A7f",
                9135:"A8d",
                9151:"A8a",
                9163:"A8c",
                9171:"A8",
                9184:"A8b",
                9190:"A8e",
                9228:"A8f"}

    gdf[["d","a"]]=gdf["ss"].str.split("-",expand=True)
    gdf["d"]=gdf.d.astype(int)
    gdf["a"]=gdf.a.astype(int)
    # now need to remove everything from the LTR region
    gdf.columns=["spliceJunction","numReads","bestSequence","donor","acceptor"]
    if args.filter:
        ltr3=[9086,9719]
        ltr5=[1,634]
        gdf=gdf[(gdf["donor"]>=ltr5[-1])&(gdf["acceptor"]<=ltr3[0])].reset_index(drop=True)
    gdf.sort_values(by="donor",inplace=True)
    gdf.reset_index(drop=True,inplace=True)
    gdf.replace({"donor": donors},inplace=True)
    gdf.replace({"acceptor": acceptors},inplace=True)

    gdf.to_csv(args.out,index=False)


def main(argv):
    parser=argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('-i',
                        '--input',
                        required=True,
                        type=str,
                        help="splice alignment")
    parser.add_argument('-o',
                        '--out',
                        required=True,
                        type=str,
                        help="output file")
    parser.add_argument('-f',
                        '--filter',
                        action="store_true",
                        help="filter LTR splice junctions")
    
    parser.set_defaults(func=sam2splice)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])

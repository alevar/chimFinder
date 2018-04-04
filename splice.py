#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import sys
import itertools
import math
import re

# ./splice.py <infp> <outfp>

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
    pre=0
    post=0
    n=0
    m_pre_tem=0
    m_pre_ref=0
    m_post_tem=0
    m_post_ref=0
    indexN=chars.index("N")
    for i in range(len(chars)):
        if i==0 and chars[i]=="S":
            pre=ints[i]
        if i==len(chars)-1 and chars[i]=="S":
            post=ints[i]
        if chars[i]=="N":
            n=ints[i]
        if i<indexN and chars[i]=="M":
            m_pre_tem=m_pre_tem+ints[i]
            m_pre_ref=m_pre_ref+ints[i]
        if i>indexN and chars[i]=="M":
            m_post_tem=m_post_tem+ints[i]
            m_post_ref=m_post_ref+ints[i]
        if i<indexN and chars[i]=="D":
            m_pre_ref=m_pre_ref+ints[i]
        if i>indexN and chars[i]=="D":
            m_post_ref=m_post_ref+ints[i]
        if i<indexN and chars[i]=="I":
            m_pre_tem=m_pre_tem+ints[i]
        if i>indexN and chars[i]=="I":
            m_post_tem=m_post_tem+ints[i]
    
    return pd.Series([pre,post,m_pre_ref,m_pre_tem,m_post_ref,m_post_tem,n])

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
    data16["Reference_end"]=data16.READ_LEN-data16.CIGAR_PRE-1+data16.POS-data16.Template_start+data16.N
    data0["Reference_start"]=data0.POS
    data0["Reference_end"]=data0.END+data0.POS-data0.CIGAR_PRE+data0.N

    data=pd.concat([data16,data0]).reset_index(drop=True)
    data.drop(["CIGAR_POST","END","CIGAR_PRE"],axis=1,inplace=True)
    return data
    
def main(argv):
    infp=os.path.abspath(argv[1])
    outfp=os.path.abspath(argv[2])
    df=pd.read_csv(infp,sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME',
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

    df=df[df["CIGAR"].str.contains("N")].sort_values(by="CIGAR").reset_index(drop=True)
    extractFlagBits(df)
    df["PRE"]=np.nan
    df["POST"]=np.nan
    df["MPRER"]=np.nan
    df["MPRET"]=np.nan
    df["MPOSTR"]=np.nan
    df["MPOSTT"]=np.nan
    df["N"]=np.nan
    df[["PRE","POST","MPRER","MPRET","MPOSTR","MPOSTT","N"]]=pd.DataFrame(df.apply(lambda row: se(row),axis=1))
    df=extractStartEnd(df)
    df["ss"]=df["Reference_start"]+df["MPRER"]
    df["se"]=df["ss"]+df["N"]
    df["spliceCoordinates"]=df["ss"].astype(str)+":"+df["se"].astype(str)
    df.drop(["ss","se"],axis=1,inplace=True)
    aggregations={
        'QNAME':{
            'count':'count',
            'reads': lambda x: set(x)
        },
        'SEQ':{
            'seq': lambda x: set(x)
        }
    }
    df=pd.DataFrame(df.groupby(by=["spliceCoordinates"])[["QNAME","SEQ","spliceCoordinates"]].agg(aggregations)).sort_values(by="count",ascending=False).reset_index()
    df["seq"]=df['seq'].str.join(";")
    df["reads"]=df['reads'].str.join(";")

    if not df is None and len(df)>0:
        df.to_csv(outfp)

if __name__=="__main__":
    main(sys.argv)
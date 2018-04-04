#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import sys
import itertools
import math
import re
import csv
import pysam

# ./splice.py <infp> <outfp>
    
def main(argv):
    infp=os.path.abspath(argv[1])
    outfp=os.path.abspath(argv[2])
    inSAM=pysam.AlignmentFile(infp,"rb")
    df=pd.DataFrame([])
    # res = collections.Counter()
    for r in inSAM.fetch():
        if 'N' in r.cigarstring:
            last_read_pos = False
            for read_loc, genome_loc in r.get_aligned_pairs():
                if read_loc is None and last_read_pos:
                    start = genome_loc
                elif read_loc and last_read_pos is None:
                    stop = genome_loc  # we are right exclusive ,so this is correct
                    df=pd.concat([df,pd.DataFrame([[str(start)+":"+str(stop),r.qname,r.query_sequence]])],axis=0)
                    del start
                    del stop
                last_read_pos = read_loc

    if not df is None and len(df)>0:
        df.columns=["spliceCoordinates","QNAME","SEQ"]
        df.reset_index(drop=True,inplace=True)
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
        df.to_csv(outfp)

if __name__=="__main__":
    main(sys.argv)
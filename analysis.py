#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
import sys
import glob
import time
import itertools
import scipy

def getReferenceLoc(row):
    if row["CIGAR"]=="*":
        return "-1"
    t=["".join(x) for _, x in itertools.groupby(row["CIGAR"], key=str.isdigit)]
    return str(row["POS"]+int(t[t.index("M")-1])-1)

def getTemplateLoc(row):
    t=["".join(x) for _, x in itertools.groupby(row["CIGAR"], key=str.isdigit)]
    start=0
    end=0
    if row["TLEN"]==0:
        return 0
    if row["TLEN"]<0:
        lenAl=sum(list(map(int,t[::2])))
        if "S" in t or "H" in t:
            start=start+int(t[t.index("S")-1])
            end=end+int(t[t.index("S")-1])
        end=end+int(t[t.index("M")-1])
        return str(start)+":"+str(end)
    if row["TLEN"]>0:
        if "S" in t:
            start=start+int(t[t.index("S")-1])
            end=end+int(t[t.index("S")-1])
        end=end+int(t[t.index("M")-1])
        return str(start)+":"+str(end)

def samFlag(flag):
	if flag & 1:
	   print('paired')
	# if not flag & 1:
	#    print('unpaired')

	if flag & 2:
	   print('read mapped in proper pair')
	# if not flag & 2:
	#    print('read not mapped in proper pair')

	if flag & 4:
	   print('read unmapped')
	# if not flag & 4:
	#    print('not read unmapped')

	if flag & 8:
	   print('mate unmapped')
	# if not flag & 8:
	#    print('not mate unmapped')

	if flag & 16:
	   print('read reverse strand')
	# if not flag & 16:
	#    print('not read reverse strand')

	if flag & 32:
	   print('mate reverse strand')
	# if not flag & 32:
	#    print('not mate reverse strand')

	if flag & 64:
	   print('first in pair')
	# if not flag & 64:
	#    print('not first in pair')

	if flag & 128:
	   print('second in pair')
	# if not flag & 128:
	#    print('not second in pair')

	if flag & 256:
	   print('not primary alignment')
	# if not flag & 256:
	#    print('primary alignment')

	if flag & 512:
	   print('read fails platform/vendor quality checks')
	# if not flag & 512:
	#    print('read passes platform/vendor quality checks')

	if flag & 1024:
	   print('read is PCR or optical duplicate')
	# if not flag & 1024:
	#    print('read is not PCR or optical duplicate')

	if flag & 2048:
	   print('supplementary alignment')
	# if not flag & 2048:
	#    print('not supplementary alignment')

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
	len(splitMates)

def main():
	dataHIV = pd.read_csv("./108PI_pos_S2.hiv.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21'])
	dataHUM = pd.read_csv("./108PI_pos_S2.hum.sam",sep="\t",comment='@',usecols=[0,1,2,3,4,5,6,7,8,9,10],names=['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','12','13','14','15','16','17','18','19','20','21'])
	dataHIV["paired"]=dataHIV["FLAG"]&1 #template having multiple segments in sequencing
	dataHIV["aligned2Mates"]=dataHIV["FLAG"]&2 #each segment properly aligned according to the aligner
	dataHIV["unmappedCurr"]=dataHIV["FLAG"]&4 #segment unmapped
	dataHIV["unmappedMate"]=dataHIV["FLAG"]&8 #next segment in the template unmapped
	dataHIV["reversedCurr"]=dataHIV["FLAG"]&16 #SEQ being reverse complemented
	dataHIV["reversedMate"]=dataHIV["FLAG"]&32 #SEQ of the next segment in the template being reverse complemented
	dataHIV["firstRead"]=dataHIV["FLAG"]&64 #the first segment in the template
	dataHIV["lastRead"]=dataHIV["FLAG"]&128 #the last segment in the template
	dataHIV["secondaryAlignment"]=dataHIV["FLAG"]&256 # secondary alignment
	dataHIV["noPassFilter"]=dataHIV["FLAG"]&512 #not passing filters, such as platform/vendor quality controls
	dataHIV["PCRdup"]=dataHIV["FLAG"]&1024 #PCR or optical duplicate
	dataHIV["suppAl"]=dataHIV["FLAG"]&2048 #supplementary alignment

	dataHUM["paired"]=dataHUM["FLAG"]&1 #template having multiple segments in sequencing
	dataHUM["aligned2Mates"]=dataHUM["FLAG"]&2 #each segment properly aligned according to the aligner
	dataHUM["unmappedCurr"]=dataHUM["FLAG"]&4 #segment unmapped
	dataHUM["unmappedMate"]=dataHUM["FLAG"]&8 #next segment in the template unmapped
	dataHUM["reversedCurr"]=dataHUM["FLAG"]&16 #SEQ being reverse complemented
	dataHUM["reversedMate"]=dataHUM["FLAG"]&32 #SEQ of the next segment in the template being reverse complemented
	dataHUM["firstRead"]=dataHUM["FLAG"]&64 #the first segment in the template
	dataHUM["lastRead"]=dataHUM["FLAG"]&128 #the last segment in the template
	dataHUM["secondaryAlignment"]=dataHUM["FLAG"]&256 # secondary alignment
	dataHUM["noPassFilter"]=dataHUM["FLAG"]&512 #not passing filters, such as platform/vendor quality controls
	dataHUM["PCRdup"]=dataHUM["FLAG"]&1024 #PCR or optical duplicate
	dataHUM["suppAl"]=dataHUM["FLAG"]&2048 #supplementary alignment

if __name__=="__main__":
	main()
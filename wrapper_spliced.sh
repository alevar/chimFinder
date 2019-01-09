#!/usr/bin/env bash

for sample in $(ls out*/alns/*.bam); do
    echo $(basename ${sample})
    sb=${sample%%.*}
    
    #extract spliced reads
    samtools view ./${sample} | awk -F '\t' '$6 ~/N/ {print}' > ./${sb}.s.tmp.sam
    
    #extract header
    samtools view -H ./${sample} > ./${sb}.h.tmp.sam
    # combine header with spliced reads
    cat ./${sb}.h.tmp.sam ./${sb}.s.tmp.sam > ./${sb}.s.sam

    # remove tmp header and alignment files
    rm ./${sb}.h.tmp.sam ./${sb}.s.tmp.sam

    # run custom pipeline
done;

#!/bin/sh
### Extract Gene Symbol and Gene ID Info from the gtf file

gtf="./Mus_musculus.GRCm39.107.gtf"


### gene_id to gene_name
grep 'gene_id' $gtf | awk -F 'gene_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'gene_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >g2s_mm39_gencode.txt
rm *_tmp

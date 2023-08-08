#!/bin/bash

source ~/.bashrc
conda activate <env_name>


filename='/data/gene_list/xaa'
while read line; do
# reading each line
python3 /data/scripts/enst_ref_seq.py $line
echo $line
done < $filename

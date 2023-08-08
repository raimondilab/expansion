#!/bin/bash

# -------------------------------------------------------------------------
# Author: Chakit Arora (PhD)
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa
# email: chakit.arora@sns.it
# -------------------------------------------------------------------------



# Activate the specified conda environment where CD-HIT is installed
source ~/.bashrc
conda activate <env_name>

# Set the input filename
filename='gene_list.txt'

# Read each line from the file 'gene_list.txt'
while read line; do

  # Use cd-hit to cluster sequences from the FASTA file
  # -d 70: maximum length difference for a sequence to be clustered
  # -c 0.6: sequence identity threshold for clustering
  # -t 3: number of threads to use
  # -i: input FASTA file
  # -o: output file for clustered sequences
  cd-hit -d 70 -c 0.6 -t 3 -i /data/fasta/${line}.fa -o /data/cd_hit/${line}.fa

  # Print the current gene symbol
  echo $line

done < $filename

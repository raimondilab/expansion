#!/bin/bash

# Activate the conda environment
source ~/.bashrc
conda activate <env_name>

# Set the input filename
filename='gene_list.txt'

# Read each line from the file 'gene_list.txt'
while read line; do

  # Use awk to process each line of the CSV file and create FASTA format entries
  # Field separator is set to '\t' for tab-separated values
  # The if condition checks if the first field does not match "symbol"
  # Then, it prints the formatted FASTA header and sequence
  awk -F '\t' '{if($1 !~ "symbol")printf(">%s_%s_%s\n%s\n", $1,$2,$3,$4)}' /data/res_transcripts/${line}_transcripts.csv | sed 's/"//g' > /data/fasta/${line}.fa

  # Print the current gene symbol
  echo $line

done < $filename

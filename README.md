![action main](https://github.com/raimondilab/precogx/actions/workflows/main.yml/badge.svg)

EXploring
Protein
AlterNative
SplIcing
cONsequence 

<div align="center">
<img src="https://github.com/raimondilab/expansion/blob/main/expansion_home.svg" alt="logo" width="400"/>
</div>
<br>

Welcome to the EXPANSION Complementary Repository!

Explore the functional implications of protein-coding alternative splice variants with [EXPANSION](https://expansion.bioinfolab.sns.it/), an integrated web-server designed to unravel the impact of alternative splicing in cancer genomics. This repository houses the data, scripts, and resources supporting the research paper. Combined with DE protein-coding transcripts, domain analysis, protein interactions, and gene enrichment, EXPANSION offers an intuitive glimpse into the effects of splice variants. Dive into Ensembl transcripts, Interpro domains, PTMs, and more, to uncover functionally significant splicing events. Analyze pre-calculated or custom DE transcript datasets effortlessly, accelerating your journey to gain insights into protein spliceforms.

## Workflow
We analyzed Ensembl transcripts  and their corresponding protein sequences. We clustered the protein sequences, aligned them, and compared them to Uniprot canonical sequences to identify differences in alternatively spliced isoforms.
<div align="center">
<img src="https://github.com/raimondilab/expansion/blob/main/expansion_workflow_about.svg" alt="logo" width="600"/>
</div>

## Structure of this Repository
There are two main folders: 

1. '/data' containing the datasets, necessary input and output files categorized in different subfolders.

2. '/scripts' containing the codes used in development of EXPANSION.


## This is the pipeline with a brief description of codes and data used:

**1. gene_list.txt:** the list of protein coding genes from [HGNC](https://www.genenames.org/download/statistics-and-files/) 

**2. tx_ensembldb.R:** R script intended to retrieve information about protein-coding transcripts for a given gene symbol using the Ensembl database and the AnnotationHub package. It then processes the retrieved data and saves it to a specified output file (/data/res_transcripts).

**3. fasta.sh:** The provided script is a Bash script designed to process a list of gene symbols from a text file (gene_list.txt), create FASTA format entries from associated data in CSV files (obtained from step-2), and save the resulting sequences in separate FASTA files (/data/fasta). The script utilizes awk as follows:
```
awk -F '\t' '{if($1 !~ "symbol")printf(">%s_%s_%s\n%s\n", $1,$2,$3,$4)}'  /data/res_transcripts/<gene_symbol>_transcripts.csv  | sed 's/"//g' > /data/fasta/<gene_symbol>.fa
```

**4. CD-HIT :** The provided script is a Bash script that utilizes the CD-HIT tool to cluster sequences from FASTA files (/data/fasta) based on sequence identity. It reads gene symbols from an input file (gene_list.txt), clusters sequences using CD-HIT, and outputs the clustered sequences into separate FASTA files (/data/cd_hit). It utilizes the following command in batch mode:
```
cd-hit -d 70 -c 0.6 -t 3 -i /data/fasta/<gene_symbol>.fa -o /data/cd_hit/<gene_symbol>.fa
```
The script uses the CD-HIT tool with specific parameters:
-d 70: Maximum length difference for a sequence to be clustered.
-c 0.6: Sequence identity threshold for clustering.
-t 3: Number of threads to use for processing.
-i: Specifies the input FASTA file path (/data/fasta/<gene_symbol>.fa).
-o: Specifies the output file path for the clustered sequences (/data/cd_hit/<gene_symbol>.fa).

**6. enst2cdhit.py:** This script essentially extracts ENST IDs and their corresponding CD-HIT cluster numbers from CD-HIT output files (/data/cd_hit), compiles them into a dictionary, and then saves the dictionary as a pickle file (enst2cdhit.pickle) for further use or analysis. Users can customize the script by modifying the input directory and filenames.
   
**7. enst_ref_seq.py:** The provided Python script processes a gene symbol by loading Ensembl transcripts, mapping them to UniProt ACs, grouping by cd-hit clusters, identifying the longest protein sequence within each cluster as the "reference" isoform, generating a fasta file (/data/res_fasta/<gene_symbol>.fa) with protein sequences in each cluster starting with the reference isoform, and finally conducting a multiple sequence alignment of these transcripts using ClustalO (/data/res_fasta/<gene_symbol>_ali.fa). 

The code uses the following pickles:

**ENST2UniAC_bk.pickle:** This pickle file maps Ensembl transcript IDs to UniProt ACs.
**enst2cdhit.pickle:** This pickle file maps Ensembl transcript IDs to cd-hit clusters.
**GN2ACmain_sprotfasta.pickle:** This pickle file maps gene symbols to canonical UniProt ACs.
**sprotfasta.pickle:** This pickle file maps gene symbols to canonical UniProt sequences.

The code also requires the following input directories:

**/data/res_transcripts:** This directory contains the file of Ensembl transcripts for the given gene symbol.
**/data/res_fasta:** This directory is where the fasta files and multiple sequence alignment files will be created.

The .hmm file is a Hidden Markov Model file that is used by ClustalO to perform the multiple sequence alignment. The HMM file contains information about the structure of the protein that is being aligned. This information can be used by ClustalO to improve the accuracy of the alignment.

**8. msa_diff_pos_annotate.py:** The code first loads two pickle files that contain information about post-translational modifications (PTMs) and binding regions. Then, it reads a file that contains a list of filenames (for eg /data/chunks/file0.txt) of multiple sequence alignments (MSAs). The code reads each MSA in turn and extracts information about deletions, insertions, divergent positions, Interpro domains, PTMs, and binding regions. From this information, it creates a table and saves it to a file. The filename of the output file is the same as the input file (<gene_symbol>_ali.fa), with the suffix _diff.txt. The output is saved in the directory '/data/msa_diff_fixed/'

**9. tx_ebseq.R:** This script conducts transcriptomic analysis through the EBSeq package in R, focusing on differential expression analysis between two conditions: GTEX and TCGA. Input gene expression data in CSV format is processed (/data/tx_count_data/), utilizing Ensembl database annotations for mapping transcript IDs to gene symbols. The outcomes encompassing differential expression results, fold changes, and annotation details are saved in separate CSV output files. For execution, the script requires a command-line argument specifying the cancer type of interest, and the gene expression data must be in the input directory, while the output files will be generated in the designated output directory. In the absence of a specified cancer type, the script defaults to performing an analysis example for the 'Colon' cancer type.



#### Contacts 
Francesco Raimondi: francesco.raimondi@sns.it<br>
Chakit Arora: chakit.arora@sns.it<br>
Natalia de Oliveira Rosa: natalia.deoliveirarosa@sns.it

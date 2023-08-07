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
<img src="https://github.com/raimondilab/expansion/blob/main/expansion_workflow_about.svg" alt="logo" width="400"/>
</div>

## The pipeline with a brief description of codes and data used:
1. gene_list : the list of protein coding genes from [HGNC](https://www.genenames.org/download/statistics-and-files/) 
2. tx_ensembldb.R: transcript to various annotation such as sequences, length, identifiers etc.
3. extracting .fa to fasta folder using the command (fasta.sh)
```
awk -F '\t' '{if($1 !~ "symbol")printf(">%s_%s_%s\n%s\n", $1,$2,$3,$4)}'  /home/carora/splice_pipeline_new3/res_transcripts/${line}_transcripts.csv  | sed 's/"//g' > /home/carora/splice_pipeline_new3/fasta/${line}.fa
```

4. CD-HIT :
```
cd-hit -d 70 -c 0.6 -t 3 -i /home/carora/splice_pipeline_new3/fasta/${line}.fa -o /home/carora/splice_pipeline_new3/cd_hit/${line}.fa
```
6. creating enst2cdhit pickle: ensembl ids to cd-hit cluster mapping
7. enst_ref_seq.py: all the pickles, clustalo, .hmm file
   performs multiple sequence alignment of sequences within a cdhit cluster. outputs _ali.fa files
8. msa_diff_pos_annotate.py
9. EBSEQ

#### Contacts 
Francesco Raimondi: francesco.raimondi@sns.it<br>
Chakit Arora: chakit.arora@sns.it<br>
Natalia de Oliveira Rosa: natalia.deoliveirarosa@sns.it

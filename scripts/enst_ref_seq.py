
# -------------------------------------------------------------------------
# Author: Chakit Arora (PhD), Francesco Raimondi
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa
# email: chakit.arora@sns.it

# USAGE: python3 enst_ref_seq.py <gene_symbol>
# Output: two MSA files in /data/res_fasta; one .fa file and one _ali.fa file
# to run in batch mode a script like splice1.sh can be used wherein the gene-list was split into
# multiple chunks for eg /data/gene_list/xaa
# -------------------------------------------------------------------------



import pickle
import pandas as pd
import os
import sys
import operator

# Read the input argument from command line (gene_symbol in caps)
inputarg = sys.argv[1]

# Load ENST to UniProt AC mapping
with open('/data/ENST2UniAC_bk.pickle', 'rb') as handle:
    enst2ac = pickle.load(handle)

# Load ENST to cd-hit cluster mapping
with open('/data/enst2cdhit.pickle', 'rb') as handle:
    enst2cdhit = pickle.load(handle)

# Load gene symbol to canonical UniProt AC mapping
with open('/data/GN2ACmain_sprotfasta.pickle', 'rb') as handle:
    gn2ac = pickle.load(handle)

# Load gene symbol to canonical UniProt sequence mapping
with open('/data/sprotfasta.pickle', 'rb') as handle:
    gn2seq = pickle.load(handle)

# Read the transcripts file generated from ensembldb
txs = pd.read_csv('/data/res_transcripts/' + inputarg + "_transcripts.csv", sep="\t")
txs = txs.fillna('-')

# Map ENST IDs to cd-hit clusters
txs['cdhit'] = txs['tx_id'].map(enst2cdhit).astype(str)

# Group transcripts by cd-hit clusters
gb = txs.groupby('cdhit')

for i in txs['cdhit'].unique():
    df = gb.get_group(i)
    df = df.fillna('-')

    # Map ENST to gene symbols
    df['uniprot'] = df['tx_id'].map(enst2ac).str.split('-').str[0]
    df['uniprot_splice_var'] = df['tx_id'].map(enst2ac)
    df = df.fillna('-')

    df.loc[df['uniprot'].isnull(), 'uniprot'] = '-'

    # Find reference isoform from gn2seq
    df1 = df[df['protein_sequence'].isin(gn2seq[inputarg])].sort_values(by=['Seq_length', 'uniprot_splice_var'], ascending=[False, False])

    if df1.shape[0] != 0:
        df.loc[df1.index[0], 'reference'] = "*"
        df.drop_duplicates(subset=['tx_id'], keep='first', inplace=True)
        df = df.fillna("")
        df = df.sort_values(by=['reference'], ascending=False).reset_index(drop=True)

        df = df[['symbol', 'tx_id', 'uniprot_splice_var', 'reference', 'protein_sequence', 'protein_domain_id', 'cdhit']]
        df["Seq_id"] = df['symbol'].astype(str) + '|' + df['tx_id'].astype(str) + '|' + df['uniprot_splice_var'] + '|' + df['reference']

        df.to_csv('/data/res_fasta/' + inputarg + '_' + str(df.iloc[0, 2]) + '_enst.txt', header=True, index=False, sep='\t')
        fastadict = dict(zip(df.Seq_id, df.protein_sequence))

        # Create fasta file
        file = open('/data/res_fasta/' + inputarg + '_' + str(df.iloc[0, 2]) + ".fa", "w")
        for key in fastadict:
            file.write(">" + key + "\n" + fastadict[key] + "\n")
        file.close()

        # Multiple sequence alignment with clustalo
        command = '/data/clustalo -i' + ' ' + '/data/res_fasta/' + inputarg + '_' + str(df.iloc[0, 2]) + '.fa' + ' ' + '-o' + ' ' + '/data/res_fasta/' + inputarg + '_' + str(df.iloc[0, 2]) + '_ali.fa' + ' --force' + ' ' + '--hmm-in=/data/7tm_1_2020.hmm'
        os.system(command)

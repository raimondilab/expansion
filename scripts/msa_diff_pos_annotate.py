
# -------------------------------------------------------------------------
# Author: Chakit Arora (PhD), Francesco Raimondi
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa
# email: chakit.arora@sns.it

# USAGE: python3 msa_diff_pos_annotate.py /data/chunks/file0.txt
# Output: a txt file in /data/msa_diff_fixed/ with MSA differences wrt canonical isoform in the
# cluster. The txt file also has PTM, domain and binding region information
# -------------------------------------------------------------------------





import pandas as pd
import os,sys,operator, math
import gzip, pickle
import numpy as np
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from scipy.stats import hypergeom


###Load PTMs from PSP - replace the pickle with a sqlite database to make the search faster
with open("/data/PSP.pickle", 'rb') as handle:
    ptms=pickle.load(handle)

####Load binding regions from IntAct
with open("/data/IntAct_BR.pickle", 'rb') as handle:
    bindregs=pickle.load(handle)

###Creating an InterproDB object for retrieval of domain info
#interpro=interprodb()

#_ali.fa filenames were split into different chunks file0.txt, file1.txt etc in /data/chunks/
chunk=sys.argv[1] #eg file0.txt

for inputarg in open(chunk, "rt"):
    inputarg=inputarg.strip("\n")
    inputarg=inputarg.replace('_ali.fa','')
    #try:
    aln = AlignIO.read('/data/res_fasta/'+inputarg+"_ali.fa", "fasta")
    sys.path.insert(0, "/data/scripts/")
    from sqldb_tools import interprodb
    interpro=interprodb() 
    #from sqldb_tools import interprodb
    for a in aln:
        unipid=a.id
        if unipid.find("*") != -1:
            ref_seq=a.seq
            ref_id=a.id
            uniac=unipid.split("|")[2]
#            if uniac.find("-1") != -1:
            uniac=uniac.split("-")[0]
            break
 

    try:
        interpro.extract([uniac])
    
        doms=interpro.make_unique()

        s1=doms[0]
        intpr={}
        for key, value in s1.items():
            intpr=value

        intpr1={}
        for key,value in intpr.items():
            for ii in range(int(value[1]),int(value[2])):
                intpr1[ii]=key+":"+value[0]
    except:
        intpr1=[]


    try:
        ptm=ptms[uniac]
    except:
        ptm=''
    try:
        brs=bindregs[uniac]
    except:
        brs=[]
    track={}
    for a in aln:
        unipid=a.id
        if unipid.find("*") == -1:
                ###Sequences imported from the MSA are linear string. This makes the parsing simpler
                #print (ref_seq)
                #print (a.seq)
            track[a.id]={"Insertions":[], "Deletions":[], "Divergent_positions":[],"InterproDomains":[], "PTMs":[], "BindRegs":[]}
                ###Looping now through the sequences to spot differences
                ###Using two counters, one for the canonical (ref)sequence and one for the target (tar) sequence
            ref_pos=0
            tar_pos=0
            ii=0
            flag=0
            for ref_aa in ref_seq:
                tar_aa=a.seq[ii]
                if ref_aa != "-":
                    ref_pos += 1
                if tar_aa != "-":
                    tar_pos += 1
                if ref_aa != "-" and tar_aa == "-":
                      #track[a.id]["del"].append(ref_pos)
                    track[a.id]["Deletions"].append(ii+1)
                    flag=1
                if ref_aa == "-" and tar_aa != "-":
                    #track[a.id]["ins"].append(ref_pos)
                    track[a.id]["Insertions"].append(ii+1)
                    flag=1
                if ref_aa != "-" and tar_aa != "-" and ref_aa != tar_aa:
                      #track[a.id]["div"].append(ref_pos)
                    track[a.id]["Divergent_positions"].append(ii+1)
                    flag=1
                if ref_pos in intpr1 and flag != 0:
                    affected_dom=intpr1[ref_pos]
                    if affected_dom not in track[a.id]["InterproDomains"]:
                        track[a.id]["InterproDomains"].append(affected_dom)
                if str(ref_pos) in ptm and flag != 0:
                    ptm_pos=str(ref_pos)+":"+ptm[str(ref_pos)]
                    if ptm_pos not in track[a.id]["PTMs"]:
                        track[a.id]["PTMs"].append(ptm_pos)
                if int(ref_pos) in brs and flag != 0:
                    brs_pos=str(ref_pos)+":"+str(brs[int(ref_pos)])
                    if brs_pos not in track[a.id]["BindRegs"]:
                        track[a.id]["BindRegs"].append(brs_pos)
                flag=0
                ii+=1


    ###Printing out the positions for each difference class, i.e. deletion, insertion, divergence
    diff=pd.DataFrame.from_dict(track, orient='index')
    diff= diff.rename_axis('Isoforms_id').reset_index()
    diff['Insertions'] = diff['Insertions'].astype(str).str.replace(r'\[]', 'nan')
    diff['Deletions'] = diff['Deletions'].astype(str).str.replace(r'\[]', 'nan')
    diff['Divergent_positions'] = diff['Divergent_positions'].astype(str).str.replace(r'\[]', 'nan')
    diff['InterproDomains'] = diff['InterproDomains'].astype(str).str.replace(r'\[]', 'nan')
    #print (diff['PTMs'])
    diff['PTMs'] = diff['PTMs'].astype(str).str.replace(r'\[]', 'nan')
    diff['BindRegs'] = diff['BindRegs'].astype(str).str.replace(r'\[]', 'nan')
    diff.to_csv('/data/msa_diff_fixed/'+inputarg+'_diff.txt', sep='\t', index=False)



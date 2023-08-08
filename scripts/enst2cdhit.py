

# -------------------------------------------------------------------------
# Author: Chakit Arora (PhD)
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa
# email: chakit.arora@sns.it
# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
# Usage: python enst2cdhit.py
# Input: The code expects a set of .clstr files located in the directory /data/cd_hit/
# Output: The code generates an output file named enst2cdhit.pickle. This output file contains a
# where keys are Ensembl transcript IDs (ENST) and values are cluster indices (clust)

# this code is used to process output files from the cd-hit tool, extract ENST IDs and their
# corresponding cd-hit clusters, and store this mapping in a dictionary that is then saved to a
# pickle file for future use.
# -------------------------------------------------------------------------



import os,sys,operator
import gzip, glob
import pickle
import glob

# Get a list of all .clstr files in the specified directory
f = glob.glob("/data/cd_hit/*.clstr")

# Create an empty dictionary to store the mapping of ENST IDs to cd-hit clusters
enst2cdhit_f = {}

# Iterate over each .clstr file
for file in f:
    k = -1
    enst2cdhit = {}
    
    # Open the .clstr file for reading line by line
    for line in open(file, "rt"):
        # Find lines that indicate the start of a new cluster
        if line.find("Cluster") != -1:
            k = k + 1  # Increment the cluster counter
            
        # Process lines containing ENST IDs
        if line.find("Cluster") != 1:
            enst = line.split('_')[1]  # Extract the ENST ID from the line
            clust = k  # Assign the current cluster to the ENST ID
            
            # Add the ENST ID and its corresponding cluster to the dictionary
            if enst not in enst2cdhit:
                enst2cdhit[enst] = clust
                
    # Update the main dictionary with the ENST to cluster mapping from the current file
    enst2cdhit_f.update(enst2cdhit)

# Save the enst2cdhit dictionary to a pickle file
with open('/data/enst2cdhit.pickle', 'wb') as handle:
    pickle.dump(enst2cdhit_f, handle, protocol=pickle.HIGHEST_PROTOCOL)

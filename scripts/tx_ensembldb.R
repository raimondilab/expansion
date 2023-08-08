#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Author: Chakit Arora, PhD
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa
# Usage : Rscript tx_ensembldb.R <gene_symbol> <output_path_filename>
# output folder in the repository: data/res_transcripts
# -------------------------------------------------------------------------

`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

# Load required libraries
library('ensembldb')
library('AnnotationHub')
# AnnotationHub::setAnnotationHubOption("PROXY", Sys.getenv("https_proxy"))

# Get command line arguments for gene symbol and output path
args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]
outputPath <- args[2]

# Create AnnotationHub object
hub <- AnnotationHub()

# Query EnsDb annotations from AnnotationHub
ahDb <- query(hub, pattern = c("Homo Sapiens", "EnsDb", 105))
ahEdb <- ahDb[[1]]

# Retrieve protein-coding transcripts for the gene
txs <- transcripts(ahEdb, filter = GeneNameFilter(c(gene)),
                   columns = c("symbol", "tx_id", "uniprot_id", "protein_sequence", "tx_biotype"),
                   return.type = "DataFrame")

# Filter for protein-coding transcripts
txs <- txs[txs$tx_biotype == "protein_coding", c("symbol", "tx_id", "uniprot_id", "protein_sequence")]

# Calculate sequence length
txs$Seq_length <- nchar(txs$protein_sequence)

# Remove version number from uniprot_id
txs$uniprot_id <- gsub("\\..*", "", txs$uniprot_id)

# Order by tx_id and remove duplicates
txs <- txs[order(txs[, 'tx_id']), ]
txs <- txs[!duplicated(txs$tx_id), ]

# Write data to file
write.table(txs, file = outputPath, sep = "\t", row.names = FALSE, append = FALSE)

# Print success message
cat("Transcripts saved to", outputPath, "\n")

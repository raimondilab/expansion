#!/usr/bin/env Rscript

# Author: Chakit Arora, PhD
# Affiliation: Bioinformatics Group, Scuola Normale Superiore in Pisa

# Script Description:
# This script performs transcriptomic analysis using the EBSeq package in R. It takes gene expression data in the form
# of a CSV file and performs differential expression analysis between two conditions (GTEX and TCGA). The script utilizes 
# Ensembl database annotations for mapping transcript IDs to gene symbols. The results of the analysis, including 
# differential expression results, fold changes, and annotation information, are stored in separate output CSV files.

# Usage:
# This script can be executed from the command line using Rscript. It requires one command-line argument, which is the
# name of the cancer type for which the analysis is to be performed. The gene expression data for the specified cancer 
# type should be present in a CSV file located in the input directory (/data/tx_count_data). The output files will be generated in the output
# directory. If no cancer type is provided, the script will run an example analysis for the 'Colon' cancer type.

# Example Usage:
# Run the script for the 'Colon' cancer type:
# $ Rscript transcript_analysis.R

# Run the script for a specific cancer type:
# $ Rscript transcript_analysis.R Breast

# Dependencies:
# This script requires the following R packages: dplyr, EnsDb.Hsapiens.v86, EBSeq.

# ------------------------
# Script starts here
# Define custom function for concatenating elements
`+` <- function(e1, e2) {
  if (is.character(e1) | is.character(e2)) {
    paste0(e1, e2)
  } else {
    base::`+`(e1, e2)
  }
}

# Set the name of the cancer to be entered while calling the script
args = commandArgs(trailingOnly = TRUE)
out_path = getwd()

# Set default cancer type to 'Colon' if no cancer type provided
if (length(args) == 0) {
  print('Running example for Colon')
  args[1] = 'Colon'
}

# Define input and output paths
input_path <- "/data/tx_count_data/"
output_path <- "/data/transcript_DE/"

# Read the data from the CSV file
data <- read.csv(paste0(input_path, args[1], "_transcript.tsv"), sep = ",")

library(dplyr)

# Remove normal samples from TCGA by excluding barcodes containing .11 at the end
data <- data %>% select(-contains(".11"))

data2 <- data[-1]
row.names(data2) <- data$sample
data2 <- data2[rowSums(data2[]) > 0, ]

data2 <- as.data.frame(t(as.matrix(data2)))

# Sort the rows by name
new_df <- data2[order(row.names(data2)), ]

# Create conditions based on the sample names
s <- c(row.names(new_df))
conditions <- rep(c('GTEX', 'TCGA'), times = c(length(grep('GTEX', s)), nrow(new_df) - length(grep('GTEX', s))))


suppressMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86

# Get all transcripts defined in Ensembl (version 86)
tx <- transcripts(edb, columns = c("tx_id", "gene_id", "gene_name"))

# Extract the transcript ids and gene names
mapping <- cbind(tx_id = tx$tx_id, name = tx$gene_name)
rownames(mapping) <- mapping[, 1]
mapping <- as.data.frame(mapping)

# Remove trailing characters after '.' in df
colnames(new_df) <- gsub("\\..*", "", colnames(new_df))

df <- as.data.frame(matrix(nrow = ncol(new_df)), ncol = 1)
colnames(df) <- 'transcript'
df$transcript <- colnames(new_df)

colnames(mapping) <- c('transcript', 'gene_symbol')

# Use transcript id to merge the two dataframes
res <- merge(df, mapping, by = "transcript", quiet = TRUE, all.x = TRUE)

new_df <- as.data.frame(t(as.matrix(new_df)))
colnames(new_df) <- NULL

### Running EBseq
library(EBSeq)

a <- which(is.na(res$gene_symbol))

IsoNames <- res$transcript[-a]
IsosGeneNames <- res$gene_symbol[-a]

new_df <- new_df[-a, ]
f <- as.matrix(new_df)
mode(f) <- "numeric"

IsoSizes <- MedianNorm(f)

NgList <- GetNg(IsoNames, IsosGeneNames)
IsoNgTrun <- NgList$IsoformNgTrun

# Options(error = recover)
IsoEBOut <- EBTest(Data = f, NgVector = IsoNgTrun, Conditions = as.factor(conditions), sizeFactors = IsoSizes, maxround = 10, Qtrm = 1, QtrmCut = 0)
IsoEBDERes <- GetDEResults(IsoEBOut, FDR = 0.05)

res_A <- IsoEBDERes$DEfound
res_B <- IsoEBDERes$PPMat
res_C <- IsoEBDERes$Status

GeneFC <- PostFC(IsoEBOut)  # logFC values

# Create output directory if it doesn't exist
dir.create(paste0(output_path, args[1]), showWarnings = FALSE)

# Write output files
write.csv(res_A, file = paste0(output_path, args[1], "/", args[1], "_ebseq_DEfound.csv"))
write.csv(res_B, file = paste0(output_path, args[1], "/", args[1], "_ebseq_PPmat.csv"))
write.csv(res_C, file = paste0(output_path, args[1], "/", args[1], "_ebseq_status.csv"))
write.csv(GeneFC, file = paste0(output_path, args[1], "/", args[1], "_ebseq_FC.csv"))

df1 <- read.csv(file = paste0(output_path, args[1], "/", args[1], "_ebseq_PPmat.csv"), sep = ",")
df2 <- read.csv(file = paste0(output_path, args[1], "/", args[1], "_ebseq_FC.csv"), sep = ",")

final <- merge(df1, df2, by = "X")
final$logFC <- log2(final$PostFC)
final <- final[order(-final$PPDE), ]

colnames(final)[1] <- 'transcript'
final <- merge(final, mapping, by = "transcript", quiet = TRUE, all.x = TRUE)

write.csv(final, file = paste0(output_path, args[1], "/", args[1], '_transcript_DE_annotated.csv'))

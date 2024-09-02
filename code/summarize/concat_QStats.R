# This script concatenates all the overlap statistics from a given dataset

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_OverlapStats.R <output file> <input files>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))


outfile = args[1]

dfOut <- matrix(NA, nrow = 1, ncol = 12)

for (i in 2:length(args)) {

  print(i)

  # Get results file name
  filename = args[i]

  # Extract dataset
  dataset <- strsplit(filename, "/")[[1]][3]

  # Extract which GWAS
  gwas <- strsplit(filename, "/")[[1]][4]

  # Extract Constrasts
  constrasts <- strsplit(filename, "/")[[1]][5]

  # Extract nsnp
  tmp <- strsplit(filename, "/")[[1]][6]
  L <- strsplit(tmp, "L-")[[1]][2]

  # Extract phenotype
  tmp <- strsplit(filename, "/")[[1]][7]
  phenotype <- strsplit(tmp, "_")[[1]][1]

  # Extract covar
  tmp <- strsplit(filename, "/")[[1]][9]
  covar <- strsplit(tmp, "_")[[1]][1]

  #Extract PC num
  tmp <- strsplit(filename, "/")[[1]][9]
  tmp2 <- strsplit(tmp, "PC")[[1]][2]
  PC <- strsplit(tmp2, ".results")[[1]][1]

  # Read in results
  df <- fread(filename)
  names_from_file <- colnames(df)
  df$dataset <- dataset
  df$gwas <- gwas
  df$contrasts <- constrasts
  df$L <- L
  df$phenotype <- phenotype
  df$covar <- covar
  df$pc <- PC

  colnames(dfOut) <- c(names_from_file, "dataset", "gwas", "contrasts", "L", "phenotype", "covar", "pc")
  dfOut <- rbind(dfOut, df)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
#colnames(dfOut) <- c(name_from_file, "dataset", "gwas", "contrasts")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")





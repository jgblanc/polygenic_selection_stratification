# This script concatenates results from selection scan into plotting format

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript concat_r2.R <output file> <input files>")}

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

  # Extract gwas
  gwas <- strsplit(filename, "/")[[1]][5]

  # Extract NSNP
  nsnp <- strsplit(filename, "/")[[1]][6]

  # Extract gwas type
  tmp <- strsplit(strsplit(filename, "/")[[1]][7], "_")[[1]][1]
  gwas_type <- strsplit(tmp, "-")[[1]][2]

  # Extract test type
  tmp <- strsplit(strsplit(filename, "/")[[1]][7], "_")[[1]][2]
  test_type <- strsplit(tmp, "-")[[1]][2]

  # Extract contrasts
  tmp <- strsplit(strsplit(filename, "/")[[1]][7], "_")[[1]][3]
  contrasts <- strsplit(tmp, ".txt")[[1]][1]

  # Read in results
  df <- fread(filename)[,1:6]
  names_from_file <- colnames(df)
  colnames(dfOut) <- c(names_from_file,"dataset", "gwas", "contrasts", "gwas_type", "test_type", "L")
  df$dataset <- dataset
  df$gwas <- gwas
  df$contrasts <- contrasts
  df$gwas_type <- gwas_type
  df$test_type <- test_type
  df$L <- nsnp
  dfOut <- rbind(dfOut, df)

}

# Remove first row
dfOut <- as.data.frame(dfOut[2:nrow(dfOut),])
colnames(dfOut) <- c(names_from_file,"dataset", "gwas", "contrasts", "gwas_type", "test_type", "L")

# Save file
print(dfOut)
fwrite(dfOut,outfile, row.names = F, col.names = T, quote = F, sep = "\t")





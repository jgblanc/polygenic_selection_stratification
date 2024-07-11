## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript standardized_r_IDs.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

inFile = args[1]
outFile = args[2]


# Initialize an empty list to store results
dfOut <- as.data.frame(matrix(NA,nrow = 1, ncol = 2))
colnames(dfOut) <- c("ID", "Var")
con <- file(inFile, open = "r")

# Read and process each line
while(TRUE) {

  line <- readLines(con, n = 1, warn = FALSE)

  if (length(line) == 0) break  # Exit loop if end of file

  # Gell fields
  fields <- strsplit(line, "\\s+")[[1]]
  dosages <- as.numeric(fields[7:length(fields)])
  ID <- fields[2]
  print(ID)


  # Convert to dosage of ALT allele
  dosages <- 2 - dosages

  # Calculate variance
  variance <- var(dosages)

  # Save results
  tmp <- c(ID, variance)
  dfOut <- rbind(dfOut, tmp)

}

# Save output
dfOut <- dfOut[3:nrow(dfOut), ]
fwrite(dfOut, outFile, row.names = F, col.names = T, quote = F, sep = "\t")


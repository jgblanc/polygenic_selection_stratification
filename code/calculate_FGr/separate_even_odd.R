# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

inFile = args[1]
outFile = args[2]
type = args[3]

# Read in SNPs
dfSNPs <- fread(inFile, header=FALSE)
print(head(dfSNPs))

if (type == "chr-all") {
   # Save all SNPs
   fwrite(dfSNPs, outFile, row.names = F, col.names = F, quote = F, sep = "\t")
   print(nrow(dfSNPs))
} else if (type == "chr-even") {
  # Pick even snps
  dfSNPs <- dfSNPs %>% separate(V1, into=c("chr", "pos"), sep = ":",  remove = FALSE)
  dfSNPs <- dfSNPs %>% filter(chr %in% c("2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22")) %>% select(V1)
  print(nrow(dfSNPs))
  fwrite(dfSNPs, outFile, row.names = F, col.names = F, quote = F, sep = "\t")
}






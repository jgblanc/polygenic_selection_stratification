# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

inFile = args[1]
outEven = args[2]
outAll = args[3]

# Read in SNPs
dfSNPs <- fread(inFile)

# Save all SNPs
fwrite(dfSNPs, outAll, row.names = F, col.names = F, quote = F, sep = "\t")
print(nrow(dfSNPs))

# Pick even snps
dfSNPs <- dfSNPs %>% separate(V1, into=c("chr", "pos"), sep = ":", keep = TRUE)
dfSNPs <- dfSNPs %>% filter(chr %in% c("2", "4", "6", "8", "10", "12", "14", "16", "18", "20", "22")) %>% select(V1)
print(nrow(dfSNPs))
fwrite(dfSNPs, outAll, row.names = F, col.names = F, quote = F, sep = "\t")







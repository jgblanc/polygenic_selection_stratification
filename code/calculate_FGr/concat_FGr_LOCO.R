## Concatenate projections
## This script reads in all the projected test vectors for each chromosome and adds them all together, leaving the focal chromosome out, and scales the final covariate

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript conccat_FGr_LOCO.R <prefix to Tm chromosomes> <chromosome number> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
snp_prefix = args[2]
chrNum = as.numeric(args[3])
outfile = args[4]
print(chrNum)

# Get list of all but the focal chromosome
chrs <- seq(1,22)
chrs <- chrs[ !chrs == chrNum]

# Find total number of SNPs
snp_nums <- fread(paste0(snp_prefix, "_", chrs[1], "_SNPcount.txt"))
for (i in 2:length(chrs)) {

  print(chrs[i])

  print(paste0("chr num ",chrs[i]))
  # Read in new chromosome
  snp_nums <- rbind(snp_nums, fread(paste0(snp_prefix,"_", chrs[i], "_SNPcount.txt")))

}
print(dim(snp_nums))
print(snp_nums)

# Add Tm for each chromosome to each other
data <- fread(paste0(Tm_prefix, "_", chrs[1], ".txt"))
data <- data[,4:ncol(data)]

for (i in 2:length(chrs)) {

  print(paste0("chr num ",chrs[i]))
  # Read in new chromosome
  filename <- paste0(Tm_prefix,"_", chrs[i], ".txt")
  tmp <- fread(filename)
  tmp <- tmp[,4:ncol(tmp)]
  data <- cbind(data, tmp)

}
data <- as.data.frame(data)
print(dim(data))

# Set parameters
L <- sum(snp_nums$nSNP)
M <- nrow(data)
print(paste0("L is ", L))
print(paste0("M is ", M))

# Compute FGrhat
FGr_hat <- apply(data, 1, sum) * (1/L)

# Construct output
dfOut <- fread(paste0(Tm_prefix,"_", chrs[1], ".txt"))
dfOut <- dfOut[,1:3]
dfOut$FGr <- scale(FGr_hat)

# Save output
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








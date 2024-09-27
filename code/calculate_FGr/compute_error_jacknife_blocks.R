## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
outfile = args[2]
snp_prefix = args[3]

# Get list of all chromosomes
chrs <- seq(1,22)

# Add Tm for each chromosome to each other
data <- fread(paste0(Tm_prefix, "_1.txt"))
data <- data[,4:ncol(data)]

for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  filename <- paste0(Tm_prefix,"_", i, ".txt")
  tmp <- fread(filename)
  tmp <- tmp[,4:ncol(tmp)]
  data <- cbind(data, tmp)

}
data <- as.data.frame(data)

# Find total number of SNPs
snp_nums <- fread(paste0(snp_prefix, "_1_SNPcount.txt"))
for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  snp_nums <- rbind(snp_nums, fread(paste0(snp_prefix,"_", i, "_SNPcount.txt")))

}


# Set parameters
L <- sum(snp_nums$nSNP)
M <- nrow(data)
print(paste0("L is ", L))
print(paste0("M is ", M))

# Compute D
FGr_hat <- (1/sqrt(L-1)) * apply(data, 1, sum)
print(var(FGr_hat))
D <- sum(FGr_hat^2)  * (1/M) * (1/(L-1))
print(D)

# Expected D
expD <- 1/(L-1)

# Compute SE for D
nblocks <- ncol(data)
allDs <- allDs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  mi <- as.numeric(snp_nums[i, 2])
  FGri <- data[,i] * (1/sqrt(mi-1))
  Di <- (t(FGri) %*% FGri) * (1/M) * (1 / (mi -1))
  allDs[i] <- (mi / (L - mi)) * (D - Di)^2

}
varD <- mean(allDs)
se <- sqrt(varD)

# Test D for significance
pval <- pnorm( D ,mean =expD, sd = se, lower.tail = FALSE)
print(pval)

# Compute SE for entries of FGr - Mair Notes
allFGrs <- matrix(NA, nrow = nrow(data), ncol = nblocks)
FGr_hat <- apply(data, 1, sum) * (1/L)
for (i in 1:nblocks) {

  mi <- as.numeric(snp_nums[i, 2])
  FGri <- data[,i] * (1/mi)
  allFGrs[,i] <- (mi / (L - mi)) * (FGri - FGr_hat)^2

}
allSigma2 <- rowMeans(allFGrs)
jkVar <- mean(allSigma2)

# Find Error
varFGr <- var(FGr_hat, na.rm = TRUE)
error <- jkVar / varFGr

# Final signal
signal <- 1 - error


# Make output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 9))
colnames(dfOut) <- c("H", "varH", "pvalH","M", "L", "jkFGr", "varFGr", "error", "signal")
dfOut[1,] <- c(D, varD, pval, M, L, jkVar, varFGr, error, signal)
print(dfOut)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








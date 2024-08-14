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
print(dim(data))

# Find total number of SNPs
snp_nums <- fread(paste0(snp_prefix, "_1_SNPcount.txt"))
for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  snp_nums <- rbind(snp_nums, fread(paste0(snp_prefix,"_", i, "_SNPcount.txt")))

}
print(dim(snp_nums))
print(snp_nums)

# Set parameters
L <- sum(snp_nums$nSNP)
M <- nrow(data)
print(paste0("L is ", L))
print(paste0("M is ", M))

# Compute D
FGr_hat <- apply(data, 1, sum) * (1/L)
D <- t(FGr_hat) %*% FGr_hat * (L^2)

# Expected D
expD <- (M-1)

# Compute SE for D
nblocks <- ncol(data)
allDs <- allDs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  mi <- as.numeric(snp_nums[i, 2])
  FGri <- data[,i] * (1/mi)
  Di <- t(FGri) %*% FGri * (mi^2)
  allDs[i] <- (mi / (L - mi)) * (D - Di)^2

}
varD <- mean(allDs)
se <- sqrt(varD)

# Test D for significance
pval <- pnorm(abs(D -expD) ,mean =0, sd = se, lower.tail = FALSE)


# Compute SE for entries of FGr - Mair Notes
allFGrs <- matrix(NA, nrow = nrow(data), ncol = nblocks)
for (i in 1:nblocks) {

  mi <- as.numeric(snp_nums[i, 2])
  FGri <- data[,i] * (1/mi)
  allFGrs[,i] <- (mi / (L - mi)) * (FGri - FGr_hat)^2

}
allSigma2 <- rowMeans(allFGrs)
jkVar <- mean(allSigma2)

# Find Error
varFGr <- var(FGr_hat, na.rm = TRUE, )
error <- jkVar / varFGr

# Comput SE for entries of FGr - Patterson notes

## Get jackknife estimate of FGr
weighted_FGrsLOCO <- matrix(NA, nrow = nrow(data), ncol = nblocks)
for (i in 1:nblocks) {
  print(i)
  mi <- as.numeric(snp_nums[i, 2])
  weighted_FGrsLOCO[,i] <- ((L - mi) / L) * (apply(data[,-i], 1, sum) * (1/(L- mi)))

}
thetaJ <- (nblocks * FGr_hat) - rowSums(weighted_FGrsLOCO)

## Get jackknife estimate of se
tau <- matrix(NA, nrow = nrow(data), ncol = nblocks)
for (i in 1:nblocks) {
  print(i)
  mi <- as.numeric(snp_nums[i, 2])
  hi <- (L / mi)
  tau[,i] <- (hi * FGr_hat) - ((hi - 1) * (apply(data[,-i], 1, sum) * (1/(L-mi))))
}
tmp <- matrix(NA, nrow = nrow(data), ncol = nblocks)
for (i in 1:nblocks) {
  print(i)
  mi <- as.numeric(snp_nums[i, 2])
  hi <- (L / mi)
  tmp[,i] <- (tau[,i] - thetaJ)^2  / (hi - 1)
}
allSigma2 <- rowMeans(tmp)
jkVar <- mean(allSigma2)

# Find Error
varFGr <- var(thetaJ)
error <- jkVar / varFGr

# Final signal
signal <- 1 - error

# Find gamma F
gammaF <- 1 - (error^2)

# Find Z
Z <- D / ((M-1) * L)

# Get traceK
trK <- (M-1) * L

# Make output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 11 ))
colnames(dfOut) <- c("D","ExpD", "varD", "pvalD","Z", "trK", "jkFGr", "varFGr", "error", "signal", "gammaF")
dfOut[1,] <- c(D, expD, varD, pval, Z, trK, jkVar, varFGr, error, signal, gammaF)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








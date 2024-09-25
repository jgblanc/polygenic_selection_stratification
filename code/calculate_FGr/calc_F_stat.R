## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyverse)
}))

betas_prefix = args[1]
var_prefix = args[2]
outfile = args[3]


# Get list of all chromosomes
chrs <- seq(1,22)

# Add betas for each chromosome to each other
data <- fread(paste0(betas_prefix, "_1.FGr.glm.linear"))
data <- data[,c(1,2,3,12)]

for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  filename <- paste0(betas_prefix,"_", i, ".FGr.glm.linear")
  tmp <- fread(filename)
  tmp <- tmp[, c(1,2,3,12)]
  data <- rbind(data, tmp)
}
print(dim(data))
colnames(data) <- c("CHR", "POS", "ID", "BETA")

# Get variance for all chromosomes
dfVar <- fread(paste0(var_prefix, "_1.txt"))

for (i in 2:22) {

  print(paste0("chr num ",i))
  # Read in new chromosome
  filename <- paste0(var_prefix,"_", i, ".txt")
  tmp <- fread(filename)
  dfVar <- rbind(dfVar, tmp)
}
print(dim(dfVar))


# Commbine Data
df <- inner_join(data,dfVar)
print(head(df))

# Calculate H
df$BETA <- df$BETA * (sqrt(df$Var))
H <- mean(df$BETA^2)
print(H)


# Calculate block jackknife for H
nblocks <- length(unique(df$block))
L <- nrow(df)
allHs <- rep(NA, nblocks)
for (i in 1:nblocks) {

  tmp <- df %>% filter(block == i)
  mi <- nrow(tmp)
  Hi <- mean(tmp$BETA^2)
  print(Hi)
  allHs[i] <- (mi / (L- mi)) * (H - Hi)^2

}
print(allHs)
varH <- mean(allHs)
se <- sqrt(varH)
expH <- 1/(100000)

# Test D for significance
pval <- pnorm( H ,mean =expH, sd = se, lower.tail = FALSE)

dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 4))
dfOut[1,] <- c(H, varH, se, pval)
colnames(dfOut) <- c("H", "varH", "seH", "pval")
print(dfOut)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








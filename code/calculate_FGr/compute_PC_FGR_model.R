## Concatenate projections and compute correlation with PCs
## This script reads in all the projected test vectors for each chromosome and adds them all together, and scales the final covariate
## It then treats this as the response variable and commputes the slope with all the PCs

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript code/calculate_FGr/compute_PC_FGR_model.R <prefix to Tm chromosomes> <file to PCs> <outfile> <type of E/O constrasts>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
PC_file = args[2]
outfile = args[3]
test_type = args[4]


# Get list of correct chromosomes
if (test_type == "all") {
  chrs <- seq(1,22)
} else if (test_type == "odd" ) {
  chrs <- c(1,3,5,7,9,11,13,15,17,19,21)
} else if (test_type == "even") {
  chrs <- c(2,4,6,8,10,12,14,16,18,20,22)
}
print(chrs)

# Add FGr together for all the blocks on chrtype chromosomes
df <- fread(paste0(Tm_prefix, "_", chrs[1], ".txt"))
for (i in 2:length(chrs)) {

  new <- fread(paste0(Tm_prefix, "_", chrs[i], ".txt"))
  df$FGr <- df$FGr + new$FGr

}

df$FGr <- scale(df$FGr)

# Join with PCs
PCs <- fread(PC_file)
dfCombine <- inner_join(PCs, df)
print(nrow(dfCombine))
print(head(dfCombine[,3]))
print(head(dfCombine[,42]))

# Construct data frame to collate results
dfOut <- matrix(NA, nrow = 40, ncol = 5)
dfOut[,1] <- seq(1,40)

# Compute multiple R^2 and rho(PC, FGr)
for (i in 1:nrow(dfOut)) {


  # Compute variance explained
  mod <- lm(dfCombine$FGr ~ . ,data=dfCombine[,3:(i+2)])
  r2 <- cor(dfCombine$FGr, fitted(mod))^2
  print(r2)

  # Compute correlation
  #name <- paste0("PC_", i)
  name <- paste0("PC", i)
  print(name)
  c1 <- dfCombine[[name]]
  print(length(c1))
  print(length(dfCombine$FGr))
  ct <- cor.test(as.numeric(dfCombine$FGr),as.numeric(c1))
  rho <- ct$estimate
  lc <- ct$conf.int[1]
  uc <- ct$conf.int[2]

  # Collect output
  dfOut[i,2] <- r2
  dfOut[i,3] <- rho
  dfOut[i,4] <- lc
  dfOut[i,5] <- uc

}


# Save output
dfOut <- as.data.frame(dfOut)
colnames(dfOut) <- c("PC", "r2", "rho", "lc", "up")
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








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
outfileData = args[5]

# Get test type
test_type <- strsplit(test_type, "-")[[1]][2]
print(test_type)

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
tmp <- fread(paste0(Tm_prefix, "_", chrs[1], ".txt"))
df <- tmp[,1:3]
df$FGr <- rowSums(tmp[,4:ncol(tmp)])
for (i in 2:length(chrs)) {

  tmp <- fread(paste0(Tm_prefix, "_", chrs[i], ".txt"))
  df$FGr <- df$FGr + rowSums(tmp[,4:ncol(tmp)])

}
df$FGr <- scale(df$FGr)

# Join with PCs
PCs <- fread(PC_file)
PC_IDs <- PCs[,1:2]
PC_nums <- PCs[,3:ncol(PCs)]

# Scale the PCs to have variance 1
PC_nums <- scale(PC_nums)
#print(paste0("The var of the PC cols are ", apply(PC_nums, 2, var)))

# Combine Dataframes
dfPCs <- cbind(PC_IDs, PC_nums)
dfCombine <- inner_join(df, dfPCs)
print(colnames(dfCombine))


# Construct data frame to collate results
totalPCs <-  (ncol(PCs) - 2)
print(totalPCs)
dfOut <- matrix(NA, nrow = totalPCs, ncol = 8)
dfOut[,1] <- seq(1,totalPCs)

# Compute multiple R^2 and rho(PC, FGr)
for (i in 1:nrow(dfOut)) {


  # Compute Bi
  name <- paste0("PC", i)
  print(name)
  c1 <- dfCombine[[name]]
  B <- cov(dfCombine$FGr, c1)
  B2 <- B^2

  # Compute correlation
  c1 <- dfCombine[[name]]
  ct <- cor.test(as.numeric(dfCombine$FGr),as.numeric(c1))
  lc <- ct$conf.int[1]
  uc <- ct$conf.int[2]


  # Collect output
  dfOut[i,2] <- B
  dfOut[i,3] <- B^2
  dfOut[i,4] <- lc
  dfOut[i,5] <- uc
  dfOut[i,6] <- sum(dfOut[1:i, 3])

}


# Save output
dfOut <- as.data.frame(dfOut)
colnames(dfOut) <- c("PC", "B", "B2", "lc", "uc", "r2")
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")

# Save other output
fwrite(dfCombine, outfileData, row.names = F, col.names = T, quote = F, sep = "\t")







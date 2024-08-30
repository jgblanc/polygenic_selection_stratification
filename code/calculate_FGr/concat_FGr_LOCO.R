## Concatenate projections
## This script reads in all the projected test vectors for each chromosome and adds them all together, leaving the focal chromosome out, and scales the final covariate

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript conccat_FGr_LOCO.R <prefix to Tm chromosomes> <chromosome number> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_prefix = args[1]
chrNum = as.numeric(args[2])
outfile = args[3]
print(chrNum)

# Get list of all but the focal chromosome
chrs <- seq(1,22)
chrs <- chrs[ !chrs == chrNum]

# Add Tm for each chromosome to each other
df <- fread(paste0(Tm_prefix, "_", chrs[1], ".txt"))
for (i in 2:length(chrs)) {

  new <- fread(paste0(Tm_prefix, "_", chrs[i], ".txt"))
  df$FGr <- df$FGr + new$FGr

}

df$FGr <- scale(df$FGr)

# Save output
fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








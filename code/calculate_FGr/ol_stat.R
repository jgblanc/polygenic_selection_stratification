## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

infile = args[1]
outfile = args[2]


df <- fread(infile)
df <- df[, c(1,2,4)]
df$FGr <- scale(df$FGr)
print(head(df))

fwrite(df, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








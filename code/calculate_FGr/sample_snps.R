# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

outAll = args[1]
outSample = args[2]
nsnp = as.numeric(strsplit(args[3], "L-")[[1]][2])

df <- fread(args[4], header=FALSE)
print(df)
for (i in 5:length(args)) {
  tmp <- fread(args[i], header=FALSE)
  df <- rbind(df, tmp)
}
fwrite(df, outAll ,row.names=F,quote=F,sep="\t", col.names = T)
print(nrow(df))

df <- df %>% sample_n(nsnp)
fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = T)




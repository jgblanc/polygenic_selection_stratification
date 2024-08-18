# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript get_IDs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

outAll = args[1]
outSample = args[2]
nsnp = strsplit(args[3], "L-")[[1]][2]

if (nsnp == "all") {
   df <- fread(args[26], header=FALSE)
   print(args[26])
   for (i in 27:length(args)) {
       tmp <- fread(args[i], header=FALSE)
       df <- rbind(df, tmp)
       }
   fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = T)
   fwrite(df, outAll ,row.names=F,quote=F,sep="\t", col.names = T)
} else if (nsnp == "pruneall") {
   df <- fread(args[4], header=FALSE)
   for (i in 5:25) {
       tmp <- fread(args[i], header=FALSE)
       df <- rbind(df, tmp)
       }
   fwrite(df, outAll ,row.names=F,quote=F,sep="\t", col.names = T)
   fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = T)
} else {
  nsnp = as.numeric(nsnp)
  df <- fread(args[4], header=FALSE)
  for (i in 5:25) {
      tmp <- fread(args[i], header=FALSE)
      df <- rbind(df, tmp)
      }
  fwrite(df, outAll ,row.names=F,quote=F,sep="\t", col.names = T)
  df <- df %>% sample_n(nsnp)
  fwrite(df, outSample ,row.names=F,quote=F,sep="\t", col.names = T)
}




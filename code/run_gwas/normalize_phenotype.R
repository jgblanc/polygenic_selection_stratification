# Script to format categorical covariates

args=commandArgs(TRUE)

if(length(args)<3){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

pheno = args[1]
ids = args[2]
outfile = args[3]

# Read in phenotypes
dfPheno <- fread(pheno) %>% drop_na()
colnames(dfPheno) <- c("#FID", "IID", "Raw")

# Read in IDs
dfIDs <- fread(ids)

# Join files
df <- inner_join(dfPheno, dfIDs, by = c("#FID", "IID"))

# Replace batch with array
df <- df %>% mutate(Value = qnorm((rank(Raw,na.last=NA)-0.5)/length(Raw))) %>% select("#FID", "IID", "Value")

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")


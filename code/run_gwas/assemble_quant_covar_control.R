# Script to format quantitative covariates when not including FGr

args=commandArgs(TRUE)

if(length(args)<4){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

age = args[1]
pc = args[2]
pcNum = as.numeric(args[3])
outfile = args[4]

# Read in dataframes
dfAge <- fread(age)
dfPc <- fread(pc)

# Subset to correct number of PCs
dfPc <- dfPc[,1:(pcNum + 2)]
print(ncol(dfPc))

# Join files
df <- inner_join(dfAge, dfPc, by = c("#FID", "IID"))

# Drop rows with NA
df <- drop_na(df)

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")


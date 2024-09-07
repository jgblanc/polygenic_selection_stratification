# Script to format categorical covariates

args=commandArgs(TRUE)

if(length(args)<3){stop("Provide path to covariate files")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

sex = args[1]
batch = args[2]
outfile = args[3]

# Read in phenotypes
dfSex <- fread(sex)
head(dfSex)
dfBatch <- fread(batch)
head(dfBatch)

# Join files
df <- inner_join(dfSex, dfBatch, by = c("#FID", "IID"))

# Replace batch with array
df <- df %>% mutate(genotype_measurement_batch_22000 = replace(genotype_measurement_batch_22000, genotype_measurement_batch_22000 < 0, 0)) %>%
  mutate(genotype_measurement_batch_22000 = replace(genotype_measurement_batch_22000, genotype_measurement_batch_22000 > 0, 1))

# Save file
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")


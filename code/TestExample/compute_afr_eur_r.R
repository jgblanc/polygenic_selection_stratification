# This script divides the WBS into GWAS an test panels

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript compute_afr-wbs.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

afr_file = args[1]
wbs_file  = args[2]
outfile = args[3]


## Read in freq files
dfAFR <- fread(afr_file)
dfWBS <- fread(wbs_file)

# Make new file
dfOut <- dfAFR %>% select("#CHROM", "ID", "REF", "ALT")
dfOut$r <- dfAFR$ALT_FREQS - dfWBS$ALT_FREQS

## Save
fwrite(dfOut, outfile ,row.names=F,quote=F,sep="\t", col.names = T)



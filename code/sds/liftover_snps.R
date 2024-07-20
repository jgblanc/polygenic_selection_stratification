## Overlapping SNPs
## This script takes a .afreq file from UKBB, subsets it to MAF > 1%, and intersects it with SDS snps (already limited to MAF > 5% in UK10K)
## It then flips the SDS value if the derived allele is the reference - this ensures that the SDS value is for the ALT allele which matches the rest of the pipeline


args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript liftover_snps.R <ukbb.freq> <SDS> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ukbb_file = args[1]
sds_file = args[2]
outfile_r = args[3]
outfile_overlap = args[4]

# Read in UKBB freq file
ukbb <- fread(ukbb_file)
ukbb <- ukbb %>% separate(ID, c("chr", "POS"), remove = FALSE) %>% filter(ukbb$ALT_FREQS > 0.01 & ukbb$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")

# Read in SDS
sds <- fread(sds_file)
sds <- sds %>% select("lifted_chr", "alt", "AA", "DA", "norm_SDS", "lifted_start")

# Wrangle columns
colnames(sds)[1] <- "#CHROM"
sds$`#CHROM` <- as.character(sds$`#CHROM`)
sds <- sds %>%
  mutate(`#CHROM` = str_replace(`#CHROM`, "chr", ""))
colnames(sds)[6] <- "POS"
sds$ID <- paste0(sds$`#CHROM`,":", sds$POS)
sds$`#CHROM` <- as.numeric(sds$`#CHROM`)
sds$POS <- as.numeric(sds$POS)

# Join files by chromosome and position
ukbb$POS <- as.numeric(ukbb$POS)
df <- inner_join(ukbb, sds, by = c("#CHROM" = "#CHROM", "POS" = "POS"))

# Get rid of rows where any alleles don't match
df_filter <- df %>% filter(ALT == alt)
df_filter2 <- df_filter %>% filter((DA == REF & AA == ALT) | (DA == ALT & AA == REF))

# If derived and alternate allele don't match, flip SDS
df_flipped <- df_filter2 %>% mutate(SDS = case_when(DA == REF ~ (-1 * norm_SDS), DA == ALT ~ norm_SDS))

# Select correct columns
df_out <- df_flipped %>% select("#CHROM", "ID.x", "REF", "ALT", "SDS")
colnames(df_out) <- c("#CHROM", "ID", "REF", "ALT", "r")

# Save output
fwrite(df_out,outfile_r, row.names = F, col.names = T, quote = F, sep = "\t")

# Get list of only snps
dfSNP <- df_out %>% select("ID")
fwrite(dfSNP ,outfile_overlap, row.names = F, col.names = T, quote = F, sep = "\t")

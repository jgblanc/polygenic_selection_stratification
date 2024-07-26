## Overlapping SNPs
## This script generates the set of constrasts for the Le et al data, as well as the overlapping SNPs

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript flip_overlapping_snps.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ukbb_file = args[1]
data_file = args[2]
r_outpath = args[3]
snp_outpath = args[4]
chr_num = args[5]
print(chr_num)


# Read in UKBB freq file
ukbb <- fread(ukbb_file)
ukbb <- ukbb %>% separate(ID, c("tmp", "POS"), remove = FALSE) %>% filter(ukbb$ALT_FREQS > 0.01 & ukbb$ALT_FREQS < 0.99) %>% select("#CHROM",	"ID","POS","REF",	"ALT")
colnames(ukbb)[1] <- "CHROM"
ukbb <- ukbb %>% separate(ID, c("tmp", "POS"), remove = FALSE)
ukbb$POS <- as.numeric(ukbb$POS)

# Read in contrasts
data <- fread(data_file)
data <- data %>% filter(CHROM == as.integer(chr_num))
print(nrow(data))

# Join files by chromosome and position
df <- inner_join(ukbb, data)
head(nrow(df))

# Change allele frequency to alternate allele freq by taking 1 -
df_flipped <- df %>% mutate(EUROPEAN_NEOLITHIC = 1 - EUROPEAN_NEOLITHIC, BRONZE_AGE = 1 - BRONZE_AGE, HISTORICAL = 1 - HISTORICAL,
                            ANATOLIA_NEOLITHIC = 1 -ANATOLIA_NEOLITHIC, MESOLITHIC = 1 - MESOLITHIC, STEPPE = 1 - STEPPE)

# Calculate statistic
df_EN <- df_flipped %>% mutate(expected = (0.84*ANATOLIA_NEOLITHIC) + (0.16*MESOLITHIC)) %>% mutate(ss = EUROPEAN_NEOLITHIC - expected)
df_BA <- df_flipped %>% mutate(expected = (0.52*STEPPE) + (0.48*EUROPEAN_NEOLITHIC)) %>% mutate(ss = BRONZE_AGE - expected)
df_H <- df_flipped %>% mutate(expected = (0.15*EUROPEAN_NEOLITHIC) + (0.85*BRONZE_AGE)) %>% mutate(ss = HISTORICAL - expected)

# Select correct columns and save output

df_EN <- df_EN %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_EN) <- c("#CHROM", "ID", "REF", "ALT", "r")
fwrite(df_EN, paste0(r_outpath,"EuropeanNeolithic_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")

df_BA <- df_BA %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_BA) <- c("#CHROM", "ID", "REF", "ALT", "r")
fwrite(df_BA, paste0(r_outpath,"BronzeAge_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")

df_H <- df_H %>% select("CHROM", "ID", "REF", "ALT", "ss")
colnames(df_H) <- c("#CHROM", "ID", "REF", "ALT", "r")
fwrite(df_H, paste0(r_outpath,"Historical_chr", chr_num, ".rvec"), row.names = F, col.names = T, quote = F, sep = "\t")

# Save set of overlapping SNPs
dfSNP <- df_flipped %>% select("ID")
fwrite(dfSNP, paste0(snp_outpath,"Historical/overlappingSNPs_chr", chr_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfSNP, paste0(snp_outpath,"BronzeAge/overlappingSNPs_chr", chr_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfSNP, paste0(snp_outpath,"EuropeanNeolithic/overlappingSNPs_chr", chr_num, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")

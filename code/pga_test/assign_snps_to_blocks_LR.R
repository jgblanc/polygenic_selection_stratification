# Script to assign snps to LD blocks and harmonize SNP directions
# .rvec is altFreq - altFreq
# A1 = ALT

args=commandArgs(TRUE)

if(length(args)<5){stop("Provide <ld block> <vilma format> <rvec> <output>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

ldFile =args[1]
ssFile = args[2]
rFile = args[3]
outfile = args[4]
pt = as.numeric(args[5])
#snpsfile = args[5]

# Read in LD block file
ld <- fread(ldFile)

# Read in summary stats
ss <- fread(ssFile)
colnames(ss)[4] <- "ALT"
colnames(ss)[5] <- "REF"
colnames(ss)[2] <- "ID"
print(head(ss))

# Read in r file
rvec <- fread(rFile)

# Combine r and ss
df <- inner_join(rvec, ss, by = c("ID", "REF", "ALT"))
df <- df %>% separate(ID, c("chr", "POS"), remove = FALSE)
print(head(df))
print(nrow(df))

# Assign SNPs to blocks

assign_SNP_to_block <- function(CHR, BP, block = ld) {

  # Filter blocks based on snp
  block_chr <- block %>% filter(chr == CHR)
  first_start <- as.numeric(block_chr[1, "start"])
  block_bp <- block_chr %>% filter( (start < BP & stop >= BP) | BP == first_start)

  # Assign
  block_num <- as.numeric(block_bp[,"block_number"])
  return(block_num)
}

# Remove SNPs below threshold
df <- df %>% filter(P <= pt)
print(nrow(df))

# Add block info - takes a while
df_blocks <- df %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[3]), as.numeric(params[4]))))

# Pick minimum p-value per block below a threshold
df_minP <- df_blocks %>% group_by(block) %>% slice_min(P, with_ties = F)
print(nrow(df_minP))

# Combine
#snps <- fread(snpsfile, header = F)
#colnames(snps) <- c("CHR","ID", "block")
#df <- inner_join(df, snps)

# Format output
out <- df_minP %>% select("CHR", "ID", "POS", "REF", "ALT", "r", "BETA","block", "P")
#out <- df %>% select("CHR", "ID", "POS", "REF", "ALT", "r", "BETA","block", "P")

# Save output
fwrite(out, outfile, row.names = F, col.names = F, quote = F, sep = "\t")






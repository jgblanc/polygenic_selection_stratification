## Project Tvec per chromosome
## This script projects the test vector from the test to the gwas panel using genotypes for a single chromosome

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript calc_FGr_blocks.R <gwas panel prefix> <output directory> <contrasts> <overlap snps> <output file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

gwas_prefix = args[1]
out_prefix = args[2]
r_file = args[3]
snps_file = args[4]
ldFile = args[5]
outfile = args[6]
outfile_snps = args[7]
gwas_IDs = args[8]
chr_num = as.numeric(args[9])
print(paste0("This Chr num is ", chr_num))


# Read in GWAS individuals
dfGWAS_IDs <- fread(gwas_IDs)
dfGWAS_IDs <- dfGWAS_IDs %>% select("#FID", "IID", "POP")
m <- nrow(dfGWAS_IDs)

# Read in and format r
r <- fread(r_file)
print(head(r))
dfSnps <- fread(snps_file)
colnames(dfSnps) <- "ID"
print(head(dfSnps))
r <- inner_join(r, dfSnps)
r <- r %>% dplyr::select("ID", "ALT", "r")
colnames(r) <- c("ID", "A1", "BETA")
print(head(r))


# Separate ID into CHR and BP
r_blocks <- r %>% separate("ID", into = c("CHR", "BP"), sep = ":", remove = FALSE) %>% filter(CHR == chr_num)
print(paste0("r has", nrow(r), " rows"))

# Read in LD block file
ld <- fread(ldFile)

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

# Add block info - takes a while
r_blocks <- r_blocks %>%
  mutate(block = apply(., MARGIN = 1, FUN = function(params)assign_SNP_to_block(as.numeric(params[2]), as.numeric(params[3])))) %>%
  drop_na()
print(paste0("Now r blocks has", nrow(r_blocks), " rows"))

# Select only SNPs that have a block
r <- r %>% filter(ID %in% r_blocks$ID)
print(paste0("Now r has", nrow(r), " rows"))

# Save as scoring weights
fwrite(r, paste0(out_prefix, "_scoringWeights.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

# Loop through blocks and calculate FGr for each block
numBlocks <- length(unique(r_blocks$block))
dfSNPs <- as.data.frame(matrix(NA, ncol = 2, nrow = numBlocks))
colnames(dfSNPs) <- c("Block", "nSNP")
for (i in 1:numBlocks) {

  # Get block num
  block_num <- unique(r_blocks$block)[i]
  print(paste0("The block number is ", block_num))

  # Select only snps on that block
  selected_snps <- r_blocks %>% filter(block == block_num) %>% select("ID")
  snps_file =  paste0(out_prefix,"_SNPs_", block_num, ".txt")
  fwrite(selected_snps, snps_file, row.names = F, col.names = T, quote = F, sep = "\t")

  # Save number of SNPs
  nsnp_in_block <- nrow(selected_snps)
  dfSNPs[i,1] <- block_num
  dfSNPs[i,2] <- nsnp_in_block

  # Compute FGr
  cmd_b <- paste("sh code/calculate_FGr/GWAS_score.sh",
                 gwas_prefix,
                 paste0(out_prefix, "_scoringWeights.txt"),
                 paste0(out_prefix,".gxt_tmp"), snps_file, gwas_IDs, sep = " ")
  print(cmd_b)
  system(cmd_b)

  # Read in FGr
  dfFGr = fread(paste0(out_prefix, ".gxt_tmp.sscore"))
  if (i == 1) {
    dfIDs <- dfFGr[,1:1]
  }
  FGr = as.matrix(dfFGr$BETA_SUM)

  # Format output
  col_name <- paste0("block_", block_num)
  dfIDs[[col_name]] <- FGr
  print(head(dfIDs))

  # Remove tmp files
  cmd <- paste("rm", paste0(out_prefix,"_SNPs_", block_num, ".txt"), paste0(out_prefix,".gxt_tmp*") , sep = " ")
  system(cmd)

}

# Combine IDs and FGr
colnames(dfIDs)[1] <- "IID"
dfOut <- inner_join(dfGWAS_IDs, dfIDs)

# Save output
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")
fwrite(dfSNPs, outfile_snps, row.names = F, col.names = T, quote = F, sep = "\t")


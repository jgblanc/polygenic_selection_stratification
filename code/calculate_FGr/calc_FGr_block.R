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
outfile = args[4]
outfile_snps = args[5]
gwas_IDs = args[6]
chr_num = as.numeric(args[7])
print(paste0("This Chr num is ", chr_num))


# Read in GWAS individuals
dfGWAS_IDs <- fread(gwas_IDs)
dfGWAS_IDs <- dfGWAS_IDs %>% select("#FID", "IID", "POP")

# Read in and format r
r <- fread(r_file)
L <- nrow(r)
print(paste0("L is ", L))
r <- r %>% filter(CHR == chr_num) %>% dplyr::select("ID", "ALT", "r")
colnames(r) <- c("ID", "A1", "BETA")
print(paste0("r had rows: ", nrow(r)))

# Get r with block info
r_blocks <- fread(r_file)
r_blocks <- r_blocks %>% filter(CHR == chr_num) %>% dplyr::select("ID", "ALT", "r", "block")
colnames(r_blocks) <- c("ID", "A1", "BETA", "block")
print(paste0("r blocks had rows: ", nrow(r)))

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
  print(paste0("The number of SNPs in block is ", nsnp_in_block))
  print(paste0("The variance of FGr block is ", var(FGr)))
  print(paste0("This should have variance 1 ", var(FGr * (1/sqrt(nsnp_in_block)))))
  print(paste0("Does this have var 1? ", var(FGr *  (1/sqrt(nsnp_in_block)) * (sqrt(100000/L)))))


  # Format output
  col_name <- paste0("block_", block_num)
  dfIDs[[col_name]] <- FGr


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


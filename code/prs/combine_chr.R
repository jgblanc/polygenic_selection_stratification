# Script to combine all chromosome from vilma output and readjust SNPs IDs to only include CHR:POS

args=commandArgs(TRUE)

if(length(args)<5){stop("Provide <number of total chrosomesome> <first prefix> <second prefix> <output> <list of chromosome numbers> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

chr_num = as.numeric(args[1])
p1 = args[2]
p2 = args[3]
outfile = args[4]


print(chr_num)
print(p1)
print(p2)

# Read in first data frame
fn <- paste0(p1, ".", args[5], p2)
print(fn)
df <- fread(fn)


if (chr_num > 1) {

  for (i in 2:chr_num) {

    print(paste0("Adding chr num ", args[4+i]))

    # Read in data frame
    fn <- paste0(p1, ".", args[4+i], p2)
    print(fn)
    tmp <- fread(fn)

    # Add to existing df
    df$NAMED_ALLELE_DOSAGE_SUM <- df$NAMED_ALLELE_DOSAGE_SUM + tmp$NAMED_ALLELE_DOSAGE_SUM
    df$SCORE1_SUM <- df$SCORE1_SUM + tmp$SCORE1_SUM

  }

}



# Save output
fwrite(df, outfile,col.names=T,row.names=F,quote=F,sep="\t")


## Find the estimation error of FGr using a Jackknife approach

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_error_jacknife.R <prefix to Tm chromosomes> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

betas_in = args[1]
infile_var = args[2]
outfile = args[3]


# Add Tm for each chromosome to each other
data <- fread(betas_in)
data <- data[,c(1,2,3,12)]
dfVar <- fread(infile_var)

df <- inner_join(dfVar, data)
print(head(df))

df$BETA <- df$BETA * (sqrt(df$Var)) 

D <- mean(df$BETA^2)
print(D)

dfOut <- as.data.frame(D)
fwrite(dfOut, outfile, row.names = F, col.names = T, quote = F, sep = "\t")








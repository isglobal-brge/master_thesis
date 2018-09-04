
sink("bigsnp_imputation.txt") # We used sink function to create a log file with all output of the code below
options(stringsAsFactors = F)

files <- read.table("files_to_impute.txt", header = T) # It reads the files to impute

library(bigsnpr)

NCORES <- 16 # Modify according to your computer/server 

# For each file to impute we read it and impute it with the function snp_fastImpute and save the output into
# a new file in the imputed directory

for (i in files) 
{
  bot$sendMessage(paste("Reading",i))
  
  bedfile <- paste(i,".bed",sep = "")
  tmpfile <- tempfile()
  snp_readBed(bedfile, backingfile = tmpfile)
  
  AGP <- snp_attach(paste0(tmpfile, ".rds"))
  
  
  G <- AGP$genotypes
  CHR <- AGP$map$chromosome
  POS <- AGP$map$physical.pos
  
  infos  <- tryCatch(snp_fastImpute(G, infos.chr = CHR, ncores = NCORES), error = function(e) NULL)
  if (is.null(infos))
  {
    print(paste("Error: ", i ,"Failed Imputation"))
  }
  else
  {
    G$code256 <- bigsnpr:::CODE_IMPUTE_PRED
    AGP$genotypes <- G
      
    output <- paste("./imputed/",i,"_imputed.bed", sep = "")
    snp_writeBed(AGP, output, ind.row = rows_along(G), ind.col = cols_along(G))
  }
}

sink()
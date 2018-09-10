setwd("/PRS/analysis/database_score")

# First, we load the necessary packages

library(snpStats)
library(SNPRelate)
library(SNPassoc)
library(rsq)

# Here we read data and extract the necessary information

SNP_file <- read.plink("../QC/SNP_file_10")

annot <- SNP_file$map
snps <- annot$snp.name    
geno <- SNP_file$genotypes

pheno <- as.numeric(!is.na(SNP_file$fam$father))

fam <- SNP_file$fam

snpgdsBED2GDS("../QC/SNP_file_10.bed", "../QC/SNP_file_10.fam", "../QC/SNP_file_10.bim", out="SNP_file_10_GDS")
genofile <- snpgdsOpen("SNP_file_10_GDS")

set.seed(123456)

# In our case in this directory we had csv files with SNPs list from the genes belonging to the different
# SFARI categories

snps.files <- dir(pattern = "snps")
cat <- sub("\\_.*", "", snps.files)

## Each category separated

rsquared <- c()
i <- 1

for (j in cat) 
{
  
# For each category we read the SNP list and intersect it with the SNPs in our data
# Then we tried 3 different prunning thresholds based on linkage desequilibrium
# Finally we computed the PRS with the function getScore and calculate the R^2 of
# the generalized linear model

  print(paste("Category: ", j, sep = ""))
  
  file.name <- paste(j, "_ensembl_snps.csv", sep = "")
  snp.a <- read.csv(file = file.name[1], header = T)

  snps.sel <- intersect(snps, snp.a[,1])
  
  ### LD 0.4
    
  snp.prune <- unlist(snpgdsLDpruning(genofile, 
                                      ld.threshold = 0.4, 
                                      snp.id = snps.sel))
  
  score <- getScore(geno, annot, snp.prune)

  a <- data.frame(score,factor(pheno))
  colnames(a) <- c("score","pheno")
  
  c <- glm(a$pheno~a$score, family="binomial")

  b <- rsq(c, adj=T, type="n")

  rsquared[i] <- b
  names(rsquared)[i] <- paste(j,"_LD_04", sep="")

  i <- i+1
  
  
  ### LD 0.2
    
  snp.prune <- unlist(snpgdsLDpruning(genofile, 
                                      ld.threshold = 0.2, 
                                      snp.id = snps.sel))
  
  score <- getScore(geno, annot, snp.prune)
  
  a <- data.frame(score,factor(pheno))
  colnames(a) <- c("score","pheno")
  
  c <- glm(a$pheno~a$score, family="binomial")

  b <- rsq(c, adj=T, type="n")

  rsquared[i] <- b
  names(rsquared)[i] <- paste(j,"_LD_02", sep="")

  i <- i+1


  ### No prunning

  score <- getScore(geno, annot, snps.sel)

  a <- data.frame(score,factor(pheno))
  colnames(a) <- c("score","pheno")

  c <- glm(a$pheno~a$score, family="binomial")

  b <- rsq(c, adj=T, type="n")

  rsquared[i] <- b
  names(rsquared)[i] <- paste(j,"_no_prun", sep="")

  i <- i+1
}

# In this case the validation is performed in the same way, you only need to run this code
# on the validation dataset

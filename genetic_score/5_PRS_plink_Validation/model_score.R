args <- commandArgs(trailingOnly=TRUE)

range <- args[1] # The range name is read

fam <- read.table("../QC/SNP_file_validation_10.fam") # Family data is read

name <- dir(pattern=paste("plink.",range,sep="")) 

score <- read.table(name, header=T) # score belonging to this range is read

sink(paste("model_",range,".txt", sep ="")) # we used sink function to save all output into a file

a <- data.frame(score$SCORE,factor(fam[,6]))
colnames(a) <- c("score","pheno")

library(rsq)

# Here we compute the glm for the phenotype against the score and obtain the rsquared of the model

c <- glm(a$pheno~a$score, family="binomial")
print(c)
rsq(c, adj=T, type="n")

sink()


# Finally we save the list of SNPs that are taken into account to compute this score

snp_score <- read.table(paste("snp_score_",range,".txt",sep=""))

bim <- read.table("../QC/SNP_file_validation_10.bim")

snp_data <- bim[bim[,2] %in% snp_score[,1],1:2]

write.table(snp_data, file=paste("snp_score_",range,".txt",sep=""), row.names=F, col.names=F)



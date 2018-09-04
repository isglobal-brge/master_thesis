fam <- read.table("../QC/SNP_file_10.fam")

library(rsq)

rsquared <- c()

files <- dir(pattern="plink.S")

j <- 1

# for each range we compute the rsquared of the phenotype against the score

for(i in files)
{
	score <- read.table(i, header=T)

	range <- strsplit(x=i, split="\\.")[[1]][2]

	a <- data.frame(score$SCORE,factor(fam[,6]))
	colnames(a) <- c("score","pheno")

	c <- glm(a$pheno~a$score, family="binomial")
	print(c)
	b <- rsq(c, adj=T, type="n")

	rsquared[j] <- b
	names(rsquared)[j] <- range

	j <- j+1
}

# we make a barplot of the rsquared of all the models

pdf("rsquared_scores_barplot.pdf")
barplot(rsquared, las=2)
dev.off()

# We save the rsquared in a file
write.table(cbind(names(rsquared), rsquared), file="rsquared_scores.txt", col.names=F, row.names=F)
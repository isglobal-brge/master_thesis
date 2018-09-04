library(bigsnpr)
set.seed(123456)

# Here we read data and extract the necessary data to further steps

bedfile <- "../SNP_file_Imputed.bed"
tmpfile <- tempfile()
snp_readBed(bedfile, backingfile = tmpfile)

SNP_file <- snp_attach(paste0(tmpfile, ".rds"))

G <- SNP_file$genotypes
CHR <- SNP_file$map$chromosome
POS <- SNP_file$map$physical.pos

NCORES <- 16

# We also read the files generated in the original PRS analysis

betas <- read.table("betas.txt", header=T)
lps <- read.table("bigsnpr_pvals_to_new_data.txt", header=T)
thrs <- read.table("bigsnpr_ranges_to_new_data.txt")

# We intersect the SNP of the previous analysis with the ones in this file

lps_ok <- lps[lps$marker.ID %in% SNP_file$map$marker.ID,]
betas_ok <- betas[betas$marker.ID %in% lps_ok$marker.ID,]

snps_intersect <- which(SNP_file$map$marker.ID %in% lps_ok$marker.ID)

# We create a new genotype file with only those SNP and compute the PRS

SNP_file_prueba <- subset(SNP_file, ind.row = rows_along(G), ind.col = snps_intersect)
SNP_file_2 <- snp_attach(SNP_file_prueba)

G_2 <- SNP_file_2$genotypes

prs <- snp_PRS(G_2, betas.keep = betas_ok$estim,
                 lpS.keep = (-1*lps_ok$lpS)+1, 
                 thr.list = -1*thrs[,1])


nb.pred <- sapply(-1*thrs[,1], function(thr) sum((-1*lps_ok$lpS)+1 > thr))

y01 <- SNP_file$fam$affection - 1

# We compute the rsquared of each model by PRS

rsquared <- c()

for(i in 1:ncol(prs))
{

    a <- data.frame(prs[,i],factor(y01))
    colnames(a) <- c("score","pheno")

    library(rsq)

    c <- glm(a$pheno~a$score, family="binomial")
    print(c)
    b <- rsq(c, adj=T, type="n")

    rsquared[i] <- b
    names(rsquared)[i] <- nb.pred[i]
}

# We plot a barplot

pdf("rsquared_bigsnpr_scores_barplot.pdf")
barplot(rsquared, las=2)
dev.off()

# We save the rsquared data into a file

write.table(cbind(names(rsquared), rsquared), file="rsquared_bigsnpr_scores.txt", col.names=F, row.names=F)

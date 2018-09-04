library(bigsnpr)

# In case you need to install the package this is the code to do so
# devtools::install_github("privefl/bigsnpr") 
set.seed(123456)

# Here we read data and extract the necessary data to further steps

bedfile <- "SNP_file_Imputed.bed"
tmpfile <- tempfile()
snp_readBed(bedfile, backingfile = tmpfile)

SNP_file <- snp_attach(paste0(tmpfile, ".rds"))


G <- SNP_file$genotypes
CHR <- SNP_file$map$chromosome
POS <- SNP_file$map$physical.pos

NCORES <- 16


# This package recommends its own protocol of prunning and clumping SNP
# The first step is to identify those SNP in long-range LD regions and exclude them

ind.excl <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)

ind.keep <- snp_pruning(G, infos.chr = CHR,
                         exclude = ind.excl,
                         ncores = NCORES)

# In the next step it computes a singular value decomposition (SVD) to use them as
# covariates for the gwas training

obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(),
                         ind.col = ind.keep,
                         ncores = NCORES, verbose=T)

# Then we split the individual into training and testing groups

ind.train <- sort(sample(nrow(G), 0.8 * nrow(G)))
ind.test <- setdiff(rows_along(G), ind.train)

# We set a variable with the phenotypes

y01 <- SNP_file$fam$affection - 1

# The next step performs a column-wise logistic regression that we will use as
# betas for each SNP to compute the scores

gwas.train <- big_univLogReg(G, y01.train = y01[ind.train],
                             ind.train = ind.train, 
                             covar.train = obj.svd$u[ind.train, ], 
                             ncores = NCORES)

# Then we perform a last clumping of SNPs

ind.keep <- snp_clumping(G, infos.chr = CHR,
                         ind.row = ind.train,
                         S = abs(gwas.train$score),
                         size = 500,
                         is.size.in.bp = TRUE,
                         infos.pos = POS,
                         ncores = NCORES)

# We create different thresholds to the score and predict the values that will be
# taken into a account by the thresholds to perform the PRS

thrs <- c(0, exp(seq(log(0.2), log(550), length.out = 50)))

lpS <- -predict(gwas.train)

# Finally, we compute the PRS with the function snp_PRS and passing it as arguments
# the genotype, the betas, the testing individuals, the SNP to keep, the values to group
# the SNP by the thresholds, and the thresholds

prs <- snp_PRS(G, betas.keep = gwas.train$estim[ind.keep],
                 ind.test = ind.test,
                 ind.keep = ind.keep,
                 lpS.keep = lpS[ind.keep], 
                 thr.list = thrs)

# We compute the aucs for each PRS and the number of predictor used in each of them
# then we select only those which contain some predictors and exclude the ones without
# any SNP in it

aucs <- apply(prs, 2, AUC, target = y01[ind.test])

nb.pred <- sapply(thrs, function(thr) sum(lpS[ind.keep] > thr))

nb.pred_ok <- nb.pred[-(grep("^0$",nb.pred)[-1])]
prs_ok <- prs[,-1*(grep("^0$",nb.pred)[-1])]

setwd("./scores/")

# We compute the rsquared of the models of the phenotype against the scores

rsquared <- c()

for(i in 1:ncol(prs_ok))
{

    a <- data.frame(prs_ok[,i],factor(y01[ind.test]))
    colnames(a) <- c("score","pheno")

    library(rsq)

    c <- glm(a$pheno~a$score, family="binomial")
    print(c)
    b <- rsq(c, adj=T, type="n")

    rsquared[i] <- b
    names(rsquared)[i] <- nb.pred_ok[i]

}

# We plot a barplot

pdf("rsquared_bigsnpr_scores_barplot.pdf")
barplot(rsquared, las=2)
dev.off()

# We save the rsquared data into a file

write.table(cbind(names(rsquared), rsquared), file="rsquared_bigsnpr_scores.txt", col.names=F, row.names=F)

# For each score we save a file with the SNP id of the SNP used in each score
for(i in 1:ncol(prs_ok))
{
   snp.ind.keep <- ind.keep[which(lpS[ind.keep] > thrs[i])] 
   snp_id <- SNP_file$map[snp.ind.keep,1:2]
   write.table(snp_id, file=paste("snp_score_",nb.pred[i],".txt",sep=""), col.names=T, row.names=F)
}


# In this section we prepared the necessary data to validate this scores

betas <- cbind(SNP_file$map[,1:2], gwas.train) # we save the betas

write.table(betas, file="betas.txt", col.names=T, row.names=F)

pvals <- cbind(SNP_file$map[,1:2], lpS) # we save the "p-values" used by the thresholds

write.table(pvals, file="pvals.txt", col.names=T, row.names=F)

# We assign "false" p-values to each score to validate them later

snps <- read.table(paste("snp_score_",nb.pred_ok,".txt",sep=""),header=T)

pvals_ok <- pvals[pvals$marker.ID %in% snps$marker.ID,]

for (i in nb.pred_ok) 
{
    snps <- read.table(paste("snp_score_",i,".txt",sep=""),header=T)
    pvals_ok$lpS[pvals_ok$marker.ID %in% snps$marker.ID] <- i
}

write.table(pvals_ok,file="bigsnpr_pvals_to_new_data.txt", col.names=T, row.names=F)

# we create a range file with the new ranges based on the values given to the "false" p-values

write.table(nb.pred_ok, file="bigsnpr_ranges_to_new_data.txt", row.names=F, col.names=F)
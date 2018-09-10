options(stringsAsFactors = F)
bim <- read.table("/PRS/validation/Original_data/SNP_file_validation.bim")
logor <- read.table("LOGOR_to_score.txt", header=T)

logor$A1 <- as.character(logor$A1)

bim <- bim[,c(2,5)]

logor$A1[logor$SNP %in% bim$V2] <- bim$V5[bim$V2 %in% logor$SNP] # We change the coding of the SNPs


write.table(logor, "logor_new_data.txt", col.names=T, row.names=F)
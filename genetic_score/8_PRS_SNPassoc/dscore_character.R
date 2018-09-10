dscore.character <- function(x, ...){
  
  snpmart <- biomaRt::useEnsembl(biomart="snp", 
                                 dataset="hsapiens_snp")
  
  
  snpInfo <- biomaRt::getBM(c("refsnp_id", "chr_name", 
                              "chrom_start", 
                              "allele", 
                              "minor_allele", 
                              "minor_allele_freq"),
                            filters = c("snp_filter"),
                            values = x, mart = snpmart)
  
  ## Insertar que solo coja chromosomas 1-23, porque sino hay snps repetidos
  ## y cuando los quiere poner como rownames da error
  
  snpInfo <- snpInfo[snpInfo[,2] %in% c(1:23),]
  snpInfo<- snpInfo[!duplicated(snpInfo$refsnp_id),]
  rownames(snpInfo) <- snpInfo$refsnp_id
  
  if (nrow(snpInfo)!=length(x)) {
    warning("Genetic score distribution has been computed with these SNPs: \n",
            paste(unique(snpInfo$refsnp_id), collapse="; "))
  }
  x.info <- x[x%in%snpInfo$refsnp_id]    
  snpInfo <- snpInfo[x.info,]
  snpInfo<- snpInfo[!duplicated(snpInfo$refsnp_id),]

  rownames(snpInfo) <- snpInfo$refsnp_id
  mafs<-snpInfo$minor_allele_freq
  mafs <- na.omit(mafs)
  ans <- dscore(mafs)
  attr(ans, "MAFs") <- snpInfo
  ans
}
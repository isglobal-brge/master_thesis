getScore <- function(geno, annot, snp.sel){
  
  x <- geno[,snp.sel]
  annot <- annot[snp.sel, ]
  snp.maf <-  dscore.character(snp.sel)
  maf <- attr(snp.maf, "MAFs")
  
  
  ff <- function(i, x, flip){
    xx <- x[,i]
    flip.i <- flip[i]
    
    ans <- as.numeric(xx) - 1
    ans[ans == -1] <- NA
    
    out <- ans
    if(flip.i) {
      out[ans==2] <- 0
      out[ans==0] <- 2
    }
    out
  }
  
  getComp <- function(x){
    out <- x
    out[x=="G"] <- "C"
    out[x=="C"] <- "G"
    out[x=="A"] <- "T"
    out[x=="T"] <- "A"
    out
  }
  
  maf.ref <- maf$minor_allele
  maf.obs <- cbind(annot$allele.1, annot$allele.2)
  eq1 <- sweep(maf.obs, 1, FUN="==", maf.ref)
  eq2 <- sweep(maf.obs, 1, FUN="==", getComp(maf.ref))
  eq <- cbind(eq1, eq2)
  id <- apply(eq, 1, function(x) which(x)[1])
  flip <- id%in%c(1,3)
  
  xx <- data.frame(lapply(1:ncol(x), ff, x=data.frame(x), flip=flip))
  colnames(xx) <- colnames(x)
  
  ans <- rowSums(xx, na.rm = T) ## insertar na.rm = 
  
  ans
  
}
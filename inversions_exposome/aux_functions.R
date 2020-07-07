########################################
########## AUXILIAR FUNCTIONS ##########
########################################

library(readxl)
library(dplyr)


load("/home/isglobal.lan/ncarreras/homews/inversionGR.rda")

#Import the excel table to get the exposure abreviations
excel_table <- read_excel("/home/isglobal.lan/ncarreras/homews/Expodata.xlsx",sheet = "Exposure_Covariate")
exp_abrev <- cbind(excel_table[,"Variable_name_TRANS"],excel_table[,"Label for tables"],excel_table[,"Group"],excel_table[,"Subgroup"])
rownames(exp_abrev) <- exp_abrev[,"Variable_name_TRANS"]
rm(excel_table)
gc(reset=TRUE)

region_inv <- function(inv){
  start <- as.numeric(start(inversionGR[inv,]))-1000000
  end <- as.numeric(end(inversionGR[inv,]))+1000000
  region <- paste(as.character(seqnames(inversionGR[inv,])),":",start,"-",end,sep="")
  return(region)
}

allsymbols_to_symbol <- function(allsymbols) {
  symbols<-strsplit(allsymbols,";")[[1]]
  symbols[symbols=="NA" | symbols==""]=NA
  symbols<-unique(symbols[!is.na(symbols)])
  return (paste(symbols, collapse=";"))
}

split_symbol <- function(topcpgs) {
  a <- paste(topcpgs$Gene_Symbol,collapse=";")
  b <- unique(strsplit(a,";")[[1]])
  c <- b[which(b!="")]
  return (c)
}

select_top_genes <- function(inversion,type){
  trans <- sig_genes_trans[which(sig_genes_trans$Inversion==inversion & sig_genes_trans$Type_analysis==type),]
  selecttrans <- data.frame(Inversion=character(),
                            Type_analysis=character(), 
                            Transcript=character(),
                            Location=character(),
                            Gene_Symbol=character(),
                            model=character(),
                            adj.p.value=numeric(),
                            DiffLevene_logFC=numeric(),
                            stringsAsFactors=FALSE)
  for (transcript in trans$Transcript){
    transcriptsDF <- trans[which(trans$Transcript==transcript),]
    bestmodel <- transcriptsDF[which.min(transcriptsDF$adj.p.value),]
    selecttrans <- rbind(selecttrans,bestmodel)
  }
  selecttrans <- selecttrans %>% distinct()
  return(selecttrans)
}

select_top_cpgs <- function(inversion,type){
  methy <- sig_cpgs_methy[which(sig_cpgs_methy$Inversion==inversion & sig_cpgs_methy$Type_analysis==type),]
  selectmethy <- data.frame(Inversion=character(),
                            Type_analysis=character(), 
                            CpG=character(),
                            Location=character(),
                            Gene_Symbol=character(),
                            model=character(),
                            adj.p.value=numeric(),
                            DiffLevene_logFC=numeric(),
                            stringsAsFactors=FALSE)
  for (cpg in methy$CpG){
    cpgsDF <- methy[which(methy$CpG==cpg),]
    bestmodel <- cpgsDF[which.min(cpgsDF$adj.p.value),]
    selectmethy <- rbind(selectmethy,bestmodel)
  }
  selectmethy <- selectmethy %>% distinct()
  return(selectmethy)
}

diff_inv <- function(dataset,model_inherit,inversion,type_analysis){
  if (type_analysis=="Mean"){
    results <- runDiffMeanAnalysis(set= dataset,model= formula(paste("~",model_inherit,"(",inversion,")+ sex + cohort"),sep=""))
    features <- getProbeResults(results, rid="DiffMean",fNames = NULL)
    topfeatures<-features[ which (features$adj.P.Val < 0.05/2.2), ]
  }
  if (type_analysis=="Variance"){
    results <- runDiffVarAnalysis(set = dataset,model = formula(paste("~",model_inherit,"(",inversion,")+ sex + cohort"),sep=""))
    features <- getProbeResults(results, rid="DiffVar",fNames = NULL)
    topfeatures<-features[ which (features$Adj.P.Val < 0.05/2.2), ]
  }
  return(topfeatures)
}

diff_inv_expo <- function(dataset, type, inversion, expo){
  if (type=="Mean"){
    resmean <- runDiffMeanAnalysis(set = dataset,model = formula(paste("~",inversion,"*",expo,"+ sex + cohort"),sep=""))
    getresmean <- getProbeResults(resmean, rid="DiffMean", fNames = NULL, coef=4)
    topfeatures <- getresmean[which(getresmean$adj.P.Val < 0.05),]
  }
  if (type=="Variance"){
    resvar <- runDiffVarAnalysis(set = dataset,model = formula(paste("~",inversion,"*",expo,"+ sex + cohort"),sep=""))
    getresvar <- getProbeResults(resvar, rid="DiffVar", fNames = NULL, coef=4)
    topfeatures <- getresvar[which(getresvar$Adj.P.Val < 0.05),]
  }
  return(topfeatures)
}

symbolToentrez_ensembl<-function(symbols,entrez_ensemble) {
  symbols<-paste(symbols, collapse=";")
  symbols<-unique(strsplit(symbols,";")[[1]])
  genes <- bitr(symbols, fromType = "SYMBOL",
                toType = entrez_ensemble,
                OrgDb = "org.Hs.eg.db")
  return (genes)
}

#Enrichment with KEGG
kegg<-function(genes_entrez){
  ans.kegg <- enrichKEGG(gene = genes_entrez$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.05)
  return(ans.kegg)
}

#Enrichment with DisGeNET
dis<-function(genes_entrez){
  ans.dis <- enricher(genes_entrez$ENTREZID, TERM2GENE=disease2gene,
                      TERM2NAME=disease2name)
  return(ans.dis)
}

#Enrichment with GO
GO <- function(genes_ensembl){
  ans.go <- enrichGO(gene = genes_ensembl$ENSEMBL, ont = "BP",
                     OrgDb ="org.Hs.eg.db",
                     keyType = "ENSEMBL",
                     readable=TRUE,
                     pvalueCutoff = 0.05)
  return(ans.go)
}

#Function to perform a multimediation analysis
multimed <- function(medi_dataset,inversion,topgenes,topcpgs){
  medi_genes <- colnames(medi_dataset)[grep("^TC", colnames(medi_dataset))]
  medi_cpgs <- colnames(medi_dataset)[grep("^cg", colnames(medi_dataset))]
  for (gene in medi_genes){
    E <- as.numeric(medi_dataset[[inversion]])
    M <- medi_dataset[,medi_cpgs]
    Y <- medi_dataset[[gene]]
    a <- medTest(E, M, Y, nperm = 1000)
    rownames(a) <- medi_cpgs
    if (nrow(a)!=0){
      for (i in 1:nrow(a)){
        de <- data.frame(Inversion=inversionGR[inversion,]$Cytogenetic.location,
                         Transcript=gene,
                         Location.Transcript=topgenes[which(topgenes$Transcript==gene),"Location"],
                         Symbol.Transcript=topgenes[which(topgenes$Transcript==gene),"Gene_Symbol"],
                         CpG=rownames(a)[i],
                         Location.CpG=topcpgs[which(topcpgs$CpG==rownames(a)[i]),"Location"],
                         Symbol.CpG=topcpgs[which(topcpgs$CpG==rownames(a)[i]),"Gene_Symbol"],
                         p.value.mediation=a[i,"p"])
        inv_cpg_gene_var <- rbind(inv_cpg_gene_var,de)
      }
    }
  }
  return(inv_cpg_gene_var)
}

#Function to study a specific mediation statistically
mediation <- function(num,medi_dataset,inversion){
  transcript <- as.character(inv_cpg_gene_var[num,"Transcript"])
  cpg <- as.character(inv_cpg_gene_var[num,"CpG"])
  print(paste(transcript,cpg,sep=" vs "))
  design1 <- formula(paste(transcript,"~ additive(",inversion,") + sex + cohort"))
  mod1 <- glm(design1, data=medi_dataset)
  design.M <- formula(paste(cpg,"~ additive(",inversion,") + sex + cohort"))
  model.M <- glm(design.M, data=medi_dataset)
  design.Y <- formula(paste(transcript,"~ additive(",inversion,") +",cpg," + sex + cohort"))
  model.Y <- glm(design.Y, data=medi_dataset)
  treatment=paste("additive(",inversion,")",sep="")
  res <- mediate(model.M, model.Y, treat=treatment, mediator=cpg)
  print(summary(res))
  inv_cpg_gene_var$Prop.Mediated[num] <- summary(res)$n0
  inv_cpg_gene_var$p.value.prop.mediated[num] <- summary(res)$n0.p
  return(inv_cpg_gene_var)
}

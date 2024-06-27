# Supplementary Figure 1 - X-Ra in different breast cancer subtypes

#Load the packages
library(GEOquery) 
library(chrXRa) 
library(tidyverse) 
library(ggpubr)

# Figure based on the dataset: GSE225845
path_here <- paste(getwd(), "/GSE225845", sep = "")

#Get the series
gse <- getGEO(ID, destdir = path_here)[[1]] 

#Get the info
pheno <- pData(phenoData(gse)) 

#Not methylation info included, download it from GEO manually or using bash: 
# wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225845/suppl/GSE225845_tumors_normalized_betas.txt.gz" -O "GSE225845_tumors_normalized_betas.txt.gz"

#Read the methy of ONLY cancer patients 
library(readr) 
met <- read_tsv("GSE225845_tumors_normalized_betas.txt") #CpGs in columns and individuals in rows 
met <- as.data.frame(met) 
#Basenames are the same in both 
rownames(met) <- met$basename 
rownames(pheno) <- pheno$`methylation id (basenames):ch1` 

#Common individuals. Careful! The dataset has the methylation matrix transposed (features in columns and individuals in rows)
common <- intersect(rownames(met), rownames(pheno)) 
met <- met[common, ] 
pheno <- pheno[common, ] 

#Only female 
fem <- pheno$`Sex:ch1` == "F" 
pheno <- pheno[fem, ] 
met <- met[fem, ] 

#Annotation neither included
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) 
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) 
common_annot <- intersect(colnames(met), rownames(annotation)) 

met2 <- met[, common_annot] 
annotation <- annotation[common_annot, ] 
paste("Identical CpGs: ", identical(colnames(met2), rownames(annotation)), sep = "") #False 

#Estimate XRa
sel <- annotation$chr == "chrX" 
metfemale <- (as.matrix(met2[, sel]))  
#Estimate XRa values 
pheno$XRa <- XRa(metfemale, colnames(metfemale)) 

plt <- ggboxplot(data = pheno[!is.na(pheno$molecular.subtype.ch1),] , 
                 x = "molecular.subtype.ch1", 
                 y = "XRa",
                 fill = "molecular.subtype.ch1")
facet(p = plt, facet.by = "race.ch1") +
  labs(title = "X-Ra values between different breast cancers (GSE225845)", x = "Breast cancer subtype", y = "X-Ra") +
  scale_fill_manual(name = "Breast cancer subtype", values = c("purple", "grey", "green2")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16), axis.title.y.right = element_blank(), legend.position = "bottom", strip.text = element_text(size = 14))
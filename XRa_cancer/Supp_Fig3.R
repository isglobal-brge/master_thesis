# Supplementary Figure 3 - X-Ra in different cancer cell lines.

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)

path_here <- paste(getwd(), "/GSE124368", sep = "")
#Get the series
gse <- getGEO("GSE124368", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse)) #All females.
annotation <- fData(gse)

sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

pheno$title <- sub("Genomic DNA from ", "Genomic DNA from\n", pheno$title)
ggplot(data = pheno, aes(x = title, y = XRa))+
  geom_point(size = 4) +
  labs(title = "Breast cancer cell lines and  X-Ra (GSE124368)", x = "Cell line", y = "X-Ra") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))
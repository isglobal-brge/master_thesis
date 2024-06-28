# Supplementary Figure 5 - Children X-Ra and BMI. 

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggpubr)
library(cowplot)

path_here <- paste(getwd(), "/BMI/GSE72556", sep = "")
#Get the series
gse <- getGEO("GSE72556", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse)) #All females.
annotation <- fData(gse)

sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

pheno$`child bmi:ch1` <- round(as.numeric(pheno$`child bmi:ch1`), 2)
pheno$child_BMI <- pheno$`child bmi:ch1`
pheno$`adult bmi:ch1`  <- round(as.numeric(pheno$`adult bmi:ch1`), 2)
pheno$adult_BMI <- pheno$`adult bmi:ch1`

#We can compare if the BMI of the mother influences the XRa values of their children.

plot1 <- ggscatter(data = pheno, 
                   x = "child_BMI", 
                   y = "XRa", 
                   add = "reg.line",
                   col = "green4", 
                   size = 3) +
  ggpubr::stat_cor(method = "spearman", 
                   cor.coef.name = "rho",
                   label.x.npc = "middle", 
                   col = "green4",
                   size = 6) +
  labs(title = "X-Ra in children's saliva and their own BMI (GSE72556)", x = "Child's BMI index", y = "X-Ra") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))

plot2 <- ggscatter(data = pheno, 
                   x = "adult_BMI", 
                   y = "XRa", 
                   add = "reg.line",
                   col = "cyan3",
                   size = 3) +
  ggpubr::stat_cor(method = "spearman", 
                   cor.coef.name = "rho", 
                   label.x.npc = "middle", 
                   col = "cyan3",
                   size = 6) +
  labs(title = "X-Ra in children's saliva and mother's BMI (GSE72556)", x = "Mother's BMI index", y = "") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))

plot_grid(plot1, plot2)

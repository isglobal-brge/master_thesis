# Figure 9 - Dynamics and associations between X-Ra and treatment outcome.

#Load the required packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(ggpubr)

#Figures based on GSE184159

path_here <- paste(getwd(), "/GSE184159", sep = "")
#Get the series
gse <- getGEO("GSE184159", destdir = path_here)[[1]]
#Get the info
pheno <- pData(phenoData(gse))

#Methylation data nor annotations included in the dataset.


annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

#Methylation, downloaded from GEO entry manually or using bash: wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184159/suppl/GSE184159_beta_detp.csv.gz" -O "GSE184159_beta_detp.csv"
met <- read.csv(paste(path_here, "/GSE184159_beta_detp.csv", sep = ""))
rownames(met) <- met$X 
met <- met[, grep("^beta", colnames(met))]
colnames(met) <- pheno$geo_accession
#We have to reorder the rows
same_rows <- intersect(rownames(met), rownames(annotation))
#Subset
met <- met[same_rows, ]
annotation <- annotation[same_rows, ]
identical(rownames(met), rownames(annotation)) 

#We write all the tissue names to identify which are the same ones.
pheno$IDs <- gsub("_.*_.", "", pheno$title)
#Select CpGs only in chrX.
sel <- as.vector(annotation$chr) == "chrX"
# Select CpGs located in the X chr.
metfemale <- t(as.matrix(met[sel, ])) 
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

## a) X-Ra comparison between treatment timepoints (GSE174159)


pheno <- pheno %>%
  mutate(time_point = pheno$`timepoint:ch1`, 
         response = pheno$`response:ch1`)

ggpaired(data = pheno, 
         x = "time_point", 
         y = "XRa", 
         id = "IDs", 
         color = "black", 
         line.color = "IDs", 
         legend = "none") +
  labs(title = "X-Ra comparison between treatment timepoints (GSE184159)", x = "Timepoint", y = "X-Ra") +
  stat_compare_means(paired = TRUE, size = 5) +
  scale_x_discrete(labels = c("Control", "Before treatment", "Midpoint", "After treatment")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.position = "None")

## b) Response comparison depending on X-Ra in Breast cancer (GSE184159)e

pheno_cancer <- pheno %>%
  select(source_name_ch1 == "Breast_Cancer") %>%
  mutate(num_time = factor(ifelse(`timepoint:ch1` == "A", 1, ifelse(`timepoint:ch1` == "B", 2, 3))))

ggscatter(data = pheno_cancer, 
          x = "num_time", 
          y = "XRa", 
          color = "response",
          palette = c("cyan4", "purple", "brown3"),
          add = "reg.line",
          point = FALSE) +
  stat_cor(aes(color = response), 
           method = "spearman", 
           cor.coef.name = "rho",
           label.x.npc = "left", 
           label.y.npc = "top",
           size = 5) +
  geom_jitter(data = pheno_cancer, aes(x = num_time, y = XRa, col = response)) +
  geom_smooth(aes(col = response, fill = response), formula = y ~ x, method = "lm", alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Before treatment", "Midpoint", "After treatment")) +
  scale_fill_manual(values = c("cyan4", "purple", "brown3")) +
  labs(title = "Response comparison depending on X-Ra in Breast cancer (GSE184159)", x = "Timepoint", y = "X-Ra") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))
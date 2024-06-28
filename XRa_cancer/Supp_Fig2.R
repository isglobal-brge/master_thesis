# Supplementary Figure 2: X-Ra in ovarian and associated tissue. 

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggstatsplot)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# The image is done using the dataset  (GSE155760)

path_here <- paste(getwd(), "/GSE155760", sep = "")

#Get the series
gse <- getGEO("GSE155760", destdir = path_here)[[1]]

#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#Test if annotations and met CpGs are in the same order
common <- intersect(rownames(met), rownames(annotation))
met <- met[common, ]
annotation <- annotation[common, ]

#Select CpGs only in chrX and estimate X-Ra.
sel <- annotation$CHR == "chrX"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

# Plot

pheno$tumour_type <- sub("from cancer-free normal control", "control", pheno$histology.ch1)
pheno$tumour_type <- sub(" ", "\n", pheno$tumour_type)
pheno$tumour_type <- sub("serous ovarian carcinoma", "serous ovarian\ncarcinoma", pheno$tumour_type)
pheno$tumour_type <- sub("poorly\ndifferentiated ovarian", "poorly differenciated\novarian", pheno$tumour_type)
ggbetweenstats(
  data = pheno, 
  x = tumour_type , 
  y = XRa, 
  bf.message = FALSE,
  results.subtitle = FALSE,
  k = 3L,
  p.adjust.method = "bonferroni", 
  centrality.label.args = list(size = 4.5, nudge_x = 0.4, segment.linetype = 4, min.segment.length = 0),
  ggsignif.args = list(textsize = 4.5, tip_length = 0.01, na.rm = TRUE)
) +
  scale_colour_discrete(name = "Tumour type") +
  labs(title = "X-Ra values between different ovarian and associated samples (GSE155760)", x = "Ovarian cancer", y = "X-Ra") + 
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16), axis.title.y.right = element_blank(), legend.position = "bottom")
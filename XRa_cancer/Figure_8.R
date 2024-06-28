# Figure 8 - Association between X-Ra and breast cancer in breast tissue

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggstatsplot)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

## All the plots based on GSE237036

path_here <- paste(getwd(), "/GSE237036", sep = "")

#Get the series
gse <- getGEO("GSE237036", destdir = path_here)[[1]]

#Get the info
pheno <- pData(phenoData(gse))

#Methylation data nor annotation available directly downloading the GEO series object.

# Download using bash or the GEO entry: wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237036/suppl/GSE237036_matrix_processed.txt.gz@ -O "GSE237036_matrix_processed.txt.gz"

met <- read.delim(paste(path_here, "/GSE237036_matrix_processed.txt.gz", sep = ""), sep = "\t") #Load the methylation data, adjust the path.
rownames(met) <- met$ID_REF
met <- met[-1] 

#Get the annotation data (EPIC array)

annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

#Reorther the methylation data to have all the annotations correctly and ready to subset.
identical(rownames(met), rownames(annotation))
same_rows <- intersect(rownames(met), rownames(annotation))
met <- met[same_rows, ]
annotation <- annotation[same_rows, ]

sel <- annotation$chr == "chrX"

## a) Escapees vs Inactivees: Blood samples in controls (GSE237036)

metfemale <- t(as.matrix(met[sel, pheno$`disease state:ch1` == "normal"])) #Subset only healthy individuals.
#CpG levels of inactive genes
ig <- IGlevels(metfemale, colnames(metfemale))
#CpG levels of active genes.
eg <- EGlevels(metfemale, colnames(metfemale))

to_plotIG <- density(ig, na.rm = TRUE)
to_plotEG <- density(eg, na.rm = TRUE)

plot(to_plotIG,
     lwd = 3, ylim = c(0, 3.1),
     main = "Escapees vs Inactivees: Blood samples in controls (GSE237036)", xlab = "CpG methylation levels",
     cex.lab = 1.4, 
     cex.axis = 1.4,
     cex.main = 1.7); 
lines(to_plotEG, lwd = 3, lty = 2, col = "red");
legend("topright", legend = c("CpGs levels in inactive genes", "CpGs levels in escapees"), col = c("black", "red"), lwd = c(3, 3), lty = c(1, 2), cex = 1.4);  
polygon(c(to_plotIG$x[to_plotIG$x < 0.2], 0.2), c(to_plotIG$y[to_plotIG$x < 0.2], 0), col = rgb(0, 0, 0, alpha = 0.5))

## a) Escapees vs Inactivees: Blood samples in TNBC (GSE237036)

metfemale <- t(as.matrix(met[sel, pheno$`disease state:ch1` == "BC"])) #Subset only affected individuals.
dim(metfemale)
#CpG levels of inactive genes
ig <- IGlevels(metfemale, colnames(metfemale))
#CpG levels of active genes.
eg <- EGlevels(metfemale, colnames(metfemale))
#Compare the distributions
to_plotIG <- density(ig, na.rm = TRUE)
to_plotEG <- density(eg, na.rm = TRUE)

plot(to_plotIG,
     lwd = 3, ylim = c(0, 3.1),
     main = "Escapees vs Inactivees: Blood samples in TNBC (GSE237036)", xlab = "CpG methylation levels",
     cex.lab = 1.4, 
     cex.axis = 1.4,
     cex.main = 1.7); lines(to_plotEG, lwd = 3, lty = 2, col = "red");
legend("topright", legend = c("CpGs levels in inactive genes", "CpGs levels in escapees"), col = c("black", "red"), lwd = c(3, 3), lty = c(1, 2), cex = 1.4);  
polygon(c(to_plotIG$x[to_plotIG$x < 0.2], 0.2), c(to_plotIG$y[to_plotIG$x < 0.2], 0), col = rgb(0, 0, 0, alpha = 0.5))

## b) X-Ra comparison cancer and control in blood samples (GSE237036)

pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)),
         `disease state:ch1` = factor(`disease state:ch1`, levels = c("normal", "BC")))

ggbetweenstats(
  data = pheno, 
  x = `disease state:ch1`, 
  y = XRa, 
  bf.message = FALSE,
  results.subtitle = FALSE,
  k = 3L,
  centrality.point.args = list(size = 5, color = "darkred"),
  centrality.label.args = list(size = 4.5, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
) +
  labs(title = "X-Ra comparison cancer and control in blood samples (GSE237036)", x = "Sample type", y = "X-Ra") +
  scale_colour_discrete(name = "Cancer type", labels = c("Control subjects", "Cancer patients")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))
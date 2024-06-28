# Figure 7 - Association between X-Ra and breast cancer in breast tissue

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggstatsplot)

## a) Escapees vs Inactivees: Breast cancer control cohort (GSE101961)

path_here <- paste(getwd(), "/GSE101961", sep = "") #Where the GSE will be downloaded.

#Get the series
gse <- getGEO("GSE101961", destdir = path_here)[[1]]

#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)

#Select CpGs only in chrX and estimate the methylation distriution levels of the samples.
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ])) 

ig <- IGlevels(metfemale, colnames(metfemale))
eg <- EGlevels(metfemale, colnames(metfemale))
# Plot both escapees and inactivees
to_plotIG <- density(ig, na.rm = TRUE)
to_plotEG <- density(eg, na.rm = TRUE)
# Final plot.
plot(to_plotIG,
     lwd = 3, ylim = c(0, 3.4),
     main = "Escapees vs Inactivees: Breast cancer control cohort (GSE101961)", xlab = "CpG methylation levels",
     cex.lab = 1.4, 
     cex.axis = 1.4,
     cex.main = 1.7); 
lines(to_plotEG, lwd = 3, lty = 2, col = "red");
legend("topright", legend = c("CpGs levels in inactive genes", "CpGs levels in escapees"), col = c("black", "red"), lwd = c(3, 3), lty = c(1, 2), cex = 1.4);  
polygon(c(to_plotIG$x[to_plotIG$x < 0.2], 0.2), c(to_plotIG$y[to_plotIG$x < 0.2], 0), col = rgb(0, 0, 0, alpha = 0.5))

## a) Escapees vs Inactivees: TNBC cohort (GSE78754)

path_here <- paste(getwd(), "/GSE78754", sep = "") #Where the GSE will be downloaded.

#Get the series
gse <- getGEO("GSE78754", destdir = path_here)[[1]]

#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)

#Select CpGs only in chrX and estimate the methylation distriution levels of the samples.
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ])) 

ig <- IGlevels(metfemale, colnames(metfemale))
eg <- EGlevels(metfemale, colnames(metfemale))

to_plotIG <- density(ig, na.rm = TRUE)
to_plotEG <- density(eg, na.rm = TRUE)

plot(to_plotIG,
     lwd = 3, ylim = c(0, 3.4),
     main = "Escapees vs Inactivees: TNBC cohort (GSE78754)", xlab = "CpG methylation levels", cex.lab = 1.4, 
     cex.axis = 1.4,
     cex.main = 1.7); 
lines(to_plotEG, lwd = 3, lty = 2, col = "red");
legend("topright", legend = c("CpGs levels in inactive genes", "CpGs levels in escapees"), col = c("black", "red"), lwd = c(3, 3), lty = c(1, 2), cex = 1.4);  
polygon(c(to_plotIG$x[to_plotIG$x < 0.2], 0.2), c(to_plotIG$y[to_plotIG$x < 0.2], 0), col = rgb(0, 0, 0, alpha = 0.5));


## b) X-Ra comparison control and cancer in breast (GSE60185)

path_here <- paste(getwd(), "/GSE60185", sep = "")

#Get the series
gse <- getGEO("GSE60185", destdir = path_here)[[1]]

#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)

#Select CpGs only in chrX and estimate X-Ra.
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)), 
         `disease state:ch1` = factor(`disease state:ch1`, levels = c("normal", "breast cancer")))
         
# Final plot

ggbetweenstats(
  data = pheno, 
  x = `disease state:ch1`, 
  y = XRa, 
  bf.message = FALSE,
  results.subtitle = FALSE,
  pairwise.display = "none",
  k = 4L, 
  centrality.label.args = list(size = 4.5, nudge_x = 0.4, segment.linetype = 4,
                               min.segment.length = 0)
) + 
  labs(title = "X-Ra comparison control and cancer in breast (GSE60185)", x = "Tissue state") +
  scale_colour_manual(name = "Breast tissue state", labels = c("Normal tissue", "Cancer tissue"), values = c("brown4", "purple")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16))
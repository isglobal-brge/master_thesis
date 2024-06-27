# Figure 10 - Deconvolution results from FlowSorted.Blood.450k/EPIC R package.

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(FlowSorted.Blood.EPIC)

# Figure generated using the dataset: (GSE60185)

path_here <- paste(getwd(), "/GSE60185", sep = "")

#Get the series
gse <- getGEO("GSE60185", destdir = path_here)[[1]]

#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse)) #All females.
annotation <- fData(gse)

#Select CpGs only in chrX and estimate X-Ra.
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)), 
         `disease state:ch1` = factor(`disease state:ch1`, levels = c("normal", "breast cancer")))

#### Deconvolution algorithm - Sometimes time consumuing.
id <- "GSE60185"
betas <- exprs(gse)

#With the betas, I should be able to use the 2nd deconvolution method.

common <- intersect(rownames(betas), IDOLOptimizedCpGs450klegacy)
betas <- betas[common, ]
common <- IDOLOptimizedCpGs450klegacy %in% rownames(betas) 

propEPIC2<-projectCellType_CP (
  betas[IDOLOptimizedCpGs450klegacy[common],],
  IDOLOptimizedCpGs450klegacy.compTable, contrastWBC=NULL, nonnegative=TRUE,
  lessThanOne=FALSE)

percEPIC2<-round(propEPIC2*100,1)
write.csv(percEPIC2, file = paste("deconvolution_", id, ".csv", sep = ""), row.names = FALSE)

#### Merge the counts with the pheno
deconv_results <- read.delim(paste(path_here, "/deconvolution_GSE60185.csv", sep = ""), sep = ",")
rownames(deconv_results) <- rownames(pheno)
#Write the % in the phenotypes dataframe
pheno$CD8T <- deconv_results$CD8T
pheno$CD4T <- deconv_results$CD4T
pheno$NK <- deconv_results$NK
pheno$Bcell <- deconv_results$Bcell
pheno$Mono <- deconv_results$Mono
pheno$Neu <- deconv_results$Neu

#Function used to extract the legend object from a ggplot and store it into an independent object.
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Plots
plot1 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8T, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease status") +
  geom_smooth(formula = y ~ x, aes(x = CD8T, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  geom_smooth(aes(x = CD8T, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  labs(title = "XRa levels and CD8T counts", x = "Cell type (CD8T) counts") +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme_light() +
  annotate(geom = "text", label = "Cancer corr = -0.1\nControl corr = -0.245*", x = 15.5, y = 0.40, size = 6) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position = "bottom")

legend <- get_legend(plot1)

plot1 <- plot1 + theme(legend.position="none")

plot2 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4T, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease status") +
  geom_smooth(formula = y ~ x, aes(x = CD4T, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  geom_smooth(aes(x = CD4T, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  labs(title = "XRa levels and CD4T counts", x = "Cell type (CD4T) counts") +
  annotate(geom = "text", label = "Cancer corr = 0.42****\nControl corr = 0.3**", x = 28, y = 0.40, size = 6) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot3 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = NK, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease status") +
  geom_smooth(formula = y ~ x, aes(x = NK, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  geom_smooth(aes(x = NK, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  annotate(geom = "text", label = "Cancer corr = 0.174****\nControl corr = 0.083", x = 20, y = 0.40, size = 6) +
  labs(title = "XRa levels and NK counts", x = "Cell type (NK) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot4 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bcell, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease Status") +
  geom_smooth(formula = y ~ x, aes(x = Bcell, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  annotate(geom = "text", label = "Cancer corr = -0.453****\nControl corr = -0.3**", x = 14.5, y = 0.40, size = 6) +
  geom_smooth(aes(x = Bcell, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  labs(title = "XRa levels and Bcell counts", x = "Cell type (Bcell) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none")

plot5 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Mono, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease status") +
  geom_smooth(formula = y ~ x, aes(x = Mono, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  geom_smooth(aes(x = Mono, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  annotate(geom = "text", label = "Cancer corr = 0.1\nControl corr = 0.14", x = 18, y = 0.40, size = 6) +
  labs(title = "XRa levels and Mono counts", x = "Cell type (Mono) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot6 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Neu, col = DiseaseState), size = 2) +
  scale_colour_discrete(name = "Disease status") +
  geom_smooth(formula = y ~ x, aes(x = Neu, y = XRa), method = "lm", col = "black", size = 2, lty = "dashed", se = F) +
  annotate(geom = "text", label = "Cancer corr = -0.055\nControl corr = 0.05", x = 25, y = 0.40, size = 6) +
  geom_smooth(aes(x = Neu, y = XRa, col = DiseaseState, fill = DiseaseState), formula = y ~ x, size = 2, method = "lm", se = T) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  labs(title = "XRa levels and Neu counts", x = "Cell type (Neu) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
library(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, blankPlot, legend, blankPlot, ncol=3, nrow = 3, widths = c(5, 5, 5), heights = c(5, 5, 1))

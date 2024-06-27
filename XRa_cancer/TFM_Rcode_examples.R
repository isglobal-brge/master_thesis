########
#Code example1: GSE184159
########

#Load packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggstatsplot)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(factoextra)
library(stats)
library(ggplot2)
library(ggfortify)
library(survival)
library(survminer)
library(minfi)
library(GenomicRanges)
library(DMRcate)
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(irr)
library(lmerTest)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
library(BiocParallel)
library(ExperimentHub)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(data.table)
#Get the series
gse <- getGEO("GSE184159", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse) #No methylation data

#Phenotypes
pheno <- pData(phenoData(gse))

annotation <- fData(gse) #There is NO annotation, we have to get it from somewhere else.

annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
dim(annotation)
#Methylation, downloaded from GEO entry manually. Already the beta values
#wget “https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184159/suppl/GSE184159_beta_detp.csv.gz” -O GSE184159_beta_detp.csv #Download the processed betas using Linux command line.
met <- read.csv("GSE184159/GSE184159_beta_detp.csv")
#Process a bit the methylation data
rownames(met) <- met$X #add rownames and eliminate X column
met <- met[, grep("^beta", colnames(met))]
colnames(met) <- pheno$geo_accession
#We have to reorder the rows
dim(met)
identical(rownames(met), rownames(annotation))
same_rows <- intersect(rownames(met), rownames(annotation))
#Subset
met <- met[same_rows, ]
annotation <- annotation[same_rows, ]
identical(rownames(met), rownames(annotation)) #Now we can use the dataset.
#We write all the tissue names to identify which are the same ones.
pheno$IDs <- gsub("_.*_.", "", pheno$title)
#Select CpGs only in chrX.
sel <- as.vector(annotation$chr) == "chrX"
metfemale <- t(as.matrix(met[sel, ]))
#X-Ra estimation
pheno$XRa <- XRa(metfemale, colnames(metfemale))

#Distribution plots: Escapees vs Inactivees
ig <- IGlevels(metfemale, colnames(metfemale))
eg <- EGlevels(metfemale, colnames(metfemale))
to_plotIG <- density(ig, na.rm = TRUE)
to_plotEG <- density(eg, na.rm = TRUE)
plot(to_plotIG,
     lwd = 3, ylim = c(0, 3.1),
     main = "X-Ra Distribution comparison (Breast tissue)", xlab = "CpG methylation levels"); lines(to_plotEG, lwd = 3, lty = 2, col = "red");legend("topright", legend = c("CpGs levels in inactive genes", "CpGs levels in escapees"), col = c("black", "red"), lwd = c(3, 3), lty = c(1, 2)); abline(v = 0.2, col = "blue"); polygon(c(to_plotIG$x[to_plotIG$x < 0.2], 0.2), c(to_plotIG$y[to_plotIG$x < 0.2], 0), col = rgb(0, 0, 0, alpha = 0.5))
#XRa dynamics plots.
pheno$time_point <- pheno$`timepoint:ch1`
pheno$response <- pheno$`response:ch1`
pheno_cancer <- pheno[pheno$source_name_ch1 == "Breast_Cancer", ]
ggpaired(data = pheno, 
         x = "time_point", 
         y = "XRa", 
         id = "IDs", 
         color = "black", 
         line.color = "IDs", 
         legend = "none") +
  labs(title = "X-Ra comparison between treatment timepoints (GSE184159)", x = "Timepoint", y = "X-Ra") +
  stat_compare_means(paired = TRUE) +
  scale_x_discrete(labels = c("Control", "Before treatment", "Midpoint", "After treatment")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 14), legend.position = "None")

pheno_cancer$num_time <- ifelse(pheno_cancer$`timepoint:ch1` == "A", 1, ifelse(pheno_cancer$`timepoint:ch1` == "B", 2, 3))

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
           label.y.npc = "top") +
  geom_jitter(data = pheno_cancer, aes(x = num_time, y = XRa, col = response)) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Before treatment", "Midpoint", "After treatment")) +
  labs(title = "Response comparison depending on X-Ra in Breast cancer (GSE184159)", x = "Timepoint", y = "X-Ra") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 14), legend.text = element_text(size = 10), legend.title = element_text(size = 12))

#Test statistical significance between conditions
dunn.test::dunn.test(pheno_cancer$XRa, pheno_cancer$`timepoint:ch1`, method = "Bonferroni")

#Linear models comparing two timepoints (A and B)
pheno$time <- ifelse(pheno$`timepoint:ch1` == "N", 0, ifelse(pheno$`timepoint:ch1` == "A", 1, ifelse(pheno$`timepoint:ch1` == "B", 2, 3)))
#We get only those with two observations
opheno <- pheno[(pheno$time == 1 | pheno$time == 2) & !pheno$IDs %in% c("CC046", "DB125", "EH093", "JR058", "LK051", "PR068", "SH196", "SH201", "SS066", "VP035"), ]

data <- matrix(opheno$XRa, 
               ncol = 2, byrow = TRUE)
colnames(data) <- c("Timepoint A", "Timepoint B")
rownames(data) <- unique(opheno$IDs)
icc(data, model = "twoway", type = "consistency", unit = "single") #ICC estimation

#Fit the linear model.
opheno$num_resp <- ifelse(opheno$`response:ch1` == "Non", 1, ifelse(opheno$`response:ch1` == "Partial", 2, 3))
opheno$num_ID <- rep(seq(1, nrow(opheno)/2), each = 2)
mod <- lmerTest::lmer(XRa ~ time + num_resp + (1 | num_ID), 
                      data = opheno)
summary(mod)


#Deconvolution: FlowSorted.Blood.EPIC
id <- "GSE184159"
betas <- read.csv("GSE184159_beta_detp.csv")
rownames(betas) <- betas$X
betas1 <- betas[, grep("^beta", colnames(betas))]
#Subset for the common probes.

common <- intersect(rownames(betas), IDOLOptimizedCpGs) #The array is an EPIC chip.
betas1 <- betas1[common, ]
common <- IDOLOptimizedCpGs %in% rownames(betas1) 
propEPIC2<-projectCellType_CP (
  betas1[IDOLOptimizedCpGs[common],],
  IDOLOptimizedCpGs.compTable, contrastWBC=NULL, nonnegative=TRUE,
  lessThanOne=FALSE)


percEPIC2<-round(propEPIC2*100,1)
write.csv(percEPIC2, file = paste("deconvolution_", id, ".csv", sep = ""), row.names = FALSE)
#Write in each individual the counts.
deconv_results <- read.delim(paste(path_here, "/deconvolution_GSE184159.csv", sep = ""), sep = ",")
rownames(deconv_results) <- rownames(pheno)
#Write the % in the phenotypes dataframe
pheno$CD8T <- deconv_results$CD8T
pheno$CD4T <- deconv_results$CD4T
pheno$NK <- deconv_results$NK
pheno$Bcell <- deconv_results$Bcell
pheno$Mono <- deconv_results$Mono
pheno$Neu <- deconv_results$Neu
#Deconvolution regression plots
pheno$Response <- factor(pheno$`response:ch1`, levels = c("Complete", "Partial", "Non"))
pheno$TimePoint <- factor(pheno$`timepoint:ch1`, levels = c("N", "A", "B", "C"))
plot1 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8T, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD8T, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD8T counts", x = "Cell type (CD8T) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot2 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4T, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD4T, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD4T counts", x = "Cell type (CD4T) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot3 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = NK, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = NK, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and NK counts", x = "Cell type (NK) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot4 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bcell, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Bcell, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Bcell counts", x = "Cell type (Bcell) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot5 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Mono, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Mono, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Mono counts", x = "Cell type (Mono) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot6 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Neu, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Neu, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Neu counts", x = "Cell type (Neu) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)
#Extra plots
new_pheno <- pheno[, c("IDs", "TimePoint", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")] %>%
  pivot_longer(cols = -c(IDs, TimePoint), 
               names_to = c("TipoCelular", "Medida"),
               names_sep = "_")
ggplot(data = new_pheno) +
  geom_boxplot(aes(x = TipoCelular, y = value, fill = factor(TimePoint))) +
  labs(title = "TILs % counts between timepoints", x = "Cell type percentage", y = "% counts") +
  theme_light() +
  scale_fill_manual(values = c("green", "purple", "blue", "yellow"), name = "Timepoint") +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  guides(shape = guide_legend(title = "Cell type", override.aes = list(col = "black")))

#Compute PCA
pca_res <- prcomp(pheno[, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")])
cor(pca_res$x, pheno$XRa)
fviz_eig(pca_res, addlabels = TRUE, title = "Eigenvalues distribution")
plot1 <- fviz_pca_var(pca_res, 
                      axes = c(1, 3),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE,     # Avoid text overlapping
                      title = "PCA - Cell counts change depending on the cell type"
)
pheno$HighXRa <- ifelse(pheno$XRa > 0.11, "1", "0")
plot2 <- autoplot(pca_res, x = 1, y = 3, data = pheno, color = "HighXRa", shape = "Response", size = 3) +
  theme_light() +
  scale_shape_manual(name = "Response",
                     values = c(15, 18, 19), 
                     labels = c("Complete", "Partial", "Non")) +
  scale_colour_manual(name = "XRa levels", 
                      labels = c("XRa < 0.11", "XRa > 0.11"), 
                      values = c("green", "purple")) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "PCA with XRa tendency depending on the timepoint") +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot_grid(plot1, plot2)
#Correlation between X-Ra and immune counts
shapiro.test(pheno$XRa); shapiro.test(pheno$CD8T); shapiro.test(pheno$CD4T); shapiro.test(pheno$NK); shapiro.test(pheno$Bcell); shapiro.test(pheno$Mono); shapiro.test(pheno$Neu) #Normality test
#Test CD8T
cor.test(pheno$CD8T, pheno$XRa, method = "spearman")
#Test CD4T
cor.test(pheno$CD4T, pheno$XRa, method = "spearman")
#Test NK
cor.test(pheno$NK, pheno$XRa, method = "spearman")
#Test Bcells
cor.test(pheno$Bcell, pheno$XRa, method = "spearman")
#Test Mono
cor.test(pheno$Mono, pheno$XRa, method = "spearman")
#Test Neu
cor.test(pheno$Neu, pheno$XRa, method = "spearman")

#Deconvolution method: FlowSorted.BloodExtended.EPIC
#wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184159&format=file&file=GSE184159%5Fraw%5Fdetp%2Ecsv%2Egz" -O GSE184159_raw_detp.csv.gz #Download the beta values using Linux command line (no IDATs available)
#Build the EPIC array to use the algorithm
ID_entry <- "GSE184159"
path_entry <- paste("/home/isglobal.lan/aalegret/GEOthings/FlowSorted_TESTING/", ID_entry, sep = "")

general_data <- fread(paste(path_entry, "/GSE184159_raw_detp.csv.gz", sep = ""), data.table = "FALSE")
rownames(general_data) <- general_data$V1
meth_data <- general_data[, grep("^meth_", colnames(general_data))]
colnames(meth_data) <- gsub("^meth_", "", colnames(meth_data))
unmeth_data <- general_data[, grep("^unmeth_", colnames(general_data))]
colnames(unmeth_data) <- gsub("^unmeth_", "", colnames(unmeth_data))
MSet <- MethylSet(Meth = as.matrix(meth_data), Unmeth = as.matrix(unmeth_data))
MSet@annotation["array"]="IlluminaHumanMethylationEPIC"
MSet@annotation["annotation"] = "ilm10b2.hg19"

hub <- ExperimentHub()
query(hub, "FlowSorted.BloodExtended.EPIC")
libraryDataGet <- function(title) {
  assign(title, ExperimentHub()[[query(
    ExperimentHub(),
    title
  )$ah_id]])
}
FlowSorted.BloodExtended.EPIC <-libraryDataGet("FlowSorted.BloodExtended.EPIC")
deconv_name <- paste("deconv_results_", ID_entry, ".txt", sep = "")
if (MSet@annotation["array"]=="IlluminaHumanMethylationEPIC"){
  ############
  #EPIC
  ############
  countsEPIC<-estimateCellCounts2(MSet, compositeCellType = "BloodExtended",
                                  processMethod = "preprocessQuantile",
                                  probeSelect = "IDOL",
                                  cellTypes = c("Bas", "Bmem", "Bnv",
                                                "CD4mem", "CD4nv",
                                                "CD8mem", "CD8nv", "Eos",
                                                "Mono", "Neu", "NK", "Treg"),
                                  referencePlatform =
                                    "IlluminaHumanMethylationEPIC",
                                  referenceset = NULL,
                                  IDOLOptimizedCpGs =IDOLOptimizedCpGs,
                                  returnAll = FALSE)
  head(countsEPIC$counts)
  countsEPIC.db<-as.data.frame(countsEPIC$counts)
  write.table(countsEPIC.db, deconv_name, sep="\t")
} else{
  ############
  #450K
  ############
  counts450K<-estimateCellCounts2(MSet, compositeCellType = "BloodExtended",
                                  processMethod = "preprocessQuantile",
                                  probeSelect = "IDOL",
                                  cellTypes = c("Bas", "Bmem", "Bnv",
                                                "CD4mem", "CD4nv",
                                                "CD8mem", "CD8nv", "Eos",
                                                "Mono", "Neu", "NK", "Treg"),
                                  referencePlatform =
                                    "IlluminaHumanMethylationEPIC",
                                  referenceset = NULL,
                                  IDOLOptimizedCpGs =IDOLOptimizedCpGs450klegacy,
                                  returnAll = FALSE)
  head(counts450K$counts)
  counts450K.db<-as.data.frame(counts450K$counts)
  write.table(counts450K.db, deconv_name, sep="\t")
}

#Write the observations of each subject
good_deconv_results <- read.delim(paste(path_here, "/deconv_results_GSE184159.txt", sep = ""), sep = "\t")
rownames(pheno) <- sub("_.*_", "_", pheno$title)
#Write the % in the phenotypes dataframe
pheno$Bas <- good_deconv_results$Bas
pheno$Bmem <- good_deconv_results$Bmem
pheno$Bnv <- good_deconv_results$Bnv
pheno$CD4mem <- good_deconv_results$CD4mem
pheno$CD4nv <- good_deconv_results$CD4nv
pheno$CD8mem <- good_deconv_results$CD8mem
pheno$CD8nv <- good_deconv_results$CD8nv
pheno$Eos <- good_deconv_results$Eos
pheno$Mono <- good_deconv_results$Mono
pheno$Neu <- good_deconv_results$Neu
pheno$NK <- good_deconv_results$NK
pheno$Treg <- good_deconv_results$Treg

#Correlation plots
plot1 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bas, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Bas, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Bas counts", x = "Cell type (Bas) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot2 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bmem, col = Response, shape = TimePoint), size = 3) +
  labs(title = "XRa levels and Bmem counts", x = "Cell type (Bmem) counts") +
  geom_smooth(formula = y ~ x, aes(x = Bmem, y = XRa), method = "lm", col = "black", size = 2) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot3 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bnv, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Bnv, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Bnv counts", x = "Cell type (Bnv) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot4 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4mem, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD4mem, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD4mem counts", x = "Cell type (CD4mem) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot5 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4nv, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD4nv, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD4nv counts", x = "Cell type (CD4nv) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot6 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8mem, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD8mem, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD8mem counts", x = "Cell type (CD8mem) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot7 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8nv, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD8nv, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and CD8nv counts", x = "Cell type (CD8nv) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot8 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Eos, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Eos, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Eos counts", x = "Cell type (Eos) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot9 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Mono, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Mono, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Mono counts", x = "Cell type (Mono) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot10 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Neu, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Neu, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Neu counts", x = "Cell type (Neu) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot11 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = NK, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = NK, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and NK counts", x = "Cell type (NK) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot12 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Treg, col = Response, shape = TimePoint), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Treg, y = XRa), method = "lm", col = "black", size = 2) +
  labs(title = "XRa levels and Treg counts", x = "Cell type (Treg) counts") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot_grid(plot1, plot2, plot3, plot4, plot5, plot6)
plot_grid(plot7, plot8, plot9, plot10, plot11, plot12)
#Extra plots
new_pheno <- pheno[, c("IDs", "TimePoint", "Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")] %>%
  pivot_longer(cols = -c(IDs, TimePoint ), 
               names_to = c("TipoCelular", "Medida"),
               names_sep = "_")
ggplot(data = new_pheno) +
  geom_boxplot(aes(x = TipoCelular, y = value, fill = factor(TimePoint))) +
  labs(title = "TILs % counts between timepoints", x = "Cell type percentage", y = "% counts") +
  theme_light() +
  scale_fill_manual(values = c("green", "purple", "blue", "yellow"), name = "Timepoint") +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  guides(shape = guide_legend(title = "Cell type", override.aes = list(col = "black")))
#Compute PCA
pca_res <- prcomp(pheno[, c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")])
cor(pca_res$x, pheno$XRa)
fviz_eig(pca_res, addlabels = TRUE, title = "Eigenvalues distribution")
#Visualization

plot1 <- fviz_pca_var(pca_res,
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE,     # Avoid text overlapping
                      title = "PCA - Cell counts change depending on the cell type"
)

pheno$HighXRa <- ifelse(pheno$XRa > 0.11, "1", "0") #Put a threshold to try to detect any clustering
plot2 <- autoplot(pca_res, data = pheno, color = "HighXRa",shape = "Response", size = 3, alpha = 0.7) +
  theme_light() +
  scale_colour_manual(name = "XRa levels", 
                      labels = c("XRa < 0.11", "XRa > 0.11"), 
                      values = c("green", "purple")) +
  scale_shape_discrete(name = "Response type") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "PCA with XRa tendency depending on the treatment timepoint") +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 12))
plot_grid(plot1, plot2)
#Estiamte correlations.
shapiro.test(pheno$XRa) #Normality test
#Test Bas
cor.test(pheno$Bas, pheno$XRa, method = "spearman")
#Test Bmem
cor.test(pheno$Bmem, pheno$XRa, method = "spearman")
#Test Bnv
cor.test(pheno$Bnv, pheno$XRa, method = "spearman")
#Test CD4mem
cor.test(pheno$CD4mem, pheno$XRa, method = "spearman")
#Test CD4nv
cor.test(pheno$CD4nv, pheno$XRa, method = "spearman")
#Test CD8mem
cor.test(pheno$CD8mem, pheno$XRa, method = "spearman")
#Test CD8nv
cor.test(pheno$CD8nv, pheno$XRa, method = "spearman")
#Test Eos
cor.test(pheno$Eos, pheno$XRa, method = "spearman")
#Test Mono
cor.test(pheno$Mono, pheno$XRa, method = "spearman")
#Test Neu
cor.test(pheno$Neu, pheno$XRa, method = "spearman")
#Test NK
cor.test(pheno$NK, pheno$XRa, method = "spearman")
#Test Treg
cor.test(pheno$Treg, pheno$XRa, method = "spearman")




########### 
#Code example 2: dataset GSE60185. Estimation of DMRs
###########


gse <- getGEO("GSE60185", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)
identical(rownames(met), rownames(annotation)) #Same order.
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
#X-Ra estimation
pheno$XRa <- XRa(metfemale, colnames(metfemale))

#DMRs using X-Ra to fit the model.
cancer_ind <- pheno$`disease state:ch1` == "breast cancer"
pheno_cancer <- pheno[cancer_ind, ]
met_cancer <- met[, cancer_ind]
#Design matrix
design <- model.matrix(~pheno_cancer$XRa)
#Eliminate those CpGs erroneously annotated.
not_NA <- !is.na(annotation$CHR) & !is.na(annotation$MAPINFO)
annotation <- annotation[not_NA, ]
met_cancer <- met[not_NA, ]
identical(rownames(met_cancer), rownames(annotation))

myAnnotation <- cpg.annotate(object = met_cancer, datatype = "array", what = "Beta", 
                             arraytype = c("450K"), 
                             analysis.type = "differential", design = design, 
                             coef = 2)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
filtered_ranges <- results.ranges[abs(results.ranges$maxdiff) > 0.02, ]

#GO enrichment
gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="GO", array.type="450K")
top_results_all <- topGSA(gst.region, n = 300)
top_results <- topGSA(gst.region[gst.region$ONTOLOGY == "BP", ], n = 10)
top_results$TERM <- factor(top_results$TERM, levels = rev(top_results$TERM))
top_results <- as_tibble(top_results)

ggplot(data = top_results, mapping = aes(y = -log10(FDR), x = TERM)) +
  geom_point(aes(col = -log10(P.DE), size = DE)) +
  geom_segment(aes(x = TERM, y = 0, 
                   xend = TERM, yend = -log10(FDR)), 
               linetype = "dashed", linewidth = 0.01, col = "grey5") +
  scale_colour_gradient(low = "cyan", high = "purple3", name = "-log10(P-val)") +
  scale_size_continuous(name = "Diff met genes") +
  labs(title = "DMRcate: GO enrichment DMRs analysis", y = "-log10(FDR)", x = "GO biological processes terms", caption = "All breast cancer samples, continuous XRa (top 10 results)") +
  theme_light() +
    theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()
#KEGG enrichment
gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="KEGG", array.type="450K")
top_results_all <- topGSA(gst.region, n = 10)
top_results$Description <- factor(top_results$Description, levels = rev(top_results$Description))
top_results <- as_tibble(top_results)

ggplot(data = top_results, mapping = aes(y = -log10(FDR), x = Description)) +
  geom_point(aes(col = -log10(P.DE), size = DE)) +
  geom_segment(aes(x = Description, y = 0, 
                   xend = Description, yend = -log10(FDR)), 
               linetype = "dashed", linewidth = 0.01, col = "grey5") +
  scale_colour_gradient(low = "cyan", high = "purple3", name = "-log10(P-val)") +
  scale_size_continuous(name = "Diff met genes") +
  labs(title = "DMRcate: KEGG enrichment DMRs analysis", y = "-log10(FDR)", x = "KEGG terms", caption = "All breast cancer samples, continuous XRa") +
  theme_light() +
    theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()

#Fit the model using Case/Control
#Design matrix
design <- model.matrix(~pheno$`disease state:ch1`)
#Eliminate those CpGs erroneously annotated.
not_NA <- !is.na(annotation$CHR) & !is.na(annotation$MAPINFO)
annotation <- annotation[not_NA, ]
met_cancer <- met[not_NA, ]
identical(rownames(met_cancer), rownames(annotation))
#Model
myAnnotation <- cpg.annotate(object = met_cancer, datatype = "array", what = "Beta", 
                             arraytype = c("450K"), 
                             analysis.type = "differential", design = design, 
                             coef = 2)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)

results.ranges <- extractRanges(DMRs)
filtered_ranges <- results.ranges[abs(results.ranges$maxdiff) > 0.02, ]

#GO enrichment
gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="GO", array.type="450K")
top_results_all <- topGSA(gst.region, n = 300)
top_results <- topGSA(gst.region[gst.region$ONTOLOGY == "BP", ], n = 10)
top_results$TERM <- factor(top_results$TERM, levels = rev(top_results$TERM))
top_results <- as_tibble(top_results)
ggplot(data = top_results, mapping = aes(y = -log10(FDR), x = TERM)) +
  geom_point(aes(col = -log10(P.DE), size = DE)) +
  geom_segment(aes(x = TERM, y = 0, 
                   xend = TERM, yend = -log10(FDR)), 
               linetype = "dashed", linewidth = 0.01, col = "grey5") +
  scale_colour_gradient(low = "cyan", high = "purple3", name = "-log10(P-val)") +
  scale_size_continuous(name = "Diff met genes") +
  labs(title = "DMRcate: GO enrichment DMRs analysis", y = "-log10(FDR)", x = "GO biological processes terms", caption = "Case/Control comparison") +
  theme_light() +
    theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()
#KEGG enrichment
gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="KEGG", array.type="450K")
top_results <- topGSA(gst.region, n = 10)

top_results$Description <- factor(top_results$Description, levels = rev(top_results$Description))
top_results <- as_tibble(top_results)

ggplot(data = top_results, mapping = aes(y = -log10(FDR), x = Description)) +
  geom_point(aes(col = -log10(P.DE), size = DE)) +
  geom_segment(aes(x = Description, y = 0, 
                   xend = Description, yend = -log10(FDR)), 
               linetype = "dashed", linewidth = 0.01, col = "grey5") +
  scale_colour_gradient(low = "cyan", high = "purple3", name = "-log10(P-val)") +
  scale_size_continuous(name = "Diff met genes") +
  labs(title = "DMRcate: KEGG enrichment DMRs analysis", y = "-log10(FDR)", x = "KEGG terms", caption = "Case/Control comparison") +
  theme_light() +
    theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()

########### 
#Code example 3: Survival analysis
###########

gse <- getGEO("GSE78754", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)
identical(rownames(met), rownames(annotation))
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ])) 
#Estimate X-Ra
pheno$XRa <- XRa(metfemale, colnames(metfemale))
#Survival analysis

#Dataset
pheno <- as_tibble(pheno)
pheno$status <- ifelse(pheno$`survival:ch1` == "alive", 1, 2) #Alive patients are censored (1), dead patients are informative, sufferend the event (2)
pheno$HighXRa <- ifelse(pheno$XRa > 0.11, "1", "0")
pheno$`month of follow up/to death:ch1` <- as.double(pheno$`month of follow up/to death:ch1`)
#Fit the survival model – Kaplan-Meier curve. We group according to different XRa levels (high and low)
sfit <- survfit(Surv(pheno$`month of follow up/to death:ch1`, pheno$status)~pheno$HighXRa, data = pheno)
sfit
summary(sfit)
#If we want to display the results in a more compact manner
summary(sfit, times = seq(0, 180, 20)) #Write the range according to our times interval (the month of follow up)

ggsurvplot(sfit, 
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("XRa < 0.11", "XRa > 0.11"),
           legend.title = "XRa values",
           title = "Kaplan-Meier curve for TNBC survival depending on X-Ra values",
           ncensor.plot = TRUE)
#Fit the survival model – CPHM. Continuous X-Ra
fit <- coxph(Surv(pheno$`month of follow up/to death:ch1`, pheno$status)~pheno$XRa , data = pheno)
summary(fit)
test.ph.a <- cox.zph(fit)
plot(test.ph.a)

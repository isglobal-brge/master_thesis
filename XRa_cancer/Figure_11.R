# Figure 11 - Deconvolution results from FlowSorted.BloodExtended.EPIC R package

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(minfi)
library(FlowSorted.Blood.EPIC)
library(FlowSorted.BloodExtended.EPIC)
library(BiocParallel)
library(ExperimentHub)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)

# Figure generated using the dataset: (GSE60185)

############################## 
# Deconvolution algorithm. Specific requirements: R versions 4.1.2

# Download the IDATs manually from GEO or using Bash:
# wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60185&format=file" -O GSE60185_RAW.tar
# tar -xvf GSE60185_RAW.tar

library(GEOquery)
#GEO dataset
path_here <- paste(getwd(), "/GSE60185/", sep = "")
ID_entry <- "GSE60185"
path_entry <- paste(path_here, ID_entry, sep = "")
gse <- getGEO(ID_entry, destdir = path_entry)[[1]]
pheno <- pData(phenoData(gse))

# Generate the SampleSheet for minfi package.

file_names <- list.files(path_entry)
sentrix_ID <- c()
sentrix_Position <- c()
GEO_entry <- c()
for (i in file_names[grep("_Grn.idat.gz$", file_names)]) {
  sentrix_ID <- c(sentrix_ID, sub(".*_(\\d+)_R\\d+C\\d+_.*", "\\1", i))
  sentrix_Position <- c(sentrix_Position, sub(".*_\\d+_(R\\d+C\\d+)_.*", "\\1", i))
  GEO_entry <- c(GEO_entry, sub("_.*", "", i))
}
idat_names <- data.frame(array = sentrix_Position, slide = sentrix_ID)
rownames(idat_names) <- GEO_entry 
#Important to test if the phenotypes and the idats are in the same order.
identical(rownames(pheno), rownames(idat_names)) #Both rownames are supposed to be the GEO accession, change the code if it is not that.

pheno$Array <- idat_names$array
pheno$Slide <- idat_names$slide

#Write the samplesheet with the phenotypes
write.csv(pheno, file = "sample_sheet.csv", row.names = FALSE)


targets <- read.metharray.sheet(path_here) #metharray.sheet reads the sample.sheet
RGSet <- read.metharray.exp(targets = targets) #Reads a metharray experiment previously stored using the funciton above.

hub <- ExperimentHub()
query(hub, "FlowSorted.BloodExtended.EPIC")
libraryDataGet <- function(title) {
  assign(title, ExperimentHub()[[query(
    ExperimentHub(),
    title
  )$ah_id]])
}
FlowSorted.BloodExtended.EPIC <-libraryDataGet("FlowSorted.BloodExtended.EPIC")
# Obtain the counts.
deconv_name <- paste("deconv_results_", ID_entry, ".txt", sep = "")
if (RGSet@annotation["array"]=="IlluminaHumanMethylationEPIC"){
  ############
  #EPIC
  ############
  countsEPIC<-estimateCellCounts2(RGSet, compositeCellType = "BloodExtended",
                                  processMethod = "preprocessNoob",
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
  counts450K<-estimateCellCounts2(RGSet, compositeCellType = "BloodExtended",
                                  processMethod = "preprocessNoob",
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

##################################

# Turn to R version 4.3

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

# Load the counts and mix them with the pheno
good_deconv_results <- read.delim(paste(path_here, "/deconv_results_GSE60185.txt", sep = ""), sep = "\t")
rownames(good_deconv_results) <- sub("_.*_.*", "", rownames(good_deconv_results))
#All are in the correct order.
#pheno <- pheno2
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


#Function used to extract the legend object from a ggplot and store it into an independent object.
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Plot

plot1 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bas, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = Bas, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Bas, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  labs(title = "XRa levels and Bas counts", x = "Cell type (Bas) counts") +
  annotate(geom = "text", label = "Cancer corr = 0.136**\nControl corr = 0.005", x = 0.12, y = 0.40, size = 6) +
  scale_colour_discrete(name = "Disease status") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position = "bottom")

legend <- get_legend(plot1)

plot1 <- plot1 + theme(legend.position="none")

plot2 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bmem, col = DiseaseState), size = 3, show.legend = FALSE) +
  labs(title = "XRa levels and Bmem counts", x = "Cell type (Bmem) counts") +
  geom_smooth(formula = y ~ x, aes(x = Bmem, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Bmem, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = -0.06\nControl corr = -0.04", x = 0.2, y = 0.40, size = 6) +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot3 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Bnv, col = DiseaseState), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Bnv, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Bnv, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = -0.232****\nControl corr = 0.143", x = 0.08, y = 0.40, size = 6) +
  labs(title = "XRa levels and Bnv counts", x = "Cell type (Bnv) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot4 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4mem, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = CD4mem, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = CD4mem, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = -0.05\nControl corr = -0.3521**", x = 0.22, y = 0.40, size = 6) +
  labs(title = "XRa levels and CD4mem counts", x = "Cell type (CD4mem) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none")

plot5 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD4nv, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = CD4nv, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = CD4nv, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = -0.06\nControl corr = 0.05", x = 0.028, y = 0.40, size = 6) +
  labs(title = "XRa levels and CD4nv counts", x = "Cell type (CD4nv) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  scale_x_continuous(breaks = c(0.00, 0.01, 0.02, 0.03)) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot6 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8mem, col = DiseaseState), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = CD8mem, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = CD8mem, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = 0.009\nControl corr = 0.11", x = 0.15, y = 0.40, size = 6) +
  labs(title = "XRa levels and CD8mem counts", x = "Cell type (CD8mem) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot7 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = CD8nv, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = CD8nv, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = CD8nv, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  labs(title = "XRa levels and CD8nv counts", x = "Cell type (CD8nv) counts") +
  annotate(geom = "text", label = "Cancer corr = -0.135****\nControl corr = -0.14", x = 0.037, y = 0.40, size = 6) +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none")

plot8 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Eos, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = Eos, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Eos, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = 0\nControl corr = 0", x = 8e-18, y = 0.40, size = 6) +
  labs(title = "XRa levels and Eos counts", x = "Cell type (Eos) counts") +
  theme_light() +
  scale_x_continuous(breaks = c(2e-18, 2e-17)) +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot9 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Mono, col = DiseaseState), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Mono, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Mono, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = -0.168***\nControl corr = -0.06", x = 0.25, y = 0.40, size = 6) +
  labs(title = "XRa levels and Mono counts", x = "Cell type (Mono) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot10 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Neu, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = Neu, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Neu, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = 0.11*\nControl corr = -0.15", x = 0.25, y = 0.40, size = 6) +
  labs(title = "XRa levels and Neu counts", x = "Cell type (Neu) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none")

plot11 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = NK, col = DiseaseState), size = 3, show.legend = FALSE) +
  geom_smooth(formula = y ~ x, aes(x = NK, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = NK, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = 0.095\nControl corr = 0.345**", x = 0.09, y = 0.40, size = 6) +
  labs(title = "XRa levels and NK counts", x = "Cell type (NK) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

plot12 <- ggplot(data = pheno, mapping = aes(y = XRa)) +
  geom_point(aes(x = Treg, col = DiseaseState), size = 3) +
  geom_smooth(formula = y ~ x, aes(x = Treg, y = XRa), method = "lm", col = "black", size = 2, linetype = "dashed", se = F) +
  geom_smooth(formula = y ~ x, aes(x = Treg, y = XRa, col = DiseaseState, fill = DiseaseState), method = "lm", size = 2, se = T) +
  annotate(geom = "text", label = "Cancer corr = 0.2358****\nControl corr = 0.41***", x = 0.25, y = 0.40, size = 6) +
  labs(title = "XRa levels and Treg counts", x = "Cell type (Treg) counts") +
  theme_light() +
  scale_fill_manual(values = c("brown4", "purple")) +
  scale_colour_manual(values = c("brown4", "purple")) +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 18), title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22), legend.position="none", axis.title.y = element_text(color = "white"))

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, blankPlot, blankPlot, blankPlot, ncol=3, nrow = 3, widths = c(5, 5, 5), heights = c(5, 5, 1))
grid.arrange(plot7, plot8, plot9, plot10, plot11, plot12, blankPlot, legend, blankPlot, ncol=3, nrow = 3, widths = c(5, 5, 5), heights = c(5, 5, 1))
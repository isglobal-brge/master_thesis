# Supplementary Figure 7 - Results from the DMRs analysis Case/Control.

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(DMRcate)
library(missMethyl)

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


# Use all the cohort this time.
#Design matrix
design <- model.matrix(~pheno$`disease state:ch1` )
#Eliminate those CpGs erroneously annotated.
not_NA <- !is.na(annotation$CHR) & !is.na(annotation$MAPINFO)
annotation <- annotation[not_NA, ]
met <- met[not_NA, ]
identical(rownames(met), rownames(annotation))

myAnnotation <- cpg.annotate(object = met, datatype = "array", what = "Beta", 
                             arraytype = c("450K"), 
                             analysis.type = "differential", design = design, 
                             coef = 2)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)

results.ranges <- extractRanges(DMRs)

filtered_ranges <- results.ranges[abs(results.ranges$maxdiff) > 0.02, ]

# a) GO enrichment plot.
gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="GO", array.type="450K")

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
  labs(title = "DMRcate: GO enrichment DMRs analysis", y = "-log10(FDR)", x = "GO biological processes terms", caption = "Compared groups: Cancer/Control") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()

# b) KEGG enrichment plot.

gst.region <- goregion(filtered_ranges, all.cpg=rownames(met), 
                       collection="KEGG", array.type="450K")
top_results_all <- topGSA(gst.region, n = 10)
top_results$Description <- factor(top_results$Description, levels = rev(top_results$Description))
top_results <- as_tibble(top_results)

top_results$Description <- factor(top_results$Description, levels = rev(top_results$Description))
top_results <- as_tibble(top_results)

ggplot(data = top_results, mapping = aes(y = -log10(FDR), x = Description)) +
  geom_point(aes(col = -log10(P.DE), size = DE)) +
  geom_segment(aes(x = Description, y = 0, 
                   xend = Description, yend = -log10(FDR)), 
               linetype = "dashed", linewidth = 0.01, col = "grey5") +
  scale_colour_gradient(low = "cyan", high = "purple3", name = "-log10(P-val)") +
  scale_size_continuous(name = "Diff met genes") +
  labs(title = "DMRcate: KEGG enrichment DMRs analysis", y = "-log10(FDR)", x = "KEGG biological processes terms", caption = "Compared groups: Cancer/Control") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12)) +
  coord_flip()
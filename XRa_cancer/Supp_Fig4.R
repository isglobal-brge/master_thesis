# Supplementary Figure 4 - X-Ra values compared to the menopausal stage and Age.

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(ggstatsplot)

path_here <- paste(getwd(), "/GSE67919", sep = "")
#Get the series
gse <- getGEO("GSE67919", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse)) #All females.
annotation <- fData(gse)

sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

# Plots

pheno$`menopausal status:ch1` <- factor(pheno$`menopausal status:ch1`, levels = c("pre", "peri", "post"))
pheno$menop <- pheno$`menopausal status:ch1`
pheno$age_surg <- as.numeric(pheno$`age at surgery:ch1`)
plot4 <- ggbetweenstats(
  data = pheno, 
  x = `menopausal status:ch1`, 
  y = XRa, 
  bf.message = FALSE,
  results.subtitle = FALSE,
  pairwise.display = "none", 
  k = 3L,
  centrality.label.args = list(size = 4.5, nudge_x = 0.4, segment.linetype = 4, min.segment.length = 0)
) +
  labs(title = "X-Ra values by menopausal stage (GSE67919)", x = "Menopausal status", y = "X-Ra") +
  scale_colour_discrete(name = "Menopausal status", labels = c("Pre-menopausal", "Peri-menopausal", "Post-menopausal")) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16), axis.title.y.right = element_blank())

plot5 <- ggscatter(data = na.omit(pheno),
                   x = "age_surg",
                   y = "XRa",
                   color = "menop", 
                   add = "reg.line") +
  stat_cor(aes(color = menop), 
           method = "spearman", 
           cor.coef.name = "rho") +
  stat_cor(method = "spearman", 
           cor.coef.name = "rho", label.x.npc = "middle", label.y.npc = "top") +
  scale_colour_discrete(name = "Menopausal status", labels = c("Pre-menopausal", "Peri-menopausal", "Post-menopausal")) +
  labs(title = "X-Ra distribution with age (GSE67919)", x = "Age", y = "X-Ra") +
  geom_smooth(aes(x = age_surg, y = XRa), col = "black", lty = 2, method = "lm", se = F) +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", linewidth = 1), axis.text = element_text(size = 14), axis.title = element_text(size = 16), title = element_text(size = 17), legend.text = element_text(size = 14), legend.title = element_text(size = 16), axis.title.y.right = element_blank())

plot_grid(plot4, plot5, nrow = 2)
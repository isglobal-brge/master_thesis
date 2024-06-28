# Supplementary Figure 6  - Children X-Ra and BMI

path_here <- paste(getwd(), "/BMI/GSE108213", sep = "")
#Get the series
gse <- getGEO("GSE108213", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse))
annotation <- fData(gse)

sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ]))
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)))

ggbetweenstats(
  data = pheno, 
  x = `race:ch1`, 
  y = XRa, 
  bf.message = FALSE,
  results.subtitle = FALSE,
  k = 3L
) +
  labs(title = "XRa values separated by ethnicity (GSE108213)", x = "Ethnicity") + 
  scale_colour_discrete(name = "Ethnicity") +
  theme_light() +
  theme(axis.line =  element_line(colour = "black", size = 1), axis.text = element_text(size = 10), axis.title = element_text(size = 12), title = element_text(size = 15), legend.text = element_text(size = 10), legend.title = element_text(size = 12), axis.title.y.right = element_blank())
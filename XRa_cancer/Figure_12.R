# Figure 12 - Survival analysis results.

#Load the necessary packages
library(GEOquery)
library(tidyverse)
library(chrXRa)
library(survival)
library(survminer)

# The figure is done using the dataset (GSE78754)

path_here <- paste(getwd(), "/GSE78754", sep = "")
#Get the series
gse <- getGEO("GSE78754", destdir = path_here)[[1]]
#Get the info
met <- exprs(gse)
pheno <- pData(phenoData(gse)) #All females.
annotation <- fData(gse)

# Estimate X-Ra
sel <- annotation$CHR == "X"
metfemale <- t(as.matrix(met[sel, ])) 
pheno <- pheno %>%
  mutate("XRa" = XRa(metfemale, colnames(metfemale)), 
         "status" = ifelse(`survival:ch1` == "alive", 1, 2),
         `age at diagnosis:ch1` = as.numeric(`age at diagnosis:ch1`)) %>% #Alive patients are censored (1), dead patients are informative, sufferend the event (2)
  mutate("HighXRa" = factor(ifelse(XRa > 0.11, "1", "0")))
#Fit the survival model. We group according to different XRa levels (high and low)
sfit <- survfit(Surv(pheno$`month of follow up/to death:ch1`, pheno$status) ~ pheno$HighXRa, data = pheno)
summary(sfit)
summary(sfit, times = seq(0, 180, 20))
ggsurvplot(sfit, 
           pval = TRUE,
           risk.table = TRUE,
           legend.labs = c("XRa < 0.11", "XRa > 0.11"),
           legend.title = "XRa values",
           title = "Kaplan-Meier curve for TNBC survival depending on X-Ra values",
           ncensor.plot = TRUE)
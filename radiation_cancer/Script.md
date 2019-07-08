Script
================
Elia Palma Rojo

  - [Exploratory analysis](#exploratory-analysis)
  - [VOOM-LIMMA](#voom-limma)
      - [GO](#go)
      - [KEGG](#kegg)
  - [DESeq](#deseq)
      - [GO](#go-1)
      - [KEGG](#kegg-1)
      - [WikiPathways](#wikipathways)
      - [MSigDb analysis](#msigdb-analysis)

``` r
library(tidyverse)
library(SummarizedExperiment)
```

The information is loaded into R studio. Only those samples belonging to
normal tissue (at a distance greater than 2cm from tumor margin) and
with information about second tumor development and radiotherapy
administration are selected.

``` r
load("data/rse_gene.Rdata")
# Pick only information about solid tissue normal (constitutional expression)
sel <- grep("-11A",rse_gene$cgc_sample_id,value=TRUE)
# Select those with new tumor events information and with or without radiation therapy
mask <- rse_gene$cgc_sample_id%in%sel & 
  rse_gene$cgc_case_new_tumor_event_after_initial_treatment%in%c("NO","YES") &
  rse_gene$xml_radiation_therapy%in%c("NO", "YES") 
rse_subset<-rse_gene[,mask]
rse_subset$xml_radiation_therapy<-droplevels(rse_subset$xml_radiation_therapy,exclude=c("0","1","2","70","80","90","pT2a","pT2c"))
rse_subset$cgc_case_new_tumor_event_after_initial_treatment<-as.factor(rse_subset$cgc_case_new_tumor_event_after_initial_treatment)
# Create a new file with the specific information
save(rse_subset,file="rse_subset.Rdata")
```

# Exploratory analysis

Data preparation

``` r
load("rse_subset.Rdata")

# Tumor stage variable
stage <- rse_subset$gdc_cases.diagnoses.tumor_stage
ids.early <- grep(paste("stage i$", "stage ia$", "stage ib$" ,
                        "stage ic$", "stage ii$", "stage iia$",
                        "stage iib$" , "stage iic$",
                        sep="|"), stage)
ids.late <- grep(paste("stage iii$", "stage iiia$", "stage iiib$" ,
                        "stage iiic$", "stage iv$", "stage iva$",
                        "stage ivb$" , "stage ivc$",
                        sep="|"), stage)
colData(rse_subset)$GROUP <- rep(NA, ncol(rse_subset))
colData(rse_subset)$GROUP[ids.early] <- "early"
colData(rse_subset)$GROUP[ids.late] <- "late"

#New column combining radiation and tumor information
colData(rse_subset)$radiation_tumor <- rep(NA, ncol(rse_subset))
rse_subset$radiation_tumor[rse_subset$xml_radiation_therapy=="YES" & rse_subset$cgc_case_new_tumor_event_after_initial_treatment=="YES"]<-"YES" 
rse_subset$radiation_tumor[rse_subset$xml_radiation_therapy!="YES" | rse_subset$cgc_case_new_tumor_event_after_initial_treatment!="YES"]<-"NO"
 

rse_subset$radiation_tumor<-as.factor(rse_subset$radiation_tumor)
rse_subset$gdc_cases.demographic.gender<-as.factor(rse_subset$gdc_cases.demographic.gender)
rse_subset$gdc_cases.demographic.ethnicity<-as.factor(rse_subset$gdc_cases.demographic.ethnicity)
rse_subset$cgc_case_primary_site<-as.factor(rse_subset$cgc_case_primary_site)
rse_subset$gdc_cases.demographic.ethnicity<-droplevels(rse_subset$gdc_cases.demographic.ethnicity,exclude="not reported")
rse_subset$xml_tissue_source_site<-droplevels(rse_subset$xml_tissue_source_site)
rse_subset$GROUP<-as.factor(rse_subset$GROUP)

# Generation of the variable phenotype (pheno) and count data (counts)
pheno<-colData(rse_subset)
counts<-assay(rse_subset)

#Subset of variables for the exploratory analysis
sub<-subset(pheno,select=c("gdc_cases.demographic.gender","gdc_cases.demographic.ethnicity","GROUP","cgc_case_age_at_diagnosis","cgc_case_primary_site","radiation_tumor"))

sub<-as.data.frame(sub)
str(sub)
```

    'data.frame':   104 obs. of  6 variables:
     $ gdc_cases.demographic.gender   : Factor w/ 2 levels "female","male": 1 2 2 1 1 2 2 2 2 1 ...
     $ gdc_cases.demographic.ethnicity: Factor w/ 2 levels "hispanic or latino",..: 2 NA NA 2 NA NA 1 2 NA 2 ...
     $ GROUP                          : Factor w/ 2 levels "early","late": 1 2 2 2 1 NA 1 2 2 1 ...
     $ cgc_case_age_at_diagnosis      : int  53 69 60 62 67 50 51 74 75 42 ...
     $ cgc_case_primary_site          : Factor w/ 14 levels "Bile Duct","Bladder",..: 13 4 13 6 4 9 4 11 6 13 ...
     $ radiation_tumor                : Factor w/ 2 levels "NO","YES": 1 1 1 1 1 1 1 1 1 1 ...

For the exploratory analysis, a contingency table is constructed in
order to determine the number of individuals belonging to the subgroup
of patients who had received radiotherapy and developed a second tumor.

``` r
library(compareGroups)
# Table representation of 2 variables of interest
mytable<-as.matrix(table(pheno$cgc_case_new_tumor_event_after_initial_treatment,pheno$xml_radiation_therapy))
names(dimnames(mytable))<-c("new tumor event","radiation therapy")
mytable
```

``` 
               radiation therapy
new tumor event NO YES
            NO  52  35
            YES 13   4
```

``` r
write.csv2(mytable,"2xtable.csv")
```

The gender, the ethnicity, the primary cancer status, the age at
diagnosis, and the first tumor tissue site were graphically represented
separately for the mentioned group and the rest of the patients. For
that purpose, the package compareGroups was used.

``` r
# Results of the exploratory analysis
results<-compareGroups(radiation_tumor~.,data=sub,max.xlev=32,na.action = na.pass)
results.table<-createTable(results)
results
```

``` 


-------- Summary of results by groups of 'radiation_tumor'---------


  var                             N   p.value method            selection
1 gdc_cases.demographic.gender    104 0.618   categorical       ALL      
2 gdc_cases.demographic.ethnicity  82 0.305   categorical       ALL      
3 GROUP                            93 0.003** categorical       ALL      
4 cgc_case_age_at_diagnosis       104 0.689   continuous normal ALL      
5 cgc_case_primary_site           104 0.478   categorical       ALL      
-----
Signif. codes:  0 '**' 0.05 '*' 0.1 ' ' 1 
```

``` r
results.table
```

``` 

--------Summary descriptives table by 'radiation_tumor'---------

__________________________________________________________________ 
                                     NO          YES     p.overall 
                                    N=100        N=4               
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
gdc_cases.demographic.gender:                              0.618   
    female                       49 (49.0%)   3 (75.0%)            
    male                         51 (51.0%)   1 (25.0%)            
gdc_cases.demographic.ethnicity:                           0.305   
    hispanic or latino            6 (7.69%)   1 (25.0%)            
    not hispanic or latino       72 (92.3%)   3 (75.0%)            
GROUP:                                                     0.003   
    early                        70 (78.7%)   0 (0.00%)            
    late                         19 (21.3%)   4 (100%)             
cgc_case_age_at_diagnosis        54.9 (16.1) 52.0 (12.8)   0.689   
cgc_case_primary_site:                                     0.478   
    Bile Duct                     6 (6.00%)   0 (0.00%)            
    Bladder                       3 (3.00%)   0 (0.00%)            
    Head and Neck                 2 (2.00%)   1 (25.0%)            
    Kidney                       28 (28.0%)   0 (0.00%)            
    Liver                         2 (2.00%)   0 (0.00%)            
    Lung                          2 (2.00%)   0 (0.00%)            
    Nervous System                2 (2.00%)   0 (0.00%)            
    Pancreas                      3 (3.00%)   0 (0.00%)            
    Prostate                      3 (3.00%)   0 (0.00%)            
    Skin                          1 (1.00%)   0 (0.00%)            
    Stomach                      11 (11.0%)   0 (0.00%)            
    Thymus                        2 (2.00%)   0 (0.00%)            
    Thyroid                      34 (34.0%)   3 (75.0%)            
    Uterus                        1 (1.00%)   0 (0.00%)            
¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ 
```

``` r
plot(results,bivar=TRUE)
```

![](Script_files/figure-gfm/exploratory%20analysis%203-1.png)<!-- -->![](Script_files/figure-gfm/exploratory%20analysis%203-2.png)<!-- -->![](Script_files/figure-gfm/exploratory%20analysis%203-3.png)<!-- -->![](Script_files/figure-gfm/exploratory%20analysis%203-4.png)<!-- -->![](Script_files/figure-gfm/exploratory%20analysis%203-5.png)<!-- -->

# VOOM-LIMMA

Using the edgeR package, a DGEList object is generated by means of the
read counts matrix. After that, a filtering of the data with zero or low
counts is performed. The following step consisted in normalizing the
data contained in the DGEList with the TMM method. The common dispersion
(squared coefficient of variation) is also estimated.

``` r
library(edgeR)

dge<-DGEList(counts=counts)
dim(dge)
```

    [1] 58037   104

``` r
#Remove rows with 0 or low counts
keep<-filterByExpr(dge)
dge<-dge[keep,,keep.lib.sizes=FALSE]
dim(dge)
```

    [1] 34854   104

``` r
#Scale normalization to RNA-seq read counts with TMM
dge<-calcNormFactors(dge)
dge<-estimateCommonDisp(dge)
dge$common.dispersion
```

    [1] 1.131766

``` r
#coefficient of variation, proportion by which true expression levels may vary between replicate samples.
sqrt(dge$common.dispersion)
```

    [1] 1.063845

A full model matrix with the adjustment variables and the variables of
interest and a null model with only the adjustment variables are
created. After that, the voom transformation is applied to the
normalized DGEList using as a design the full model matrix and a plot
with the mean-variance trend is generated.

``` r
library(limma)

#Create a full model matrix and a null model that contains all and only the adjustment variables respectively
mod<-model.matrix(~xml_radiation_therapy*cgc_case_new_tumor_event_after_initial_treatment,data=pheno)
mod0<- model.matrix(~xml_radiation_therapy+cgc_case_new_tumor_event_after_initial_treatment,data=pheno)
v<-voom(counts=dge,plot=TRUE,design=mod)
```

![](Script_files/figure-gfm/voom-1.png)<!-- -->

In this step, the number of latent factors is estimated and the sva
function is used to estimate the surrogate variables. The full and null
matrices are then modified including the unknown batch found with SVA.

``` r
library(sva)
#Surrogate variable analysis
n.sv<-num.sv(v$E,mod,method="leek")
n.sv
```

    [1] 99

``` r
svobj<-sva(v$E, mod, mod0, method="irw",n.sv=n.sv)
```

``` 
Number of significant surrogate variables is:  99 
Iteration (out of 5 ):1  2  3  4  5  
```

``` r
modSv<-cbind(mod,svobj$sv)
mod0Sv<-cbind(mod0,svobj$sv)
```

The FC and standard errors are estimated by fitting a linear model for
each gene using as a design the new generated matrices. The empirical
Bayes smoothing is then applied to the standard errors.

``` r
#Limma
fit <- lmFit(v,modSv)
fit <- eBayes(fit)
#Results
deGenes<-topTable(fit, coef=4,sort.by="p",number=34854,lfc=log2(1.5))
head(deGenes)
```

``` 
                      logFC   AveExpr        t      P.Value adj.P.Val
ENSG00000280296.1  176.6423 -7.286362 15.86637 0.0001294953 0.4692611
ENSG00000105048.16 132.4833 -4.010106 14.87664 0.0001648899 0.4692611
ENSG00000261175.5  165.2482 -3.812365 13.84425 0.0002158488 0.4692611
ENSG00000181626.11 156.1828 -6.261007 13.74778 0.0002215677 0.4692611
ENSG00000255462.1  108.7544 -3.554885 13.73299 0.0002224610 0.4692611
ENSG00000130957.4  131.5549 -3.796820 13.58673 0.0002315473 0.4692611
                           B
ENSG00000280296.1  -4.436971
ENSG00000105048.16 -4.336914
ENSG00000261175.5  -4.387594
ENSG00000181626.11 -4.438832
ENSG00000255462.1  -4.300342
ENSG00000130957.4  -4.396563
```

Given that all the obtained adjusted P values were higher than 0.05, the
genes shown here are those with a p value (not adjusted) smaller than
0.001

``` r
mask.deGenes<-deGenes$P.Value<0.001
masked.deGenes<-deGenes[mask.deGenes,]
#Obtain ENSEMBL names
deGenes_splt<-rownames(masked.deGenes)
deGenes_splt<-unlist(strsplit(deGenes_splt, "[.]"))
deGenes_splt<-subset(deGenes_splt, grepl("^EN", deGenes_splt))
deGenes_splt
```

``` 
 [1] "ENSG00000280296" "ENSG00000105048" "ENSG00000261175"
 [4] "ENSG00000181626" "ENSG00000255462" "ENSG00000130957"
 [7] "ENSG00000241853" "ENSG00000166869" "ENSG00000114638"
[10] "ENSG00000265933" "ENSG00000236180" "ENSG00000148513"
[13] "ENSG00000230438" "ENSG00000234393" "ENSG00000240163"
[16] "ENSG00000227857" "ENSG00000224647" "ENSG00000272021"
[19] "ENSG00000211753" "ENSG00000280274" "ENSG00000115112"
[22] "ENSG00000177776" "ENSG00000249201" "ENSG00000042832"
[25] "ENSG00000251257" "ENSG00000170807" "ENSG00000185038"
[28] "ENSG00000280085" "ENSG00000227913" "ENSG00000226339"
[31] "ENSG00000257905" "ENSG00000275850" "ENSG00000259645"
[34] "ENSG00000167080" "ENSG00000249609" "ENSG00000236611"
[37] "ENSG00000278716" "ENSG00000258091" "ENSG00000231439"
[40] "ENSG00000279186" "ENSG00000248909" "ENSG00000279591"
[43] "ENSG00000231063" "ENSG00000115541" "ENSG00000259438"
[46] "ENSG00000256553" "ENSG00000113070" "ENSG00000200378"
[49] "ENSG00000187416" "ENSG00000007038" "ENSG00000274374"
[52] "ENSG00000229391" "ENSG00000232629" "ENSG00000279281"
```

``` r
geneUniverse<-rownames(fit)
geneUniverse_splt<-unlist(strsplit(geneUniverse, "[.]"))
geneUniverse_splt<-subset(geneUniverse_splt, grepl("^EN", geneUniverse_splt))
```

List of the gene universe and all the found differentially expressed
genes. These lists are then annotated in Entrez IDs. The gene symbol for
the gene sets are also obtained.

``` r
library(org.Hs.eg.db)
geneUniverse <- unlist(mget(geneUniverse_splt, org.Hs.egENSEMBL2EG, 
                            ifnotfound=NA))
deGenes_2 <- unlist(mget(deGenes_splt, org.Hs.egENSEMBL2EG, 
                      ifnotfound=NA))
```

ENSEMBL, ENTREZ and SYMBOL IDs for gene set and gene universe
respectively

``` r
library(clusterProfiler)
eg<-bitr(deGenes_splt,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
eg
```

``` 
           ENSEMBL  ENTREZID        SYMBOL
2  ENSG00000105048      7138         TNNT1
3  ENSG00000261175 102724344     LINC02188
4  ENSG00000181626    342850       ANKRD62
6  ENSG00000130957      8789          FBP2
8  ENSG00000166869     63928          CHP2
9  ENSG00000114638      7348         UPK1B
10 ENSG00000265933    400643     LINC00668
12 ENSG00000148513     91074      ANKRD30A
13 ENSG00000230438    221756    SERPINB9P1
16 ENSG00000227857 101929626  LOC101929626
21 ENSG00000115112     29842       TFCP2L1
23 ENSG00000249201 101928857 CTD-3080P12.3
24 ENSG00000042832      7038            TG
26 ENSG00000170807    442721         LMOD2
27 ENSG00000185038    339766        MROH2A
34 ENSG00000167080    124872      B4GALNT2
36 ENSG00000236611 102724679     LINC02556
39 ENSG00000231439 100132169        WASIR2
44 ENSG00000115541      3336         HSPE1
47 ENSG00000113070      1839         HBEGF
49 ENSG00000187416    375612        LHFPL3
50 ENSG00000007038     10942        PRSS21
52 ENSG00000229391      3128      HLA-DRB6
53 ENSG00000232629      3120      HLA-DQB2
```

``` r
eg_universe<-bitr(geneUniverse_splt,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(eg_universe)
```

``` 
          ENSEMBL ENTREZID   SYMBOL
1 ENSG00000000003     7105   TSPAN6
2 ENSG00000000005    64102     TNMD
3 ENSG00000000419     8813     DPM1
4 ENSG00000000457    57147    SCYL3
5 ENSG00000000460    55732 C1orf112
6 ENSG00000000938     2268      FGR
```

## GO

Enrichment analysis with Gene Ontology (GO) and html report

``` r
library(GOstats)
params <- new("GOHyperGParams", geneIds=deGenes_2,
              universeGeneIds=geneUniverse, ontology = 'BP',
              annotation="org.Hs.eg.db",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")
hgOver <- hyperGTest(params)
#html file with all the results
htmlReport(hgOver, file="goBPlimma.html")
#Top signalling pathways
head(summary(hgOver))
```

``` 
      GOBPID       Pvalue  OddsRatio     ExpCount Count Size
1 GO:0045214 0.0005325730   71.27727 0.0347019194     2   42
2 GO:0051543 0.0008262362        Inf 0.0008262362     1    1
3 GO:0051545 0.0008262362        Inf 0.0008262362     1    1
4 GO:0030239 0.0011588306   47.45758 0.0512266429     2   62
5 GO:0005986 0.0016518422 1310.00000 0.0016524724     1    2
6 GO:0051542 0.0016518422 1310.00000 0.0016524724     1    2
                                                 Term
1                              sarcomere organization
2          regulation of elastin biosynthetic process
3 negative regulation of elastin biosynthetic process
4                                  myofibril assembly
5                        sucrose biosynthetic process
6                        elastin biosynthetic process
```

Gene classification based on GO distribution

``` r
ggo<-groupGO(gene=eg$ENTREZID,
             OrgDb = org.Hs.eg.db,
             ont="CC",
             readable = TRUE)
head(ggo)
```

``` 
                   ID          Description Count GeneRatio
GO:0016020 GO:0016020             membrane    10     10/24
GO:0005576 GO:0005576 extracellular region     6      6/24
GO:0005623 GO:0005623                 cell    12     12/24
GO:0009295 GO:0009295             nucleoid     0      0/24
GO:0019012 GO:0019012               virion     0      0/24
GO:0030054 GO:0030054        cell junction     1      1/24
                                                                                      geneID
GO:0016020               FBP2/CHP2/UPK1B/TFCP2L1/B4GALNT2/HSPE1/HBEGF/LHFPL3/PRSS21/HLA-DQB2
GO:0005576                                                  FBP2/UPK1B/TG/HSPE1/HBEGF/PRSS21
GO:0005623 TNNT1/FBP2/CHP2/UPK1B/ANKRD30A/TFCP2L1/LMOD2/B4GALNT2/HSPE1/HBEGF/PRSS21/HLA-DQB2
GO:0009295                                                                                  
GO:0019012                                                                                  
GO:0030054                                                                              FBP2
```

``` r
barplot(ggo,drop=TRUE,showCategory=12)
```

![](Script_files/figure-gfm/GO%20classification-1.png)<!-- -->

GO over-representation analysis

``` r
ego<-enrichGO(gene =deGenes_splt,
              universe = geneUniverse_splt,
              keyType = 'ENSEMBL',
              OrgDb = org.Hs.eg.db,
              ont = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.01,
              qvalueCutoff  = 0.05,
        readable = TRUE)
head(ego)
```

``` 
                   ID                                Description GeneRatio
GO:0005865 GO:0005865              striated muscle thin filament      2/14
GO:0036379 GO:0036379                                myofilament      2/14
GO:0030017 GO:0030017                                  sarcomere      3/14
GO:0030669 GO:0030669 clathrin-coated endocytic vesicle membrane      2/14
GO:0044449 GO:0044449                     contractile fiber part      3/14
GO:0030016 GO:0030016                                  myofibril      3/14
             BgRatio       pvalue    p.adjust      qvalue           geneID
GO:0005865  31/16679 0.0003000327 0.005216861 0.002965374      TNNT1/LMOD2
GO:0036379  32/16679 0.0003198814 0.005216861 0.002965374      TNNT1/LMOD2
GO:0030017 188/16679 0.0004681344 0.005216861 0.002965374 TNNT1/FBP2/LMOD2
GO:0030669  40/16679 0.0005011131 0.005216861 0.002965374   HBEGF/HLA-DQB2
GO:0044449 205/16679 0.0006026632 0.005216861 0.002965374 TNNT1/FBP2/LMOD2
GO:0030016 208/16679 0.0006287079 0.005216861 0.002965374 TNNT1/FBP2/LMOD2
           Count
GO:0005865     2
GO:0036379     2
GO:0030017     3
GO:0030669     2
GO:0044449     3
GO:0030016     3
```

``` r
barplot(ego,drop=TRUE,showCategory = 8)
```

![](Script_files/figure-gfm/GO%20over-representation-1.png)<!-- -->

``` r
emapplot(ego)
```

![](Script_files/figure-gfm/GO%20over-representation-2.png)<!-- -->

``` r
goplot(ego)
```

![](Script_files/figure-gfm/GO%20over-representation-3.png)<!-- -->

``` r
#In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot extracts the complex association.
cnetplot(ego,categorySize="pvalue")
```

![](Script_files/figure-gfm/GO%20over-representation-4.png)<!-- -->

## KEGG

``` r
library(KEGG.db)
params.kegg <- new("KEGGHyperGParams", geneIds=deGenes_2,
                   universeGeneIds=geneUniverse,
                   annotation="org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   testDirection="over")
hg_kegg <- hyperGTest(params.kegg)
summary(hg_kegg)
```

``` 
  KEGGID     Pvalue OddsRatio   ExpCount Count Size
1  00030 0.01968332  69.56000 0.01982463     1   26
2  05320 0.02568084  52.61616 0.02592451     1   34
3  00051 0.02717591  49.59048 0.02744949     1   36
4  05014 0.03907463  33.92810 0.03964926     1   52
5  00010 0.04719144  27.84946 0.04803660     1   63
                                 Term
1           Pentose phosphate pathway
2          Autoimmune thyroid disease
3     Fructose and mannose metabolism
4 Amyotrophic lateral sclerosis (ALS)
5        Glycolysis / Gluconeogenesis
```

# DESeq

The package BiocParallel is employed to parallelize the work and speed
up the process. A DESeqDataSet object is constructed with the count
matrix, the phenotypic information and using as a design the interaction
between the variables “radiation therapy” and “new tumor event after
initial treatment”.

``` r
library(BiocParallel)
library(DESeq2)

register(MulticoreParam(10))
identical(colnames(counts),rownames(pheno))
```

    [1] TRUE

``` r
dds<-DESeqDataSetFromMatrix(countData = counts,colData=pheno,design=~xml_radiation_therapy*cgc_case_new_tumor_event_after_initial_treatment)
```

To reduce the memory size of the object and increase the speed of the
transformation and testing functions, a pre-filtering is performed
keeping only those rows with at least 10 reads.

``` r
nrow(dds)
```

    [1] 58037

``` r
dds<-dds[rowSums(counts(dds))>=10,]
nrow(dds)
```

    [1] 57185

Differential expression analysis including the SVA ifnromation

``` r
dds<-DESeq(dds,parallel=TRUE,BPPARAM=MulticoreParam(10))

dds2 <- dds
temp<-colData(dds2)
temp2 <- cbind(temp, svobj$sv)
colData(dds2) <- temp2

design(dds2)
```

    ~xml_radiation_therapy * cgc_case_new_tumor_event_after_initial_treatment

``` r
mm1 <- paste0(design(dds2), collapse="")
mm2 <- paste0(paste0("V", 1:99), collapse="+")
design(dds2) <- formula(paste(mm1, mm2, sep="+"))

dds2 <- DESeq(dds2,parallel=TRUE,BPPARAM=MulticoreParam(10))
```

Differential expression analysis
results

``` r
res<-results(dds2,name="xml_radiation_therapyYES.cgc_case_new_tumor_event_after_initial_treatmentYES",lfcThreshold =log2(1.5),alpha=0.05)
# Show the log2 fold changes. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
plotMA(res)
```

![](Script_files/figure-gfm/results%20deseq-1.png)<!-- -->

``` r
#Shrinked
resLFC <- lfcShrink(dds, coef="xml_radiation_therapyYES.cgc_case_new_tumor_event_after_initial_treatmentYES", type="apeglm")
resLFC
```

    log2 fold change (MAP): xml radiation therapyYES.cgc case new tumor event after initial treatmentYES 
    Wald test p-value: xml radiation therapyYES.cgc case new tumor event after initial treatmentYES 
    DataFrame with 57185 rows and 5 columns
                               baseMean        log2FoldChange
                              <numeric>             <numeric>
    ENSG00000000003.14 487605.898337197 -2.04429789677525e-06
    ENSG00000000005.5  3164.29641211176  7.24983475965346e-07
    ENSG00000000419.12 181885.892135295 -3.51980762294671e-06
    ENSG00000000457.13 124605.499340207  4.29833759209134e-07
    ENSG00000000460.16 54070.5686729281 -5.69774122001906e-06
    ...                             ...                   ...
    ENSG00000283695.1  15.9109226716692 -2.66786245040828e-07
    ENSG00000283696.1  3835.53885164004 -9.11558377482549e-08
    ENSG00000283697.1  7374.51174638862 -7.67276158900929e-07
    ENSG00000283698.1  72.9255848358186 -6.68122006652495e-07
    ENSG00000283699.1   11.042106713358  7.29546370981828e-08
                                     lfcSE               pvalue
                                 <numeric>            <numeric>
    ENSG00000000003.14 0.00144269196166245    0.221004593442356
    ENSG00000000005.5  0.00144269416838909    0.388178433933817
    ENSG00000000419.12 0.00144266807887319    0.439062725026935
    ENSG00000000457.13 0.00144265365109118    0.937050238682581
    ENSG00000000460.16 0.00144267942389674    0.147049847997961
    ...                                ...                  ...
    ENSG00000283695.1   0.0014426950158536    0.403321334806188
    ENSG00000283696.1  0.00144268833976172    0.967521500734405
    ENSG00000283697.1  0.00144268547263694     0.77327680775513
    ENSG00000283698.1  0.00144269511371146 5.71314798041982e-05
    ENSG00000283699.1  0.00144269492508373    0.799600931329605
                                      padj
                                 <numeric>
    ENSG00000000003.14   0.746758077158952
    ENSG00000000005.5    0.914980325078434
    ENSG00000000419.12   0.945517044059442
    ENSG00000000457.13   0.999737909061608
    ENSG00000000460.16   0.621789693009265
    ...                                ...
    ENSG00000283695.1    0.925803058867808
    ENSG00000283696.1    0.999737909061608
    ENSG00000283697.1    0.999737909061608
    ENSG00000283698.1  0.00137360182308898
    ENSG00000283699.1    0.999737909061608

``` r
plotMA(resLFC,ylim=c(-3,3))
abline(h=c(-1,1),col="dodgerblue",lwd=2)
```

![](Script_files/figure-gfm/results%20deseq-2.png)<!-- -->

``` r
# Order the results according to padj (from lower to higher)
res2<-subset(res,padj<0.05)
res2<-subset(res2,na.rm=TRUE)
res2<-res2[order(res2$padj),]
dim(res2)
```

    [1] 304   6

``` r
head(res2)
```

    log2 fold change (MLE): xml radiation therapyYES.cgc case new tumor event after initial treatmentYES 
    Wald test p-value: xml radiation therapyYES.cgc case new tumor event after initial treatmentYES 
    DataFrame with 6 rows and 6 columns
                               baseMean    log2FoldChange            lfcSE
                              <numeric>         <numeric>        <numeric>
    ENSG00000113070.7  326512.969308293  22.3031905031618 2.54545601554449
    ENSG00000115112.7  928692.848156493 -22.4595415631781 2.68461202167874
    ENSG00000229391.7  402184.902592179 -23.7505515714418 2.91429584014344
    ENSG00000242534.2  182607.336415292 -24.0739143896344 3.01565779561627
    ENSG00000116990.10 86802.4444174346 -26.7695241134057 3.44970916012581
    ENSG00000211685.3  95603.0679479471 -30.2497392212123 4.07019677929685
                                    stat               pvalue
                               <numeric>            <numeric>
    ENSG00000113070.7   8.53215607333722 1.43645083606164e-17
    ENSG00000115112.7  -8.14813421299452 3.69582044434701e-16
    ENSG00000229391.7  -7.94894902282141 1.88100429308918e-15
    ENSG00000242534.2   -7.7889977845159 6.75427910426246e-15
    ENSG00000116990.10 -7.59036788241291 3.18998497674841e-14
    ENSG00000211685.3    -7.288290549337 3.13912430612722e-13
                                       padj
                                  <numeric>
    ENSG00000113070.7  4.07363092598721e-13
    ENSG00000115112.7  5.24048859906185e-12
    ENSG00000229391.7   1.7781133582572e-11
    ENSG00000242534.2  4.78861502794448e-11
    ENSG00000116990.10 1.80929567911216e-10
    ENSG00000211685.3  1.48370710329103e-09

``` r
write.csv2(res2,"all_ENSEMBL_genes.csv")
```

Differential expressed genes

``` r
deGenes.deseq<-rownames(res2)
deGenes.deseq<-unlist(strsplit(deGenes.deseq, "[.]"))
deGenes.deseq<-subset(deGenes.deseq, grepl("^EN", deGenes.deseq))
geneUniverse.deseq<-rownames(dds2)
geneUniverse.deseq2<-unlist(strsplit(geneUniverse.deseq, "[.]"))
geneUniverse.deseq2<-subset(geneUniverse.deseq2, grepl("^EN", geneUniverse.deseq2))
```

ENSEMBL, ENTREZ and UNIPROT IDs for gene set and gene
universe

``` r
eg2<-bitr(deGenes.deseq,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(eg2)
```

``` 
           ENSEMBL ENTREZID   SYMBOL
1  ENSG00000113070     1839    HBEGF
2  ENSG00000115112    29842  TFCP2L1
3  ENSG00000229391     3128 HLA-DRB6
5  ENSG00000116990     4610     MYCL
7  ENSG00000128422     3872    KRT17
10 ENSG00000132437     1644      DDC
```

``` r
dim(eg2)
```

    [1] 248   3

``` r
#Report gene names
write.csv2(eg2, "all_gene_names.csv")

eg2_universe<-bitr(geneUniverse.deseq2,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
```

``` r
#Plot the results with some of the genes up and down regulated
library(ggpubr)
ggmaplot(res, main = "radiation therapy and new tumor event",
   fdr = 0.05, fc = 2, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames = as.vector(rownames(res)),
   legend = "top", top = 20,
   font.label = c("bold", 11),
   font.legend = "bold",
   font.main = "bold",
   ggtheme = ggplot2::theme_minimal())
```

![](Script_files/figure-gfm/GO%20gene%20names%20up%20and%20down-1.png)<!-- -->

``` r
#Upregulated genes
ur.genes<-rownames(subset(res2,log2FoldChange>1))
ur.genes<-unlist(strsplit(ur.genes, "[.]"))
ur.genes<-subset(ur.genes, grepl("^EN", ur.genes))
ur.genes2<-bitr(ur.genes,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(ur.genes2)
```

``` 
          ENSEMBL ENTREZID   SYMBOL
1 ENSG00000113070     1839    HBEGF
2 ENSG00000080824     3320 HSP90AA1
3 ENSG00000173467   155465     AGR3
4 ENSG00000133063     1118    CHIT1
6 ENSG00000092421    57556   SEMA6A
7 ENSG00000180481   144321 GLIPR1L2
```

``` r
write.csv2(ur.genes2, "upregulated_genes.csv")

#Downregulated genes
dr.genes<-rownames(subset(res2,log2FoldChange<(-1)))
dr.genes<-unlist(strsplit(dr.genes, "[.]"))
dr.genes<-subset(dr.genes, grepl("^EN", dr.genes))
dr.genes2<-bitr(dr.genes,fromType="ENSEMBL",toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(dr.genes2)
```

``` 
           ENSEMBL ENTREZID   SYMBOL
1  ENSG00000115112    29842  TFCP2L1
2  ENSG00000229391     3128 HLA-DRB6
4  ENSG00000116990     4610     MYCL
6  ENSG00000128422     3872    KRT17
9  ENSG00000132437     1644      DDC
10 ENSG00000118194     7139    TNNT2
```

``` r
write.csv2(dr.genes2, "downregulated_genes.csv")
```

Data preparation for the enrichment
analysis

``` r
geneUniverse.deseq <- unlist(mget(geneUniverse.deseq2, org.Hs.egENSEMBL2EG, 
                            ifnotfound=NA))
deGenes.deseq2 <- unlist(mget(deGenes.deseq, org.Hs.egENSEMBL2EG, 
                      ifnotfound=NA))
deGenes.ur<- unlist(mget(ur.genes, org.Hs.egENSEMBL2EG, 
                      ifnotfound=NA))
deGenes.dr<- unlist(mget(dr.genes, org.Hs.egENSEMBL2EG, 
                      ifnotfound=NA))
```

## GO

Gene ontology analysis. Ontology: - BP:biological process - MF:
molecular function - CC: cellular component

All genes

``` r
params2 <- new("GOHyperGParams", geneIds=deGenes.deseq2,
              universeGeneIds=geneUniverse.deseq, ontology = 'BP',
              annotation="org.Hs.eg.db",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")
hgOver2 <- hyperGTest(params2)
htmlReport(hgOver2, file="go_deseq.html")
head(summary(hgOver2))
```

``` 
      GOBPID       Pvalue OddsRatio  ExpCount Count Size
1 GO:0007166 1.075586e-09  2.582019  34.30210    70 2845
2 GO:0050896 2.642684e-09  2.319223 105.40209   148 8742
3 GO:0006952 5.128623e-09  2.959050  18.44717    46 1530
4 GO:0006950 1.506320e-08  2.295599  45.24983    81 3753
5 GO:0006955 1.641017e-08  2.676796  23.24585    52 1928
6 GO:0048870 1.654308e-08  2.835497  19.15854    46 1589
                                     Term
1 cell surface receptor signaling pathway
2                    response to stimulus
3                        defense response
4                      response to stress
5                         immune response
6                           cell motility
```

Upregulated genes

``` r
params2.ur <- new("GOHyperGParams", geneIds=deGenes.ur,
              universeGeneIds=geneUniverse.deseq, ontology = 'BP',
              annotation="org.Hs.eg.db",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")
hgOver2.ur <- hyperGTest(params2.ur)
htmlReport(hgOver2.ur, file="go_deseq_ur.html")
head(summary(hgOver2.ur))
```

``` 
      GOBPID       Pvalue OddsRatio  ExpCount Count Size
1 GO:0048870 4.578248e-07  4.028746 7.3069762    23 1589
2 GO:0051674 4.578248e-07  4.028746 7.3069762    23 1589
3 GO:0040011 4.356168e-06  3.482297 8.3232391    23 1810
4 GO:0006928 9.255738e-06  3.236718 9.3624944    24 2036
5 GO:0006986 1.220666e-05 10.196522 0.7725437     7  168
6 GO:0048675 1.557761e-05 12.660287 0.5334231     6  116
                                       Term
1                             cell motility
2                      localization of cell
3                                locomotion
4 movement of cell or subcellular component
5              response to unfolded protein
6                            axon extension
```

Downregulated genes

``` r
params2.dr <- new("GOHyperGParams", geneIds=deGenes.dr,
              universeGeneIds=geneUniverse.deseq, ontology = 'BP',
              annotation="org.Hs.eg.db",
              pvalueCutoff=0.05, conditional=FALSE,
              testDirection="over")
hgOver2.dr <- hyperGTest(params2.dr)
htmlReport(hgOver2.dr, file="go_deseq_dr.html")
head(summary(hgOver2.dr))
```

``` 
      GOBPID       Pvalue OddsRatio   ExpCount Count Size
1 GO:0006952 7.113238e-11  4.183565 11.4115074    37 1530
2 GO:0006955 9.723283e-10  3.601923 14.3799910    40 1928
3 GO:0045087 1.908765e-09  4.883363  6.1607223    25  826
4 GO:0002376 7.756602e-09  3.051332 20.9733064    48 2812
5 GO:0006950 3.067639e-08  2.754469 27.9917564    56 3753
6 GO:0010273 6.053122e-08 69.097656  0.1118775     5   15
                          Term
1             defense response
2              immune response
3       innate immune response
4        immune system process
5           response to stress
6 detoxification of copper ion
```

General classification

``` r
ggo2<-groupGO(gene=eg2$ENTREZID,
             OrgDb = org.Hs.eg.db,
             ont="CC",
             level=4,
             readable = TRUE)
head(ggo2)
```

``` 
                   ID                                         Description
GO:0001533 GO:0001533                                  cornified envelope
GO:0042383 GO:0042383                                          sarcolemma
GO:0044459 GO:0044459                                plasma membrane part
GO:0097524 GO:0097524                               sperm plasma membrane
GO:0070056 GO:0070056                      prospore membrane leading edge
GO:0070057 GO:0070057 prospore membrane spindle pole body attachment site
           Count GeneRatio
GO:0001533     1     1/246
GO:0042383     3     3/246
GO:0044459    59    59/246
GO:0097524     1     1/246
GO:0070056     0     0/246
GO:0070057     0     0/246
                                                                                                                                                                                                                                                                                                                                                                                              geneID
GO:0001533                                                                                                                                                                                                                                                                                                                                                                                      KRT1
GO:0042383                                                                                                                                                                                                                                                                                                                                                                          SLC8A1/SGCG/ANK1
GO:0044459 HBEGF/SEMA6A/TM4SF1/NGFR/GZMB/ERRFI1/LGR6/B3GNT3/HLA-DPB1/SEMA6D/CDHR5/CLEC4M/CCR7/MTTP/BDKRB1/CACNA1I/CXCL9/IGSF21/LRRC4/HLA-DPA1/MERTK/DLGAP1/SLC8A1/CFTR/CCR2/OLR1/NSG1/CCR5/SLC44A4/SLCO4A1/ATP6V0D2/TSPAN33/DLC1/SEMA3B/SCN7A/CHRNA2/HSP90AB1/GABRB3/HLA-H/XCR1/HLA-DRA/TGM3/CEACAM7/CEACAM6/PSD3/S1PR1/FCER1A/KCNMB1/HYAL2/CMKLR1/BST2/FXYD3/HLA-G/SGCG/TYROBP/ANK1/CALHM6/MLC1/ENG
GO:0097524                                                                                                                                                                                                                                                                                                                                                                                  HSP90AB1
GO:0070056                                                                                                                                                                                                                                                                                                                                                                                          
GO:0070057                                                                                                                                                                                                                                                                                                                                                                                          
```

``` r
barplot(ggo2,drop=TRUE,showCategory=20)
```

![](Script_files/figure-gfm/GO%20deseq%20classification-1.png)<!-- -->

GO over-representation analysis

``` r
ego2<-enrichGO(gene =deGenes.deseq,
              universe = geneUniverse.deseq2,
              keyType = 'ENSEMBL',
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff = 0.01,
              qvalueCutoff  = 0.05,
        readable = TRUE)
head(ego2)
```

``` 
                   ID                          Description GeneRatio
GO:0045926 GO:0045926        negative regulation of growth    16/216
GO:0034341 GO:0034341         response to interferon-gamma    14/216
GO:0060326 GO:0060326                      cell chemotaxis    16/216
GO:0010273 GO:0010273         detoxification of copper ion     5/216
GO:1990169 GO:1990169        stress response to copper ion     5/216
GO:0061687 GO:0061687 detoxification of inorganic compound     5/216
             BgRatio       pvalue     p.adjust       qvalue
GO:0045926 258/17863 1.044190e-07 0.0001756638 0.0001549827
GO:0034341 195/17863 1.117099e-07 0.0001756638 0.0001549827
GO:0060326 291/17863 5.328927e-07 0.0004224462 0.0003727113
GO:0010273  15/17863 6.716156e-07 0.0004224462 0.0003727113
GO:1990169  15/17863 6.716156e-07 0.0004224462 0.0003727113
GO:0061687  16/17863 9.672997e-07 0.0005070263 0.0004473337
                                                                                                 geneID
GO:0045926  SEMA6A/MAP2/SEMA6D/BDKRB1/MT1A/SCGB3A1/MT2A/MT1X/CDA/MT1E/SEMA3B/AGT/MT1H/HYAL2/BST2/SEMA5A
GO:0034341  HLA-DPB1/MT2A/HLA-DPA1/HSP90AB1/HLA-H/KIF16B/HLA-DRA/KYNU/CCL2/SLC30A8/BST2/HLA-G/EVL/SOCS1
GO:0060326 HBEGF/SMOC2/HSPB1/CAMK1D/CCR7/CHGA/CXCL9/CXCL2/CCR2/CCR5/XCR1/S1PR1/CCL2/CMKLR1/NR4A1/SEMA5A
GO:0010273                                                                     MT1A/MT2A/MT1X/MT1E/MT1H
GO:1990169                                                                     MT1A/MT2A/MT1X/MT1E/MT1H
GO:0061687                                                                     MT1A/MT2A/MT1X/MT1E/MT1H
           Count
GO:0045926    16
GO:0034341    14
GO:0060326    16
GO:0010273     5
GO:1990169     5
GO:0061687     5
```

``` r
barplot(ego2,drop=TRUE,showCategory = 10)
```

![](Script_files/figure-gfm/deseq%20GO%20over-representation-1.png)<!-- -->

``` r
emapplot(ego2)
```

![](Script_files/figure-gfm/deseq%20GO%20over-representation-2.png)<!-- -->

``` r
goplot(ego2)
```

![](Script_files/figure-gfm/deseq%20GO%20over-representation-3.png)<!-- -->

``` r
#In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function is used to extract the complex associations.
cnetplot(ego2,categorySize="pvalue",colorEdge=TRUE)
```

![](Script_files/figure-gfm/deseq%20GO%20over-representation-4.png)<!-- -->

## KEGG

All genes

``` r
params.kegg2 <- new("KEGGHyperGParams", geneIds=deGenes.deseq2,
                   universeGeneIds=geneUniverse.deseq,
                   annotation="org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   testDirection="over")
hg_kegg2 <- hyperGTest(params.kegg2)
summary(hg_kegg2)
```

``` 
   KEGGID       Pvalue OddsRatio  ExpCount Count Size
1   04612 8.876578e-05  7.618280 1.0786806     7   69
2   05330 1.874136e-04 11.046512 0.5471568     5   35
3   05332 2.455269e-04 10.352471 0.5784230     5   37
4   04940 4.018411e-04  9.195736 0.6409552     5   41
5   05310 8.588496e-04 10.931034 0.4377255     4   28
6   05320 1.017217e-03  7.344961 0.7816526     5   50
7   04621 1.994429e-03  6.227512 0.9067171     5   58
8   04060 2.620155e-03  2.964370 4.1427590    11  265
9   05416 4.015144e-03  5.229790 1.0630476     5   68
10  04145 8.711825e-03  3.255828 2.3449579     7  150
11  04062 9.058614e-03  2.954936 2.9546470     8  189
12  05150 9.129931e-03  5.330518 0.8285518     4   53
13  05323 1.238254e-02  3.907807 1.3913417     5   89
14  05145 1.567064e-02  3.191271 2.0322969     6  130
15  04514 1.622138e-02  3.165176 2.0479299     6  131
16  04610 2.233205e-02  4.007073 1.0786806     4   69
17  00980 2.340978e-02  3.945664 1.0943137     4   70
18  00380 2.735509e-02  4.974650 0.6565882     3   42
19  04260 3.183740e-02  3.562903 1.2037451     4   77
20  04672 3.460149e-02  4.508721 0.7191204     3   46
                                           Term
1           Antigen processing and presentation
2                           Allograft rejection
3                     Graft-versus-host disease
4                      Type I diabetes mellitus
5                                        Asthma
6                    Autoimmune thyroid disease
7           NOD-like receptor signaling pathway
8        Cytokine-cytokine receptor interaction
9                             Viral myocarditis
10                                    Phagosome
11                  Chemokine signaling pathway
12              Staphylococcus aureus infection
13                         Rheumatoid arthritis
14                                Toxoplasmosis
15               Cell adhesion molecules (CAMs)
16          Complement and coagulation cascades
17 Metabolism of xenobiotics by cytochrome P450
18                        Tryptophan metabolism
19                   Cardiac muscle contraction
20 Intestinal immune network for IgA production
```

Upregulated genes

``` r
params.kegg2.ur <- new("KEGGHyperGParams", geneIds=deGenes.ur,
                   universeGeneIds=geneUniverse.deseq,
                   annotation="org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   testDirection="over")
hg_kegg2.ur <- hyperGTest(params.kegg2.ur)
summary(hg_kegg2.ur)
```

``` 
  KEGGID       Pvalue OddsRatio  ExpCount Count Size
1  04621 0.0005629719 12.122751 0.3885930     4   58
2  05110 0.0054480876  9.364379 0.3617935     3   54
3  04010 0.0081582674  3.830673 1.7955678     6  268
4  04612 0.0107407058  7.217172 0.4622917     3   69
5  04141 0.0236761661  3.990062 1.1054802     4  165
6  00524 0.0330647225 38.013158 0.0334994     1    5
7  00520 0.0406597843  6.740306 0.3215942     2   48
8  04110 0.0494568240  3.898760 0.8307851     3  124
                                         Term
1         NOD-like receptor signaling pathway
2                   Vibrio cholerae infection
3                      MAPK signaling pathway
4         Antigen processing and presentation
5 Protein processing in endoplasmic reticulum
6         Butirosin and neomycin biosynthesis
7 Amino sugar and nucleotide sugar metabolism
8                                  Cell cycle
```

``` r
htmlReport(hg_kegg2.ur, file="kegg_deseq_ur.html")
```

Downregulated genes

``` r
params.kegg2.dr <- new("KEGGHyperGParams", geneIds=deGenes.dr,
                   universeGeneIds=geneUniverse.deseq,
                   annotation="org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   testDirection="over")
hg_kegg2.dr <- hyperGTest(params.kegg2.dr)
summary(hg_kegg2.dr)
```

``` 
   KEGGID       Pvalue OddsRatio  ExpCount Count Size
1   05330 1.239246e-05 20.351064 0.3126611     5   35
2   05332 1.641687e-05 19.072473 0.3305274     5   37
3   04940 2.747272e-05 16.941489 0.3662601     5   41
4   05320 7.309994e-05 13.531915 0.4466586     5   50
5   05310 9.896837e-05 19.947917 0.2501288     4   28
6   05150 1.200451e-03  9.727891 0.4734582     4   53
7   04060 2.183284e-03  3.899540 2.3672908     8  265
8   05416 3.025398e-03  7.428385 0.6074558     4   68
9   04612 3.190581e-03  7.312821 0.6163889     4   69
10  04260 4.739379e-03  6.502283 0.6878543     4   77
11  05145 5.768769e-03  4.803404 1.1613125     5  130
12  00380 6.039452e-03  8.995290 0.3751933     3   42
13  04062 6.361214e-03  3.981468 1.6883697     6  189
14  04672 7.788527e-03  8.152824 0.4109260     3   46
15  04145 1.042956e-02  4.126192 1.3399759     5  150
16  04610 2.329336e-02  5.290353 0.6163889     3   69
17  05140 2.418582e-02  5.210478 0.6253221     3   70
18  04514 2.886336e-02  3.702100 1.1702457     4  131
19  05322 3.102382e-02  3.614744 1.1970452     4  134
20  05323 4.461513e-02  4.045800 0.7950524     3   89
21  05414 4.587017e-02  3.998593 0.8039856     3   90
                                           Term
1                           Allograft rejection
2                     Graft-versus-host disease
3                      Type I diabetes mellitus
4                    Autoimmune thyroid disease
5                                        Asthma
6               Staphylococcus aureus infection
7        Cytokine-cytokine receptor interaction
8                             Viral myocarditis
9           Antigen processing and presentation
10                   Cardiac muscle contraction
11                                Toxoplasmosis
12                        Tryptophan metabolism
13                  Chemokine signaling pathway
14 Intestinal immune network for IgA production
15                                    Phagosome
16          Complement and coagulation cascades
17                                Leishmaniasis
18               Cell adhesion molecules (CAMs)
19                 Systemic lupus erythematosus
20                         Rheumatoid arthritis
21                       Dilated cardiomyopathy
```

``` r
htmlReport(hg_kegg2.dr, file="kegg_deseq_dr.html")
```

Over representation analysis for all, upregulated and downregulated
genes respectively

``` r
kk<-enrichKEGG(gene=eg2$ENTREZID,
               organism="hsa",
               universe=eg2_universe$ENTREZID,
               pvalueCutoff=0.05)
head(kk)
```

``` 
               ID                         Description GeneRatio BgRatio
hsa04612 hsa04612 Antigen processing and presentation     7/121 69/7751
hsa04978 hsa04978                  Mineral absorption     6/121 53/7751
hsa05330 hsa05330                 Allograft rejection     5/121 35/7751
hsa05332 hsa05332           Graft-versus-host disease     5/121 37/7751
hsa04940 hsa04940            Type I diabetes mellitus     5/121 41/7751
hsa05310 hsa05310                              Asthma     4/121 28/7751
               pvalue   p.adjust     qvalue
hsa04612 9.192207e-05 0.01390496 0.01286463
hsa04978 1.611472e-04 0.01390496 0.01286463
hsa05330 1.904367e-04 0.01390496 0.01286463
hsa05332 2.494163e-04 0.01390496 0.01286463
hsa04940 4.079770e-04 0.01819577 0.01683442
hsa05310 8.658295e-04 0.02875238 0.02660121
                                     geneID Count
hsa04612 3320/3115/3113/3326/3122/3305/3135     7
hsa04978      4489/4502/4501/6546/4493/4496     6
hsa05330           3002/3115/3113/3122/3135     5
hsa05332           3002/3115/3113/3122/3135     5
hsa04940           3002/3115/3113/3122/3135     5
hsa05310                3115/3113/3122/2205     4
```

``` r
diffexp.genes2<-kk$geneID
barplot(kk,showCategory = 10)
```

![](Script_files/figure-gfm/deseq%20KEGG%20over-representation-1.png)<!-- -->

``` r
#Upregulated
kk.ur<-enrichKEGG(gene=ur.genes2$ENTREZID,
               organism="hsa",
               universe=eg2_universe$ENTREZID,
               pvalueCutoff=0.05)
head(kk.ur)
```

``` 
               ID             Description GeneRatio BgRatio       pvalue
hsa04657 hsa04657 IL-17 signaling pathway      5/46 93/7751 0.0002074347
           p.adjust     qvalue                     geneID Count
hsa04657 0.03070034 0.02925921 3320/2920/727897/3326/6347     5
```

``` r
diffexp.genes2.ur<-kk.ur$geneID
barplot(kk.ur,showCategory = 10)
```

![](Script_files/figure-gfm/deseq%20KEGG%20over-representation-2.png)<!-- -->

``` r
#Downregulated
kk.dr<-enrichKEGG(gene=dr.genes2$ENTREZID,
               organism="hsa",
               universe=eg2_universe$ENTREZID,
               pvalueCutoff=0.05)
head(kk.dr)
```

``` 
               ID                Description GeneRatio BgRatio
hsa04978 hsa04978         Mineral absorption      6/75 53/7751
hsa05330 hsa05330        Allograft rejection      5/75 35/7751
hsa05332 hsa05332  Graft-versus-host disease      5/75 37/7751
hsa04940 hsa04940   Type I diabetes mellitus      5/75 41/7751
hsa05320 hsa05320 Autoimmune thyroid disease      5/75 50/7751
hsa05310 hsa05310                     Asthma      4/75 28/7751
               pvalue    p.adjust      qvalue
hsa04978 1.073107e-05 0.001252072 0.001104246
hsa05330 1.918843e-05 0.001252072 0.001104246
hsa05332 2.537984e-05 0.001252072 0.001104246
hsa04940 4.233853e-05 0.001566526 0.001381573
hsa05320 1.118651e-04 0.003311208 0.002920269
hsa05310 1.388535e-04 0.003425054 0.003020674
                                geneID Count
hsa04978 4489/4502/4501/6546/4493/4496     6
hsa05330      3002/3115/3113/3122/3135     5
hsa05332      3002/3115/3113/3122/3135     5
hsa04940      3002/3115/3113/3122/3135     5
hsa05320      3002/3115/3113/3122/3135     5
hsa05310           3115/3113/3122/2205     4
```

``` r
diffexp.genes2.dr<-kk.dr$geneID
barplot(kk.dr,showCategory = 10)
```

![](Script_files/figure-gfm/deseq%20KEGG%20over-representation-3.png)<!-- -->

## WikiPathways

Is a more focused approach than GO. WikiPathways is a continuously
updated pathway database.

``` r
library(magrittr)
gene<-eg2$ENTREZID
gene.up<-dr.genes2$ENTREZID
gene.down<-ur.genes2$ENTREZID

#Open WikiPathways (updated pathway database) file
wp2gene <- read.gmt("wikipathways-20190410-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")

#TERM2GENE
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)

#TERM2NAME
wpid2name <- wp2gene %>% dplyr::select(wpid, name)

#Enricher function
ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
ewp <- setReadable(ewp, org.Hs.eg.db, keytype = "ENTREZID")
#To add the gene symbols we are going to use the package DOSE
library(DOSE)
ewp<-DOSE::setReadable(ewp, org.Hs.eg.db, keytype = "ENTREZID")
#Representation
head(ewp)
```

``` 
           ID         Description GeneRatio BgRatio       pvalue
WP3529 WP3529    Zinc homeostasis     6/113 37/6613 3.289482e-05
WP3286 WP3286  Copper homeostasis     6/113 54/6613 2.886334e-04
WP3624 WP3624       Lung fibrosis     6/113 63/6613 6.703881e-04
WP2328 WP2328 Allograft Rejection     7/113 90/6613 8.201469e-04
          p.adjust     qvalue
WP3529 0.006118436 0.00571331
WP3286 0.026842907 0.02506553
WP3624 0.038136829 0.03561164
WP2328 0.038136829 0.03561164
                                                geneID Count
WP3529                MT1A/MT2A/MT1X/MT1E/MT1H/SLC30A8     6
WP3286                  MT1JP/MT1A/MT2A/MT1X/MT1E/MT1H     6
WP3624             MT2A/CXCL2/MUC5B/CCR2/SERPINA1/CCL2     6
WP2328 GZMB/HLA-DPB1/CXCL9/HLA-DPA1/HLA-DRA/C1QC/HLA-G     7
```

``` r
barplot(ewp)
```

![](Script_files/figure-gfm/enricher%20clusterprofiler-1.png)<!-- -->

``` r
dotplot(ewp, showCategory = 20)
```

![](Script_files/figure-gfm/enricher%20clusterprofiler-2.png)<!-- -->

``` r
emapplot(ewp, showCategory = 20)
```

![](Script_files/figure-gfm/enricher%20clusterprofiler-3.png)<!-- -->

## MSigDb analysis

MSigDB is a collection of annotated gene sets. The enrichment is
performed with two different collections, the oncogenic signatures and
the immunologic signatures, and with all the genes, the upregulated ones
and the downregulated separately.

``` r
#Molecular Signatures Database. 
library(msigdbr)
#To retrieve all human gene sets
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
```

``` 
         gs_name  gs_id gs_cat gs_subcat human_gene_symbol species_name
1 AAACCAC_MIR140 M12609     C3       MIR             ABCC4 Homo sapiens
2 AAACCAC_MIR140 M12609     C3       MIR             ACTN4 Homo sapiens
  entrez_gene gene_symbol sources
1       10257       ABCC4    <NA>
2          81       ACTN4    <NA>
```

Oncogenic collection: Gene sets represent signatures of cellular
pathways which are often dis-regulated in cancer. The majority of
signatures were generated directly from microarray data from NCBI GEO or
from internal unpublished profiling experiments which involved
perturbation of known cancer genes. In addition, a small number of
oncogenic signatures were curated from scientific publications.

``` r
#Oncogenic collection
m_t2g <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
```

    # A tibble: 6 x 2
      gs_name      entrez_gene
      <chr>              <int>
    1 AKT_UP.V1_DN      137872
    2 AKT_UP.V1_DN         249
    3 AKT_UP.V1_DN         271
    4 AKT_UP.V1_DN       51129
    5 AKT_UP.V1_DN         362
    6 AKT_UP.V1_DN       25842

``` r
em<-enricher(gene,TERM2GENE=m_t2g)
head(em)
```

``` 
                                                         ID
KRAS.600.LUNG.BREAST_UP.V1_DN KRAS.600.LUNG.BREAST_UP.V1_DN
KRAS.600_UP.V1_DN                         KRAS.600_UP.V1_DN
RPS14_DN.V1_UP                               RPS14_DN.V1_UP
KRAS.300_UP.V1_DN                         KRAS.300_UP.V1_DN
KRAS.LUNG.BREAST_UP.V1_DN         KRAS.LUNG.BREAST_UP.V1_DN
                                                Description GeneRatio
KRAS.600.LUNG.BREAST_UP.V1_DN KRAS.600.LUNG.BREAST_UP.V1_DN    14/167
KRAS.600_UP.V1_DN                         KRAS.600_UP.V1_DN    14/167
RPS14_DN.V1_UP                               RPS14_DN.V1_UP    10/167
KRAS.300_UP.V1_DN                         KRAS.300_UP.V1_DN     8/167
KRAS.LUNG.BREAST_UP.V1_DN         KRAS.LUNG.BREAST_UP.V1_DN     8/167
                                BgRatio       pvalue    p.adjust
KRAS.600.LUNG.BREAST_UP.V1_DN 289/11250 0.0001010267 0.008637786
KRAS.600_UP.V1_DN             289/11250 0.0001010267 0.008637786
RPS14_DN.V1_UP                192/11250 0.0005773916 0.032911323
KRAS.300_UP.V1_DN             143/11250 0.0012983332 0.048535110
KRAS.LUNG.BREAST_UP.V1_DN     145/11250 0.0014191553 0.048535110
                                   qvalue
KRAS.600.LUNG.BREAST_UP.V1_DN 0.007763107
KRAS.600_UP.V1_DN             0.007763107
RPS14_DN.V1_UP                0.029578659
KRAS.300_UP.V1_DN             0.043620351
KRAS.LUNG.BREAST_UP.V1_DN     0.043620351
                                                                                              geneID
KRAS.600.LUNG.BREAST_UP.V1_DN 50489/8424/3848/338/623/8911/3963/79852/6332/6664/9844/64167/3779/3613
KRAS.600_UP.V1_DN              3872/50489/3848/338/57118/6236/8911/3963/6332/6664/3779/684/7704/3613
RPS14_DN.V1_UP                                    3128/9934/7045/7450/8942/64231/4680/6347/7940/7305
KRAS.300_UP.V1_DN                                             3848/338/57118/6236/8911/3963/6664/684
KRAS.LUNG.BREAST_UP.V1_DN                                     338/623/8911/3963/6332/6664/64167/3779
                              Count
KRAS.600.LUNG.BREAST_UP.V1_DN    14
KRAS.600_UP.V1_DN                14
RPS14_DN.V1_UP                   10
KRAS.300_UP.V1_DN                 8
KRAS.LUNG.BREAST_UP.V1_DN         8
```

``` r
em.up<-enricher(gene.up,TERM2GENE=m_t2g)
head(em.up)
```

``` 
                                 ID       Description GeneRatio   BgRatio
RPS14_DN.V1_UP       RPS14_DN.V1_UP    RPS14_DN.V1_UP     9/105 192/11250
KRAS.600_UP.V1_DN KRAS.600_UP.V1_DN KRAS.600_UP.V1_DN    10/105 289/11250
                        pvalue   p.adjust     qvalue
RPS14_DN.V1_UP    7.463861e-05 0.01171826 0.01044941
KRAS.600_UP.V1_DN 3.591599e-04 0.02819405 0.02514119
                                                             geneID Count
RPS14_DN.V1_UP        3128/9934/7045/7450/8942/64231/4680/7940/7305     9
KRAS.600_UP.V1_DN 3872/50489/3848/338/57118/3963/6664/3779/684/7704    10
```

``` r
em.down<-enricher(gene.down,TERM2GENE=m_t2g)
head(em.down)
```

``` 
                                                 ID
LEF1_UP.V1_DN                         LEF1_UP.V1_DN
HINATA_NFKB_IMMU_INF           HINATA_NFKB_IMMU_INF
GCNP_SHH_UP_EARLY.V1_DN     GCNP_SHH_UP_EARLY.V1_DN
KRAS.LUNG.BREAST_UP.V1_UP KRAS.LUNG.BREAST_UP.V1_UP
                                        Description GeneRatio   BgRatio
LEF1_UP.V1_DN                         LEF1_UP.V1_DN      7/62 190/11250
HINATA_NFKB_IMMU_INF           HINATA_NFKB_IMMU_INF      3/62  17/11250
GCNP_SHH_UP_EARLY.V1_DN     GCNP_SHH_UP_EARLY.V1_DN      6/62 169/11250
KRAS.LUNG.BREAST_UP.V1_UP KRAS.LUNG.BREAST_UP.V1_UP      5/62 145/11250
                                pvalue    p.adjust      qvalue
LEF1_UP.V1_DN             7.872657e-05 0.006924768 0.006587303
HINATA_NFKB_IMMU_INF      1.025892e-04 0.006924768 0.006587303
GCNP_SHH_UP_EARLY.V1_DN   3.223137e-04 0.014504116 0.013797288
KRAS.LUNG.BREAST_UP.V1_UP 1.190478e-03 0.040178643 0.038220620
                                                       geneID Count
LEF1_UP.V1_DN             4071/2947/2920/1080/4973/80736/2888     7
HINATA_NFKB_IMMU_INF                           2920/4973/6347     3
GCNP_SHH_UP_EARLY.V1_DN       4071/10461/80820/4929/6347/8788     6
KRAS.LUNG.BREAST_UP.V1_UP            1839/2920/7869/1135/3976     5
```

Immunologic collection: Immunologic signatures collection (also called
ImmuneSigDB) is composed of gene sets that represent cell types, states,
and perturbations within the immune system. The signatures were
generated by manual curation of published studies in human and mouse
immunology.

``` r
#Immunologic collection
m_t2g2 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g2)
```

    # A tibble: 6 x 2
      gs_name                             entrez_gene
      <chr>                                     <int>
    1 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN          20
    2 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN       10057
    3 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN       25864
    4 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN          34
    5 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN          54
    6 GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN       51205

``` r
em2<-enricher(gene,TERM2GENE=m_t2g2)
head(em2)
```

``` 
                                                                                                                                                                           ID
GSE10325_BCELL_VS_MYELOID_DN                                                                                                                     GSE10325_BCELL_VS_MYELOID_DN
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                                                         GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                                                   GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                                                                 GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                                                         GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP
                                                                                                                                                                  Description
GSE10325_BCELL_VS_MYELOID_DN                                                                                                                     GSE10325_BCELL_VS_MYELOID_DN
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                                                         GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                                                   GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                                                                 GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                                                         GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP
                                                                                       GeneRatio
GSE10325_BCELL_VS_MYELOID_DN                                                              16/229
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP    13/229
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                        13/229
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                     12/229
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                            11/229
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                                11/229
                                                                                         BgRatio
GSE10325_BCELL_VS_MYELOID_DN                                                           200/20652
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 200/20652
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                     200/20652
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                  200/20652
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                         199/20652
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                             200/20652
                                                                                             pvalue
GSE10325_BCELL_VS_MYELOID_DN                                                           8.612534e-10
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 3.872615e-07
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                     3.872615e-07
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                  2.532746e-06
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                         1.448152e-05
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                             1.517695e-05
                                                                                           p.adjust
GSE10325_BCELL_VS_MYELOID_DN                                                           3.836023e-06
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 5.749543e-04
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                     5.749543e-04
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                  2.820212e-03
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                         4.224883e-03
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                             4.224883e-03
                                                                                             qvalue
GSE10325_BCELL_VS_MYELOID_DN                                                           3.556523e-06
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 5.330621e-04
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                     5.330621e-04
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                  2.614727e-03
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                         3.917051e-03
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                             3.917051e-03
                                                                                                                                                                    geneID
GSE10325_BCELL_VS_MYELOID_DN                                                           1839/4610/64170/7045/978/729230/1234/64231/5265/2205/7940/5199/3779/1240/10462/3613
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP              3002/25939/433/4283/1234/219972/6039/6347/684/197259/53827/3135/441168
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                                        1839/4071/1728/4283/2920/4501/4973/1234/8942/1087/2180/9651/8651
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                                           6236/4502/3113/4501/8544/3301/978/8942/4496/2180/7305/3613
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                                                       3320/54206/2947/871/4283/4502/2920/3301/6446/3164/286
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                                                   2947/25939/9934/55283/4502/729230/2330/340348/1901/5199/55092
                                                                                       Count
GSE10325_BCELL_VS_MYELOID_DN                                                              16
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP    13
GSE360_HIGH_DOSE_B_MALAYI_VS_M_TUBERCULOSIS_MAC_DN                                        13
GSE24634_IL4_VS_CTRL_TREATED_NAIVE_CD4_TCELL_DAY10_DN                                     12
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN                                            11
GSE24142_EARLY_THYMIC_PROGENITOR_VS_DN2_THYMOCYTE_FETAL_UP                                11
```

``` r
em2.up<-enricher(gene.up,TERM2GENE=m_t2g2)
head(em2.up)
```

``` 
                                                                                                                                                                           ID
GSE10325_BCELL_VS_MYELOID_DN                                                                                                                     GSE10325_BCELL_VS_MYELOID_DN
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                                                                     GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                                                   GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                                                                               GSE2706_2H_VS_8H_LPS_STIM_DC_DN
GSE29618_MONOCYTE_VS_MDC_DN                                                                                                                       GSE29618_MONOCYTE_VS_MDC_DN
                                                                                                                                                                  Description
GSE10325_BCELL_VS_MYELOID_DN                                                                                                                     GSE10325_BCELL_VS_MYELOID_DN
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                                                                     GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                                                   GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                                                                               GSE2706_2H_VS_8H_LPS_STIM_DC_DN
GSE29618_MONOCYTE_VS_MDC_DN                                                                                                                       GSE29618_MONOCYTE_VS_MDC_DN
                                                                                       GeneRatio
GSE10325_BCELL_VS_MYELOID_DN                                                              14/141
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP    12/141
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                              11/141
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                             10/141
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                           10/141
GSE29618_MONOCYTE_VS_MDC_DN                                                               10/141
                                                                                         BgRatio
GSE10325_BCELL_VS_MYELOID_DN                                                           200/20652
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 200/20652
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                           200/20652
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                          200/20652
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                        200/20652
GSE29618_MONOCYTE_VS_MDC_DN                                                            200/20652
                                                                                             pvalue
GSE10325_BCELL_VS_MYELOID_DN                                                           9.943699e-11
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 1.307722e-08
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                           1.310619e-07
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                          1.190705e-06
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                        1.190705e-06
GSE29618_MONOCYTE_VS_MDC_DN                                                            1.190705e-06
                                                                                           p.adjust
GSE10325_BCELL_VS_MYELOID_DN                                                           3.865116e-07
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 2.541558e-05
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                           1.698125e-04
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                          5.785338e-04
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                        5.785338e-04
GSE29618_MONOCYTE_VS_MDC_DN                                                            5.785338e-04
                                                                                             qvalue
GSE10325_BCELL_VS_MYELOID_DN                                                           3.335849e-07
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP 2.193532e-05
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                           1.465593e-04
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                          4.993128e-04
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                        4.993128e-04
GSE29618_MONOCYTE_VS_MDC_DN                                                            4.993128e-04
                                                                                                                                                          geneID
GSE10325_BCELL_VS_MYELOID_DN                                                           4610/64170/7045/978/729230/1234/64231/5265/2205/7940/5199/3779/1240/10462
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP         3002/25939/433/4283/1234/219972/6039/684/197259/53827/3135/441168
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                                          3128/50489/3115/25939/9934/3113/3122/9844/2205/10462/51466
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                                               3002/3115/9934/3113/7045/729230/3122/64231/6039/7305
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                                        1236/84966/4502/4501/51365/100128590/8942/5265/23362/3135
GSE29618_MONOCYTE_VS_MDC_DN                                                                                3128/50489/3115/25939/3113/1234/3122/2205/10462/51466
                                                                                       Count
GSE10325_BCELL_VS_MYELOID_DN                                                              14
GSE19888_ADENOSINE_A3R_INH_PRETREAT_AND_ACT_BY_A3R_VS_TCELL_MEMBRANES_ACT_MAST_CELL_UP    12
GSE29618_MONOCYTE_VS_MDC_DAY7_FLU_VACCINE_DN                                              11
GSE17301_ACD3_ACD28_VS_ACD3_ACD28_AND_IFNA5_STIM_CD8_TCELL_DN                             10
GSE2706_2H_VS_8H_LPS_STIM_DC_DN                                                           10
GSE29618_MONOCYTE_VS_MDC_DN                                                               10
```

``` r
em2.down<-enricher(gene.down,TERM2GENE=m_t2g2)
head(em2.down)
```

``` 
                                                                                                   ID
GSE20715_0H_VS_48H_OZONE_LUNG_DN                                     GSE20715_0H_VS_48H_OZONE_LUNG_DN
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN         GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP                 GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP                 GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN                 GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN
                                                                                          Description
GSE20715_0H_VS_48H_OZONE_LUNG_DN                                     GSE20715_0H_VS_48H_OZONE_LUNG_DN
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN         GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP                 GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP                 GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN                 GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN
                                                   GeneRatio   BgRatio
GSE20715_0H_VS_48H_OZONE_LUNG_DN                        8/88 200/20652
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN          7/88 199/20652
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN      7/88 200/20652
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP              7/88 200/20652
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP              7/88 200/20652
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN              7/88 200/20652
                                                         pvalue
GSE20715_0H_VS_48H_OZONE_LUNG_DN                   2.228078e-06
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN     2.278456e-05
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN 2.353013e-05
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP         2.353013e-05
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP         2.353013e-05
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN         2.353013e-05
                                                      p.adjust      qvalue
GSE20715_0H_VS_48H_OZONE_LUNG_DN                   0.006608480 0.006529442
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN     0.009970054 0.009850811
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN 0.009970054 0.009850811
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP         0.009970054 0.009850811
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP         0.009970054 0.009850811
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN         0.009970054 0.009850811
                                                                                       geneID
GSE20715_0H_VS_48H_OZONE_LUNG_DN                   54206/3315/55466/871/27101/3301/3326/10598
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN             3320/54206/2947/871/2920/3301/3164
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN       54206/80031/6236/2920/4929/2180/3164
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP                  3320/1728/871/9229/4583/1080/2888
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP                 3320/54206/2947/871/2920/3301/3164
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN                 4071/4133/2947/1728/4929/6332/2180
                                                   Count
GSE20715_0H_VS_48H_OZONE_LUNG_DN                       8
GSE42021_CD24HI_TREG_VS_CD24HI_TCONV_THYMUS_DN         7
GSE26343_UNSTIM_VS_LPS_STIM_NFAT5_KO_MACROPHAGE_DN     7
GSE360_L_MAJOR_VS_B_MALAYI_LOW_DOSE_MAC_UP             7
GSE42021_TREG_PLN_VS_CD24HI_TREG_THYMUS_UP             7
GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN             7
```

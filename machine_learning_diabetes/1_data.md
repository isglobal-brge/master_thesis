
# Data

**Contents:**

- NA treatment
   - NA function

```R

require(data.table);require(parallel)


```

```R

gera1 <- fread("/scratch/dguinon/gds_model_GERA1.txt")
gera2 <- fread("/scratch/dguinon/gds_model_GERA2.txt")

gera1 <- gera1[,-1]
gera2 <- gera2[,-1]

```

## NA treatment

First look at the missing data in the dataset

```R

na_column1 <- sapply(gera1, function(x) sum(is.na(x))) #suma por columnas
sum(na_column1)

na_column2 <- sapply(gera2, function(x) sum(is.na(x))) #suma por columnas
sum(na_column2)

GERA1 <- fread("/scratch/dguinon/GERA1.sample")
GERA2 <- fread("/scratch/dguinon/GERA2.sample")

```
## NA function

NA imputation function:

```R

ff <- function(x, i){
  mean(x==i, na.rm=TRUE)
}

impute_function <- function(data){
  
  p0 <- mclapply(data, ff, i=0, mc.cores=1)
  p1 <- mclapply(data, ff, i=1, mc.cores=1)
  p2 <- mclapply(data, ff, i=2, mc.cores=1)
  
  p <- data.frame(p0=unlist(p0),
                  p1=unlist(p1),
                  p2=unlist(p2))
  
  NA_indices <- mclapply(data,  function(x) which(is.na(x)), mc.cores=1)
  
  #
  #replace missing genotypes by sampling from (0,1,2) based on probabilities given in table p
  #
  
  for (i in 1:length(NA_indices)){
    inds <- NA_indices[[i]]
    n <- length(inds)
    if (n!=0) {
      data[inds, i] <- sample(c(0:2), n,
                              replace = TRUE,
                              prob = p[i, ])
    }
  }
  return(data)
}

gera1_imputed <- impute_function(gera1)
gera2_imputed <- impute_function(gera2)

#setwd("/scratch/dguinon/")
write.table(gera1_imputed, "gera1_imputed_d.txt", row.names = FALSE, quote = FALSE)
write.table(gera2_imputed, "gera2_imputed_d.txt", row.names = FALSE, quote = FALSE)

```

Binding data and exporting final datasets

```R

gera1_def  <- cbind(gera1_imputed, GERA1$AFF)
names(gera1_def)[names(gera1_def)=="GERA1$AFF"] <- "affected"
gera2_def <- cbind(gera2_imputed, GERA2$AFF)
names(gera2_def)[names(gera2_def)=="GERA2$AFF"] <- "affected"

write.table(gera1_def, "gera1_def1.txt", row.names = FALSE, quote = FALSE)
write.table(gera2_def, "gera2_def1.txt", row.names = FALSE, quote = FALSE)

############################################################################################

gera1_imputed <- fread("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera/gera1_imputed_d.txt")
gera2_imputed <- fread("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera/gera2_imputed_d.txt")

names(gera1_imputed)[names(gera1_imputed)=="cascon"] <- "affected"
names(gera2_imputed)[names(gera2_imputed)=="cascon"] <- "affected"

setwd("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera/") # directory

write.table(gera1_imputed, "gera1_def_1.txt", row.names = FALSE, quote = FALSE)
write.table(gera2_imputed, "gera2_def_1.txt", row.names = FALSE, quote = FALSE)

```


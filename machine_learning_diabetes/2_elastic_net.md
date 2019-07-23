
# Elastic net 

**Contents:**

- Data wrangling
- Elastic net model
- Results


```R

require(data.table);require(biglasso)

``` 

## Data wrangling

Load the train and test sets

```R

gera1 <- fread("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera/gera1_def_1.txt")
gera2 <- fread("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera/gera2_def_1.txt")

gera1_x <- gera1[, -8574]
gera2_x <- gera2[, -8574]

```

Data need be converted into big.matrix format for biglasso implementation.

```R

X <- as.big.matrix(gera1_x)
dim(X)

Y <- as.matrix(gera1$affected)
dim(Y)

#############################

Z <- as.big.matrix(gera2_x)
dim(Z)

Z_y <- as.matrix(gera2$affected)

```

## Elastic net model

```R
time.cvfit <- system.time(
  cvfit1 <- cv.biglasso(X, Y, family = 'binomial', alpha = 0.5, penalty='enet', seed = 1234, nfolds = 10, ncores = 4)
)
print(time.cvfit)# to inspect the elapsed time

cvfit1$lambda.min # optimal lambda value
summary(cvfit1)

# coefficients obtained by the model
coefs<-as.matrix(coef(cvfit1))
coefs[which(coefs != 0)]
cvfit1$lambda.min

``` 

## Results

Plotting the model

```R

par(mfrow = c(2, 2), mar = c(3.5, 3.5, 3, 1), mgp = c(2.5, 0.5, 0))
plot(cvfit1, type = "all")

```

Predictions on the train and the test sets

```R

predict_train <- as.vector(predict(cvfit1, X, type="class", lambda= cvfit1$lambda.min))

table <- table(predict_train, Y)

accuracy <- sum(diag(table))/sum(table) #accuracy of the model

###################################################################################

predict_test <- as.vector(predict(cvfit1, Z, type="class", lambda= cvfit1$lambda.min))

table1 <- table(predict_test, Z_y)

accuracy1 <- sum(diag(table1))/sum(table1) #accuracy of the model

``` 

We can create a new train and a new test datasets with the selected variables

```R

variables_predicted <- predict(cvfit1, type = "vars", lambda= cvfit1$lambda.min)

#### select variables to create new training set

gera1_d <- as.data.frame(gera1_x)
gera_final1 <- data.frame(c(1:nrow(gera1)))
xx <- NULL

for(i in 1:length(variables_predicted)){
  xx <- as.numeric(variables_predicted[i]) 
  gera_final1[,i] <- gera1_d[, xx]
}

for(i in 1:length(variables_predicted)){
  xxx <- names(variables_predicted)[i] 
  colnames(gera_final1)[i] <- xx
}

dim(gera_final1)

#### select variables to create new test set

gera2_d <- as.data.frame(gera2_x)
gera_final2 <- data.frame(c(1:nrow(gera2)))
pp <- NULL

for(i in 1:length(variables_predicted)){
  pp <- as.numeric(variables_predicted[i]) 
  gera_final2[,i] <- gera2_d[, pp]
}

for(i in 1:length(variables_predicted)){
  ppp <- names(variables_predicted)[i] 
  colnames(gera_final2)[i] <- ppp
}

dim(gera_final2)

``` 

Finally the data is binded with the response vectors and exported for future usage

```R

gera_final1 <- cbind(gera_final1, gera1$affected)

gera_final2 <- cbind(gera_final2, gera2$affected)


names(gera_final1)[names(gera_final1)=="V2"] <- "affected"
names(gera_final2)[names(gera_final2)=="V2"] <- "affected"

setwd("C:/Users/okuti/Desktop/ESTADISTICA/tfg/deeplearning/diabetes/datos/gera") # Directory

write.table(gera_final1, "gera_final1.txt", row.names = FALSE, quote = FALSE)

write.table(gera_final2, "gera_final2.txt", row.names = FALSE, quote = FALSE)

```


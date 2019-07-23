
# XGBoost

**Contents:**

- Data
- Up sampling
- XGBoost model
- Model validation


```R

require(data.table);require(caret)


```

## Data

The response variable must be converted to a factor:

```R

gera1_def <- fread("/scratch/dguinon/gera_final1.txt")
gera2_def <- fread("/scratch/dguinon/gera_final2.txt")

gera1_def$affected <- as.factor(gera1_def$affected)
gera2_def$affected <- as.factor(gera2_def$affected)

gera1_def <- as.data.frame(gera1_def)
gera2_def <- as.data.frame(gera2_def)

``` 

## Up-sampling

Solving the imbalanced data problem with the up-sample function.

```R

set.seed(9560) # setting the seed, for it samples randomly
up_train <- upSample(x = gera1_def[, -ncol(gera1_def)],
                      y = gera1_def$affected)
table(up_train$Class) # This shows how the classes ended up
names(up_train)[names(up_train)=="Class"] <- "affected"

``` 

## XGBoost model

```R

data <- up_train
test <- gera2_def

y <- as.factor(data$affected)
data$affected <- y

test_response <- test$affected

trctrl <- trainControl(method = "cv", number = 10) # 10-fold cross validation

tune_grid <- expand.grid(nrounds = 200,
                         max_depth = c(3,6),
                         eta = 0.05,
                         gamma = 0.01,
                         colsample_bytree = c(0.5,0.75),
                         min_child_weight = 2,
                         subsample = 0.5) # tunning hyperparameters

#################################################################################

                           ###### XGBoost ######

rf_fit <- train(affected ~. -affected, data = data, method = "xgbTree",
                trControl=trctrl,
                tuneGrid = tune_grid,
                tuneLength = 10, linout = 0)
``` 

## Model validation
                
```R

test_predict <- predict(rf_fit, test)
(table1 <- table(test_response, test_predict))
(table1[1,1]+table1[2,2])/(sum(table1))


espe <- (table1[1,1]/(table1[1,1]+table1[1,2])) # Specificity

sens <- (table1[2,2]/(table1[2,2]+table1[2,1])) # Sensitivity

```



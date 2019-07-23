
# Deep learning

**Contents:**

- Data
- Tunning hyperparameters
- Final model


```R

require(data.table);require(h2o)


```

# Data

```R
gera1_def <- fread("/scratch/dguinon/gera_final1.txt")
gera2_def <- fread("/scratch/dguinon/gera_final2.txt")

h2o.init()

gera1_def <- as.data.frame(gera1_def)
gera2_def <- as.data.frame(gera2_def)

gera1_def$affected <- as.factor(gera1_def$affected)
gera2_def$affected <- as.factor(gera2_def$affected)

gera1_def <- as.h2o(gera1_def)
gera2_def <- as.h2o(gera2_def)

data <- gera1_def
test <- gera2_def

response <- "affected"
predictors <- setdiff(names(gera1_def), response)

```

# Tunning hyperparameters

```R
#set parameter space
activation_opt<- c("Rectifier", "Maxout", "Tanh", "RectifierWithDropout", "MaxoutWithDropout", "TanhWithDropout") 
hidden_opt <- list(c(50,50,50,50), c(100, 100, 50, 50, 50), c(100, 100, 100, 100, 100), c(500,500,500)) 
#layer configurations
l1_opt <- c(0,1e-3,1e-5)
l2_opt <- c(0,1e-3,1e-5)

hyper_params <- list( activation=activation_opt,
                      hidden=hidden_opt,
                      l1=l1_opt,
                      l2=l2_opt )

#set search criteria
search_criteria <- list(strategy = "RandomDiscrete", max_models=10)

#train model
dl_grid <- h2o.grid("deeplearning"
                    ,grid_id = "deep_learn"
                    ,hyper_params = hyper_params
                    ,search_criteria = search_criteria
                    ,training_frame = data
                    ,x=predictors
                    ,y=response
                    ,nfolds = 5
                    ,epochs = 100)

#get best model
d_grid <- h2o.getGrid("deep_learn",sort_by = "accuracy")
best_dl_model <- h2o.getModel(d_grid@model_ids[[1]])
h2o.performance(best_dl_model, xval = TRUE)

```

# Final model

```R

model1 <- h2o.deeplearning(
  model_id="dl_model_faster1", 
  training_frame=data,
  nfolds=10,
  x=predictors,
  y=response,
  hidden=c(500, 500, 500),                  
  epochs=1000000,                      
  score_validation_samples=10000,      
  fold_assignment = "Modulo",
  stopping_rounds=2,
  stopping_metric="misclassification",
  stopping_tolerance=0.01,
  balance_classes=TRUE,
  class_sampling_factors = c(1, 3),
  reproducible = TRUE,
  seed = 1
)
summary(model2)


pred1 <- h2o.performance(model = model1, 
                               newdata = test)

h2o.auc(pred1)
h2o.confusionMatrix(pred1)

```



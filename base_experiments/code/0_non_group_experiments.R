#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(doMC)
library(caret)
library(data.table)
library(randomForest)
library(rrBLUP)
library(gbm)
library(e1071)
library(kernlab)
library(xgboost)
library(foreach)
library(doRNG)
source("lib_general_functions.R")

FitBLUP <- function(df){
  set.seed(352)
  trait.response <- mixed.solve(df$y.train, 
                                Z = df$x.train,
                                K = NULL, 
                                SE = FALSE, 
                                return.Hinv = FALSE)
  e <- as.matrix(trait.response$u)
    
  test.trait.pred <- df$x.test %*% e
  test.blup.pred <- (test.trait.pred[, 1]) + trait.response$beta
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.blup.pred)
    
  return(test.perf.metrics$rsquared)
}

FitRF <- function(df, parallel_, ncores = NULL){
  ntrees = 1000

  if (parallel_ == TRUE){
    ncores = getDoParWorkers()
    set.seed(345)
    rf <- foreach(ntree = rep(floor(ntrees/ncores), ncores), .combine = combine, 
                  .multicombine = TRUE, .packages = 'randomForest') %dorng% {
                randomForest(x = df$x.train, y = df$y.train, ntree = ntree)
          }   
    test.rf.pred <- predict(rf, df$x.test)
    test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.rf.pred) 
    
    return(test.perf.metrics$rsquared)
  
  } else {
    set.seed(345)
    rf.fit <- randomForest(x = df$x.train, 
                           y = df$y.train,
                           ntree = ntrees)
        
    test.rf.pred <- predict(rf.fit, df$x.test)
    test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.rf.pred)
    
    return(test.perf.metrics$rsquared)
  }   
}

FitGBM <- function(df){
  set.seed(889)
  gbmfit <- gbm.fit(x = df$x.train, 
                    y = df$y.train, 
                    distribution = "gaussian",
                    n.trees = 1500,
                    interaction.depth = 20,
                    n.minobsinnode = 30,
                    shrinkage = 0.1,
                    verbose = FALSE)
  
  best.iter <- gbm.perf(gbmfit, plot.it = FALSE, method = "OOB")
  test.gbmpredict <- predict(gbmfit, 
                             as.data.frame(df$x.test), 
                             best.iter)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.gbmpredict)

  return(test.perf.metrics$rsquared)
}

FitSVM <- function(df, parallel_){
  ctrl <- trainControl(method = "cv", number = 5, allowParallel = parallel_)
  svmgrid <- expand.grid(sigma = 2^seq(-15, 1, by = 2),
                         C = c(0.1, 1, 10, 100))
  set.seed(326)
  svmfit <- train(x = df$x.train, 
                  y = df$y.train,
                  scaled = FALSE,
                  method = "svmRadial",
                  metric = "RMSE",
                  tuneGrid = svmgrid,
                  trControl = ctrl)
    
  test.svm.pred <- predict(svmfit$finalModel, df$x.test) 
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.svm.pred)
   
  return(test.perf.metrics$rsquared)
}

FitKNNReg <- function(df, parallel_){
  ctrl <- trainControl(method = "cv", number = 10, allowParallel = parallel_)
  knngrid <- expand.grid(k = 1:30)
  set.seed(34535)
  
  knnregfit <- train(x = df$x.train, 
                     y = df$y.train, 
                     method = "knn",
                     trControl = ctrl,
                     tuneGrid = knngrid)

  test.knn.pred <- predict(knnregfit$finalModel, df$x.test)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.knn.pred)

  return(test.perf.metrics$rsquared)
}

FitXGBoost <- function(df){
  if (getDoParRegistered() == TRUE){
      registerDoSEQ()
  }
  ctrl <- trainControl(method = "cv", number = 3, allowParallel = FALSE)
  xgbgrid <- expand.grid(nrounds = c(500, 1000, 1500, 2000),
                         max_depth = c(3, 5, 7, 10),
                         eta = 0.01,
                         gamma = c(0.05, 0.1),
                         subsample = 0.5,
                         colsample_bytree = 0.5,
                         min_child_weight = c(1, 5, 10))

  set.seed(47567)
  xgbfit <- train(x = df$x.train,
                  y = df$y.train, 
                  method = "xgbTree", 
                  verbose = 0,
                  eval_metric = "rmse",
                  trControl = ctrl,
                  tuneGrid = xgbgrid)

  test.xgb.pred <- predict(xgbfit$finalModel, df$x.test)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.xgb.pred)

  return(test.perf.metrics$rsquared)
}

FitMultipleModels <- function(df, parallel_){
  predictions <- c(rrBLUP = FitBLUP(df),
                   RF = FitRF(df, parallel_),
                   GBM = FitGBM(df),
                   SVM = FitSVM(df, parallel_),
                   KNN = FitKNNReg(df, parallel_),
                   XGB = FitXGBoost(df))

  return(predictions)
}

PerformExperiments <- function(trait.list, parallel_){
  genotype.data <- as.matrix(data.frame(fread(paste0("../genotype_data/",
                                                     "12k_ld_imputed.csv"), 
                                              header = TRUE), 
                                        row.names = 1)
                              )
  phenotype.data <- read.csv(file = "../phenotype_data/quantitative_traits.csv",
                             row.names = 1,
                             header = TRUE)
    
  results <- NULL
  for (trait.name in trait.list){
    trait.varieties <- GetTraitDetails(trait.name)
    x <- genotype.data[trait.varieties,]
    y <- phenotype.data[trait.varieties, trait.name]

    test.indices <- GenerateTestIndices(length(y))
    x.learn <- x[-test.indices, ]
    y.learn <- y[-test.indices]
    x.test <- x[test.indices, ]
    y.test <- y[test.indices]
        
    df <- list(x.train = x.learn,
               y.train = y.learn,
               x.test = x.test,
               y.test = y.test)

    if (parallel_ == TRUE){
        ncores <- detectCores()
        registerDoMC(cores = 10)
    }

    predictions <- FitMultipleModels(df, parallel_)
    write(c(trait.name, predictions), 
          ncolumn = length(c(trait.name, predictions)),
          append = TRUE, 
          file = "../output/non_group_base_results/results.txt")
    results <- rbind(results, predictions)   
  }
  
  rownames(results) <- trait.list
  write.csv(results,
            file = "../output/non_group_base_results/results.csv")
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD",
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", 
                "PLT_POST", "SDHT")

PerformExperiments(trait.list, TRUE)
#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(doMC)
library(caret)
library(data.table)
library(randomForest)
library(glmnet)
library(rrBLUP)
library(gbm)
library(e1071)
library(kernlab)
library(xgboost)
library(foreach)
library(doRNG)

GenerateTestIndices <- function(sample.number, percent.split = 0.25){
  sample.size <- floor(percent.split * sample.number)
  set.seed(1357)
  sample(seq_len(sample.number), size = sample.size)
}

GenerateCVFolds <- function(n, folds = 5){
  set.seed(13579)
  split(sample(1:n), rep(1:folds, length=n))
}

GetPerformanceMetrics <- function(y.true, model.pred){
  rsquared <- 1 - (sum((y.true - model.pred)^2) / sum((y.true - mean(y.true))^2))
  mse <- sum((y.true - model.pred)^2) / length(y.true)
  rmse <- sqrt(sum((y.true - model.pred)^2) / length(y.true))

  list(predictions = model.pred, rsquared = rsquared, mse = mse, rmse = rmse)
}

GetTraitDetails <- function(trait.name){
  trait.varieties <- read.table(file = paste0("../select_varieties/", 
                                              trait.name, ".txt"
                                             ),
                                sep = "\n", 
                                header = TRUE
                               )[,1]
  trait.varieties <- as.vector(trait.varieties)
}

LogFitModel <- function(trait, fold, chromosome, learner, fit.model){
  model.path <- paste0("../logs/fit_models/", trait, "/", learner, "_", fold, "_", chromosome, ".rds")
  saveRDS(fit.model, model.path)
}

LogMetaFeatures <- function(meta.features){
  trait <- meta.features$trait
  fold <- meta.features$fold
  chromosome <- meta.features$chromosome
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features
  y.test <- meta.features$test$y.test.act

  validation.path <- paste0("../logs/meta_features/", 
                            trait, "/validation_mfs_", 
                            fold, "_", chromosome, ".csv")
  test.path <- paste0("../logs/meta_features/", 
                      trait, "/test_mfs_", 
                      fold, "_", chromosome, ".csv")

  val.dat <- cbind(validation.ensemble, y.true = y.validation)
  test.dat <- cbind(test.ensemble, y.true = y.test)

  write.csv(val.dat, validation.path)
  write.csv(test.dat, test.path)

  base.validation.path <- paste0("../logs/meta_features/", 
                                 trait, "/validation_",
                                 chromosome, "_")
  base.test.path <- paste0("../logs/meta_features/", 
                           trait, "/test_",
                           chromosome, "_")

  val.mses <- c(fold, meta.features$validation$mses)
  val.rmses <- c(fold, meta.features$validation$rmses)
  val.rsquareds <- c(fold, meta.features$validation$rsquareds)

  test.mses <- c(fold, meta.features$test$mses)
  test.rmses <- c(fold, meta.features$test$rmses)
  test.rsquareds <- c(fold, meta.features$test$rsquareds)

  learners <- c("fold", names(test.rsquareds))

  if (fold == 1){     
    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rsquareds.txt"))
    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "mses.txt"))
    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rmses.txt"))

    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rsquareds.txt"))
    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "mses.txt"))
    write(learners, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rmses.txt"))

    write(val.rsquareds, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rsquareds.txt"))
    write(val.mses, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "mses.txt"))
    write(val.rmses, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rmses.txt"))

    write(test.rsquareds, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rsquareds.txt"))
    write(test.mses, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "mses.txt"))
    write(test.rmses, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rmses.txt"))
  
  } else {
    write(val.rsquareds, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rsquareds.txt"))
    write(val.mses, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "mses.txt"))
    write(val.rmses, append = TRUE, ncol = length(learners), 
          file = paste0(base.validation.path, "rmses.txt"))

    write(test.rsquareds, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rsquareds.txt"))
    write(test.mses, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "mses.txt"))
    write(test.rmses, append = TRUE, ncol = length(learners), 
          file = paste0(base.test.path, "rmses.txt"))
  }
}

FitBLUP <- function(df){
  set.seed(352)
  trait.response <- mixed.solve(df$y.train, 
                                Z = df$x.train,
                                K = NULL, 
                                SE = FALSE, 
                                return.Hinv = FALSE)
  LogFitModel(df$trait, df$fold, df$chromosome, "rrBLUP", trait.response)
  e <- as.matrix(trait.response$u)
    
  validation.trait.pred <- df$x.validation %*% e
  validation.blup.pred <- (validation.trait.pred[,1]) + trait.response$beta
  validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.blup.pred)
    
  test.trait.pred <- df$x.test %*% e
  test.blup.pred <- (test.trait.pred[,1]) + trait.response$beta
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.blup.pred)  
  perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)
  
  return(perf.metrics)
}

FitRF <- function(df, parallel_, ncores = NULL){
  ntrees = 1000

  if (parallel_ == TRUE){
    ncores = getDoParWorkers() / 2
    set.seed(345)
    rf <- foreach(ntree = rep(floor(ntrees/ncores), ncores), .combine = combine, 
                  .multicombine = TRUE, .packages = 'randomForest') %dorng% {
                  randomForest(x = df$x.train, y = df$y.train, ntree = ntree)
          }   
    LogFitModel(df$trait, df$fold, df$chromosome, "RandomForest", rf)

    validation.rf.pred <- predict(rf, df$x.validation)
    validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.rf.pred)
    test.rf.pred <- predict(rf, df$x.test)
    test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.rf.pred)
    perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)
    
    return(perf.metrics)
  
  } else {
    set.seed(345)
    rf.fit <- randomForest(x = df$x.train, 
                           y = df$y.train,
                           ntree = ntrees)
    LogFitModel(df$trait, df$fold, df$chromosome, "RandomForest", rf.fit)

    validation.rf.pred <- predict(rf.fit, df$x.validation)
    validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.rf.pred)
    test.rf.pred <- predict(rf.fit, df$x.test)
    test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.rf.pred)
    perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)

    return(perf.metrics)
  }   
}

FitGBM <- function(df){
  set.seed(889)
  gbmfit <- gbm.fit(x = df$x.train, 
                    y = df$y.train, 
                    distribution = "gaussian",
                    n.trees = 1500,
                    interaction.depth = 5,
                    n.minobsinnode = 30,
                    shrinkage = 0.1,
                    verbose = FALSE)
  LogFitModel(df$trait, df$fold, df$chromosome, "GBM", gbmfit)

  best.iter <- gbm.perf(gbmfit, plot.it = FALSE, method = "OOB")
  validation.gbmpredict <- predict(gbmfit, 
                                   as.data.frame(df$x.validation), 
                                   best.iter)
    
  validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.gbmpredict)
  test.gbmpredict <- predict(gbmfit, 
                             as.data.frame(df$x.test), 
                             best.iter)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.gbmpredict)
  perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)

  return(perf.metrics)
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
  LogFitModel(df$trait, df$fold, df$chromosome, "tSVM", svmfit)

  validation.svm.pred <- predict(svmfit$finalModel, df$x.validation) 
  validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.svm.pred)
  test.svm.pred <- predict(svmfit$finalModel, df$x.test) 
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.svm.pred)
  perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)
    
  return(perf.metrics)
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
  LogFitModel(df$trait, df$fold, df$chromosome, "tKNN", knnregfit)

  validation.knn.pred <- predict(knnregfit$finalModel, df$x.validation)
  validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.knn.pred)
  test.knn.pred <- predict(knnregfit$finalModel, df$x.test)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.knn.pred)
  perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)

  return(perf.metrics)
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
                         colsample_bytree = 0.5,
                         min_child_weight = c(1, 5, 10),
                         subsample = 0.5)

  set.seed(47567)
  xgbfit <- train(x = df$x.train,
                  y = df$y.train, 
                  method = "xgbTree",
                  verbose = 0,
                  eval_metric = "rmse",
                  trControl = ctrl,
                  tuneGrid = xgbgrid)
  LogFitModel(df$trait, df$fold, df$chromosome, "XGB", xgbfit)

  validation.xgb.pred <- predict(xgbfit$finalModel, df$x.validation)
  validation.perf.metrics <- GetPerformanceMetrics(df$y.validation, validation.xgb.pred)
  test.xgb.pred <- predict(xgbfit$finalModel, df$x.test)
  test.perf.metrics <- GetPerformanceMetrics(df$y.test, test.xgb.pred)
  perf.metrics <- list(validation = validation.perf.metrics, test = test.perf.metrics)

  return(perf.metrics)
}

FitMultipleModels <- function(df, parallel_){
  predictions <- list(rrBLUP = FitBLUP(df),
                      RF = FitRF(df, parallel_),
                      GBM = FitGBM(df),
                      SVM = FitSVM(df, parallel_),
                      KNN = FitKNNReg(df, parallel_),
                      XGB = FitXGBoost(df))

  meta.features <- list(validation = list(meta.features = NULL,
                                          rsquareds = NULL, 
                                          mses = NULL, 
                                          rmses = NULL),
                        test = list(meta.features = NULL,
                                    rsquareds = NULL, 
                                    mses = NULL, 
                                    rmses = NULL)
                       )
  
  for (learner in names(predictions)){
    for (dset in names(predictions[[learner]])){
      meta.features[[dset]]$mses <- c(meta.features[[dset]]$mses, 
                                      predictions[[learner]][[dset]]$mse
                                     )
      meta.features[[dset]]$rmses <- c(meta.features[[dset]]$rmses, 
                                       predictions[[learner]][[dset]]$rmse
                                      )
      meta.features[[dset]]$rsquareds <- c(meta.features[[dset]]$rsquareds, 
                                           predictions[[learner]][[dset]]$rsquared
                                          )
      meta.features[[dset]]$meta.features <- cbind(meta.features[[dset]]$meta.features, 
                                                   predictions[[learner]][[dset]]$predictions
                                                  ) 
    }   
  }

  for (dset in names(meta.features)){
    for (metric in names(meta.features[[dset]])){
      if (is.matrix(meta.features[[dset]][[metric]])){
        colnames(meta.features[[dset]][[metric]]) <- names(predictions)
      } else {
        names(meta.features[[dset]][[metric]]) <- names(predictions)
      } 
    }   
  }
    
  meta.features$validation[["y.valid.act"]] <- df$y.validation
  meta.features$test[["y.test.act"]] <- df$y.test
  meta.features[["trait"]] <- df$trait
  meta.features[["fold"]] <- df$fold
  meta.features[["chromosome"]] <- df$chromosome
  LogMetaFeatures(meta.features)
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
  all.markers <- colnames(genotype.data)

  for (trait.name in trait.list){
    for(chromosome in 1:12){
      chromosome.markers <- as.matrix(read.table(file = paste0("../all_snps/", chromosome, ".txt"),
                                      header = FALSE)
                                     )[,1]
      chr.markers <- intersect(all.markers, chromosome.markers)

      trait.varieties <- GetTraitDetails(trait.name)
      x <- genotype.data[trait.varieties, chr.markers]
      y <- phenotype.data[trait.varieties, trait.name]
        
      test.indices <- GenerateTestIndices(length(y))
      x.learn <- x[-test.indices, ]
      y.learn <- y[-test.indices]
      x.test <- x[test.indices, ]
      y.test <- y[test.indices]

      df <- list(trait = trait.name,
                 x.learn = x.learn, 
                 y.learn = y.learn, 
                 x.test = x.test, 
                 y.test = y.test)
            
      cv.indices <- GenerateCVFolds(length(df$y.learn)) 
      for (fold in 1:length(cv.indices)){
        if (parallel_ == TRUE){
          ncores <- detectCores()
          registerDoMC(cores = ncores)
        }

        validation.indices <- cv.indices[[fold]]
        newdf <- list(trait = df$trait,
                      fold = fold,
                      chromosome = chromosome,
                      x.train = df$x.learn[-validation.indices, ],
                      y.train = df$y.learn[-validation.indices],
                      x.validation = df$x.learn[validation.indices, ],
                      y.validation = df$y.learn[validation.indices],
                      x.test = df$x.test,
                      y.test = df$y.test)
              
        FitMultipleModels(newdf, parallel_)
      }           
    }   
  }
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD",
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", 
                "PLT_POST", "SDHT")

PerformExperiments(trait.list, TRUE)
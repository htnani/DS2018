#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(doMC)
library(caret)
library(quadprog)
library(glmnet)
library(pls)
library(gbm)

GetPerformanceMetrics <- function(y.true, model.pred){
  rsquared <- 1 - (sum((y.true - model.pred)^2) / sum((y.true - mean(y.true))^2))
  mse <- sum((y.true - model.pred)^2) / length(y.true)
  rmse <- sqrt(sum((y.true - model.pred)^2) / length(y.true))

  list(predictions = model.pred, rsquared = rsquared, mse = mse, rmse = rmse)
}

sLR <- function(meta.features){
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features

  train.df <- as.data.frame(cbind(validation.ensemble, y.validation))
  lm.fit <- lm(y.validation~. + 0, data = train.df)
  weights <- coef(lm.fit)

  validation.fhat <- rowSums(validation.ensemble %*% diag(weights))
  test.fhat <- rowSums(test.ensemble %*% diag(weights))
  fhats <- list(validation = validation.fhat, test = test.fhat)

  return(fhats)
}

GradientDescent <- function(meta.features, alpha = 0.000005, n.iters = 200000){
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features

  X <- validation.ensemble
  y <- y.validation

  cost <- function(X, y, theta){
    sum(((X %*% theta) - y)^2) / (2 * length(y))
  } 

  cost.history <- double(n.iters)
  theta.history <- list(n.iters)

  # One is added to accomodate for the intercept
  theta <- matrix(rep(0, ncol(X)), nrow = ncol(X))

  for (i in 1:n.iters){
    error <- (X %*% theta - y)
    delta <- t(X) %*% error / length(y)
    theta <- theta - alpha * delta
    cost.history[i] <- cost(X, y, theta)
    theta.history[[i]] <- theta
  }

  weights <- diag(theta[, 1])

  validation.fhat <- rowSums(validation.ensemble %*% weights)
  test.fhat <- rowSums(test.ensemble %*% weights)

  fhats <- list(validation = validation.fhat, test = test.fhat)
  
  return(fhats)
}

FitRR <- function(meta.features){
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features

  set.seed(998)
  ridge.fit <- cv.glmnet(validation.ensemble, 
                         y.validation, 
                         standardize = FALSE, 
                         alpha = 0,
                         lambda = 10^seq(10, -3, length = 100),
                         parallel = TRUE)
    
  validation.fhat <- predict(ridge.fit, 
                             s = ridge.fit$lambda.min, 
                             newx = validation.ensemble)
  test.fhat <- predict(ridge.fit, 
                       s = ridge.fit$lambda.min, 
                       newx = test.ensemble)

  fhats <- list(validation = validation.fhat, test = test.fhat)
  
  return(fhats)
}

KRLS <- function(meta.features){
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features

  ctrl <- trainControl(method = "cv", number = 10)
  krlsgrid <- expand.grid(sigma = 10^seq(0, 3),
                          lambda = seq(0,1,0.05)[-1])
  set.seed(34657)
  krlsfit <- train(x = validation.ensemble,
                   y = y.validation,
                   method = "krlsRadial",
                   print.level = 0,
                   trControl = ctrl,
                   tuneGrid = krlsgrid)

  validation.fhat <- krlsfit$finalModel$fitted[, 1]
  test.fhat <- (predict(krlsfit$finalModel, test.ensemble)$fit)[, 1]  
  fhats <- list(validation = validation.fhat, test = test.fhat)

  return(fhats)
}

FitPCR <- function(meta.features){
  validation.ensemble <- meta.features$validation$meta.features
  y.validation <- meta.features$validation$y.valid.act
  test.ensemble <- meta.features$test$meta.features

  train.df <- as.data.frame(cbind(validation.ensemble, y = y.validation))
  pcrfit <- pcr(y~., data = train.df, scale = FALSE, validation = "CV")

  modelperf <- MSEP(pcrfit)
  modelperf <- as.matrix(modelperf$val)
  comp.inx <- seq(1, 2 * ncol(validation.ensemble), 2) 
  modelperf <- modelperf[comp.inx][-1]
  minp <- min(modelperf)
  comp <- which(modelperf == minp)

  validation.fhat <- predict(pcrfit, validation.ensemble, ncomp = comp)
  test.fhat <- predict(pcrfit, test.ensemble, ncomp = comp) 
  fhats <- list(validation = validation.fhat, test = test.fhat)
  
  return(fhats)
}

# Blending Logic
# Get validation meta features per chromosome.
GetValidationSets <- function(trait, nfolds, chromosomes){
  all.folds.metafs <- NULL
  for (fold in 1:nfolds){
    fold.metafs <- NULL
    y.actual <- NULL
    for (chrNo in 1:chromosomes){
      path <- paste0("../logs/meta_features/", trait, 
                     "/validation_mfs_", fold, "_", chrNo, ".csv")
      metafs <- read.csv(file = path, header = TRUE, row.names = 1)
      metafs <- as.matrix(metafs)
      learner.count <- ncol(metafs) - 1
      # The actual value is the last column, so we remove it. 
      if (is.null(y.actual)){
        y.actual <- metafs[, learner.count + 1]
      }           
      metafs <- metafs[, 1:learner.count]
      learners <- colnames(metafs)
      learner.chrs <- sapply(learners, 
                             function(learner)paste0(learner, "_", chrNo), 
                             USE.NAMES = FALSE)
      colnames(metafs) <- learner.chrs
      fold.metafs <- cbind(fold.metafs, metafs)
    }
    all.folds.metafs[[fold]] <- cbind(fold.metafs, y.actual)
  }

  all.folds.mses <- NULL
  for (chrNo in 1:chromosomes){
    mpath <- paste0("../logs/meta_features/", trait, "/validation_",
                    chrNo, "_mses.txt")
    mses <- read.table(mpath, header = TRUE, row.names = 1)
    mses <- as.matrix(mses)
    learners <- colnames(mses)
    learner.chrs <- sapply(learners, 
                           function(learner)paste0(learner, "_", chrNo), 
                           USE.NAMES = FALSE)
    colnames(mses) <- learner.chrs
    all.folds.mses <- cbind(all.folds.mses, mses)
  }

  results <- list(mses = all.folds.mses,
                  metafs = all.folds.metafs)

  return(results)
}

# Step 1B
# Get Performance of stacked learner predictions
GetStackedPerformances <- function(all.stacked.metafs, y.actual){
  chromosomes <- length(all.stacked.metafs)
  learner.performances <- list()
  for (chrNo in 1:chromosomes){
    stacked.metafs <- all.stacked.metafs[[chrNo]]
    learners <- colnames(stacked.metafs)
    for (learner in learners){
      predictions <- stacked.metafs[, learner]
      metrics <- GetPerformanceMetrics(y.actual, predictions)
      learner.performances[[learner]] <- c(learner.performances[[learner]],
                                           metrics$rsquared
                                          )
    }
  }
  learner.performances <- do.call(cbind, learner.performances)

  return(learner.performances)
}

# Stacks all predictions for each fold for each learner with averaging.
CreateSuperMetafs <- function(unstacked.metafs){
  folds <- length(unstacked.metafs)
  learner.predictions <- NULL
  for (fold in 1:folds){
    metafs <- unstacked.metafs[[fold]]
    learners <- colnames(metafs)
    for (learner in learners){
      learner.prediction <- metafs[, learner]
      learner.predictions[[learner]] <- cbind(learner.predictions[[learner]], 
                                              learner.prediction,
                                              deparse.level = 0)
    }
  }

  learners <- names(learner.predictions)
  stacked.metafs <- NULL
  for (learner in learners) {
    # Stack all features for a particular learner with averaging
    stacked.feature <- rowMeans(learner.predictions[[learner]])
    stacked.metafs <- cbind(stacked.metafs, 
                            stacked.feature,
                            deparse.level = 0)
  }
  colnames(stacked.metafs) <- learners

  return(stacked.metafs)
}

# Step 1A
BuildSuperStacksPerChr <- function(trait, nfolds, chromosomes){
  y.actual <- NULL
  all.stacked.metafs <- NULL
  for (chrNo in 1:chromosomes){
    unstacked.metafs <- NULL 
    for (fold in 1:nfolds){
      # Read all test meta-features. 
      path <- paste0("../logs/meta_features/", trait, "/test_mfs_",
                     fold, "_", chrNo, ".csv")
      metafs <- read.csv(file = path, header = TRUE, row.names = 1)
      metafs <- as.matrix(metafs)
      learner.count <- ncol(metafs) - 1
      # The actual value is the last column, so we remove it. 
      if (is.null(y.actual)){
        y.actual <- metafs[, learner.count + 1]
      }           
      metafs <- metafs[, 1:learner.count]
      unstacked.metafs[[fold]] <- metafs
    }
  
    stacked.metafs <- CreateSuperMetafs(unstacked.metafs)
    learners <- colnames(stacked.metafs)
    learner.chrs <- sapply(learners, 
                           function(learner)paste0(learner, "_", chrNo), 
                           USE.NAMES = FALSE)
    colnames(stacked.metafs) <- learner.chrs
    all.stacked.metafs <- cbind(all.stacked.metafs, stacked.metafs)
  }
  
  results <- list(y.actual = y.actual, 
                  stacked.metafs = all.stacked.metafs)

  return(results)
}

# Step 2
# Learn weights using each validation set for each chromosome in isolation
StackPredsWithWeights <- function(trait, nfolds, chromosomes){
  super.results <- BuildSuperStacksPerChr(trait, nfolds, chromosomes)
  validation.results <- GetValidationSets(trait, nfolds, chromosomes)
  validation.metafs <- validation.results$metafs
  validation.mses <- validation.results$mses
  stacked.mfs <- super.results$stacked.metafs # Per Chromosome
  y.actual <- super.results$y.actual # For test predictions
  # The stacked.perfs above are the performances of each level 0 learner
  # on a per chromosome basis. 

  folds.results <- NULL
  for (fold in 1:nfolds){
    meta.features <- list(validation = list(meta.features = NULL, 
                                            y.valid.act = NULL, 
                                            mses = NULL), 
                          test = list(meta.features = NULL, 
                                      y.test.act = NULL),
                          trait = trait, fold = NULL)

    valid.metafs <- validation.metafs[[fold]] 
    valid.metafs <- as.matrix(valid.metafs)
    learner.count <- ncol(valid.metafs) - 1
    # The actual value is the last column, so we remove it. 
    y.valid.act <- valid.metafs[, learner.count + 1]
    valid.metafs <- valid.metafs[, 1:learner.count]

    meta.features$validation$meta.features <- valid.metafs
    meta.features$validation$y.valid.act <- y.valid.act
    meta.features$validation$mses <- unlist(validation.mses[fold, ])
    
    meta.features$test$meta.features <- stacked.mfs
    meta.features$test$y.test.act <- y.actual

    # Stack predictions with multiple learners
    nf.hat <- list(LR = sLR(meta.features),
                   GD = GradientDescent(meta.features),
                   RR = FitRR(meta.features),
                   KRLS = KRLS(meta.features),
                   PCR = FitPCR(meta.features))
        
    for (combiner in names(nf.hat)){
      folds.results[[combiner]] <- cbind(folds.results[[combiner]], 
                                         nf.hat[[combiner]]$test)
    }
  }
  results <- list(folds.results = folds.results, y.actual = y.actual)

  return(results)
}

# Step 3. Aggregate final predictions. 
AggregateFinalPredictions <- function(trait, nfolds = 5, chromosomes = 12){
  results <- StackPredsWithWeights(trait, nfolds, chromosomes)
  folds.predictions <- results$folds.results
  y.actual <- results$y.actual

  final.rsquareds <- NULL
  learners <- names(folds.predictions)
  for (learner in learners){
    learner.predictions <- rowMeans(folds.predictions[[learner]])
    res <- GetPerformanceMetrics(y.actual, learner.predictions)
    final.rsquareds <- c(final.rsquareds, res$rsquared)
  }
  names(final.rsquareds) <- learners

  return(final.rsquareds)
} 

RunExperiments <- function(trait.list, parallel_ = TRUE){
  all.results <- NULL
  for (trait in trait.list){
    if (parallel_ == TRUE){
      ncores <- detectCores()
      registerDoMC(cores = ncores)
    }

    trait.results <- AggregateFinalPredictions(trait)
    all.results <- rbind(all.results, trait.results)
  }
  rownames(all.results) <- trait.list

  write.csv(t(all.results), "../logs/results/frameworks_b_results.csv")
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD",
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", 
                "PLT_POST", "SDHT")

RunExperiments(trait.list)
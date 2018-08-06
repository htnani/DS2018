#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(doMC)
library(caret)
library(quadprog)
library(glmnet)
library(pls)
library(FSelector)

GetPerformanceMetrics <- function(y.true, model.pred){
  rsquared <- 1 - (sum((y.true - model.pred)^2) / sum((y.true - mean(y.true))^2))
  mse <- sum((y.true - model.pred)^2) / length(y.true)
  rmse <- sqrt(sum((y.true - model.pred)^2) / length(y.true))

  list(predictions = model.pred, rsquared = rsquared, mse = mse, rmse = rmse)
}

LogMetaPerformance <- function(meta.performance){
  trait <- meta.performance$trait
  fold <- meta.performance$fold 
  validation <- meta.performance$validation
  test <- meta.performance$test

  base.validation.path <- paste0("../logs/meta_performance/validation_", trait, "_")
  base.test.path <- paste0("../logs/meta_performance/test_", trait, "_")

  val.rsquareds <- c(fold, validation$rsquareds)
  val.mses <- c(fold, validation$mses)
  val.rmses <- c(fold, validation$rmses)

  test.rsquareds <- c(fold, test$rsquareds)
  test.mses <- c(fold, test$mses)
  test.rmses <- c(fold, test$rmses)
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

#' Implementation of Linear Regression without intercept using lm
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

  weights <- diag(theta[,1])

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

  validation.fhat <- krlsfit$finalModel$fitted[,1]
  test.fhat <- (predict(krlsfit$finalModel, test.ensemble)$fit)[,1]
    
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

BuildMetaFeatures <- function(trait, nfolds = 5, parallel_ = TRUE){
  if (parallel_ == TRUE){
    ncores <- detectCores()
    registerDoMC(cores = ncores)
  }   

  # Average and merge test meta features.
  learner.features <- NULL
  y.true <- NULL
  for (fold in 1:nfolds){
    test.path <- paste0("../logs/meta_features/", trait, "/test_mfs_", fold, ".csv")
    test.meta.features <- as.matrix(read.csv(file = test.path, header = TRUE, row.names = 1))
    learners <- colnames(test.meta.features)

    if (is.null(y.true)){
      y.true <- test.meta.features[, learners[length(learners)]]
    }

    learners <- learners[-length(learners)]
    for (learner in learners){
      learner.features[[learner]] <- cbind(learner.features[[learner]], test.meta.features[,learner])
    }
  }

  merged.test.mfs <- NULL
  for (learner in names(learner.features)){
    ldata <- rowMeans(learner.features[[learner]])
    merged.test.mfs <- cbind(merged.test.mfs, ldata)
  }
  colnames(merged.test.mfs) <- names(learner.features)

  merged.test.perf <- list(rsquareds = NULL, rmses = NULL, mse = NULL)
  for (learner in colnames(merged.test.mfs)){
    res <- GetPerformanceMetrics(y.true, merged.test.mfs[, learner])
    merged.test.perf$rsquareds <- c(merged.test.perf$rsquareds, res$rsquared)
    merged.test.perf$rmses <- c(merged.test.perf$rmses, res$rmse)
    merged.test.perf$mse <- c(merged.test.perf$mse, res$mse) 
  }
  for (metric in names(merged.test.perf)){
    names(merged.test.perf[[metric]]) <- colnames(merged.test.mfs)
  }

  # Handle validation meta features.
  validation.mses.path <- paste0("../logs/meta_features/", trait, "/validation_mses.txt")
  validation.rsquareds.path <- paste0("../logs/meta_features/", trait, "/validation_rsquareds.txt")
  validation.mses <- (read.table(file = validation.mses.path, header = TRUE, row.names = 1))
  validation.rsquareds <- read.table(file = validation.rsquareds.path, header = TRUE, row.names = 1)

  meta.features <- list(validation = list(meta.features = NULL, 
                                          y.valid.act = NULL, 
                                          mses = NULL), 
                        test = list(meta.features = NULL, 
                                    rsquareds = NULL,
                                    rmses = NULL, 
                                    y.test.act = NULL),
                        trait = trait, 
                        fold = NULL
                       )

  results <- NULL
  for (fold in 1:nfolds){
    validation.path <- paste0("../logs/meta_features/", trait, "/validation_mfs_", fold, ".csv")
    validation.meta.features <- as.matrix(read.csv(file = validation.path, header = TRUE, row.names = 1))

    lcount <- ncol(validation.meta.features) 
    meta.features$validation$meta.features <- validation.meta.features[, 1:(lcount - 1)]
    meta.features$validation$y.valid.act <- validation.meta.features[, lcount]
    meta.features$validation$mses <- unlist(validation.mses[fold, ])

    meta.features$test$meta.features <- merged.test.mfs
    meta.features$test$y.test.act <- y.true
    meta.features$test$rmses <- merged.test.perf$rmses
    meta.features$test$rsquareds <- merged.test.perf$rsquareds

    meta.features$fold = fold

    nf.hat <- list(LR = sLR(meta.features),
                   GD = GradientDescent(meta.features),
                   RR = FitRR(meta.features),
                   KRLS = KRLS(meta.features),
                   PCR = FitPCR(meta.features)
                  )
    for (combiner in names(nf.hat)){
      results[[combiner]] <- cbind(results[[combiner]], nf.hat[[combiner]]$test)
    }
  }

  combiners.preds <- NULL
  for (combiner in names(results)){
    all.combiner.df <- results[[combiner]]
    combiner.mean.preds <- rowMeans(all.combiner.df)
    combiners.preds <- cbind(combiners.preds, combiner.mean.preds) 
  }
  colnames(combiners.preds) <- names(results)
  
  combiners.perf <- list(rsquareds = NULL, rmses = NULL, mse = NULL)
  for (learner in colnames(combiners.preds)){
    res <- GetPerformanceMetrics(y.true, combiners.preds[, learner])
    combiners.perf$rsquareds <- c(combiners.perf$rsquareds, res$rsquared)
    combiners.perf$rmses <- c(combiners.perf$rmses, res$rmse)
    combiners.perf$mse <- c(combiners.perf$mse, res$mse) 
  }
  for (metric in names(combiners.perf)){
    names(combiners.perf[[metric]]) <- colnames(combiners.preds)
  }
  
  all.test.rsquareds <- round(combiners.perf$rsquareds, 4)
  
  return(all.test.rsquareds)
}

CombineRegressions <- function(trait.list){
  rsquareds <- NULL
  for (trait in trait.list){
    res <- BuildMetaFeatures(trait)
    rsquareds <- rbind(rsquareds, res)    
  }
  rownames(rsquareds) <- trait.list
  write.csv(t(rsquareds), "../logs/results/framework_a_results.csv")
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD",
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", 
                "PLT_POST", "SDHT")

CombineRegressions(trait.list)
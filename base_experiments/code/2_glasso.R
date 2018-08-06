#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(foreach)
library(doMC)
library(doRNG)
library(caret)
library(data.table)
library(gglasso)
source("lib_general_functions.R")

FitGlasso <- function(df, group){
  set.seed(34834) 
  grlasso <- cv.gglasso(x = df$x.train, 
                        y = df$y.train, 
                        group = group, 
                        loss = "ls",
                        pred.loss = "L1", 
                        nfolds = 5)
  grlasso.fit <- grlasso$gglasso.fit
  best.lambda <- grlasso$lambda.min

  preds <- predict(grlasso.fit, 
                   newx = df$x.test, 
                   s = best.lambda)
  rsquared <- GetPerformanceMetrics(df$y.test, preds)$rsquared
    
  return(rsquared)
}

PerformExperiments <- function(trait.list, parallel_){
	registerDoMC(cores = 6)
  
	genotype.data <- as.matrix(data.frame(fread(paste0("../genotype_data/",
                                                     "12k_ld_imputed.csv"), 
                                              header = TRUE), 
                                        row.names = 1)
                              )
  phenotype.data <- read.csv(file = "../phenotype_data/quantitative_traits.csv",
                             row.names = 1,
                             header = TRUE)
  dataset.markers <- colnames(genotype.data) 
  group <- GetChromosomeGroupings(dataset.markers)
  results <- NULL

  # We grouped the execution of the trait experiments into two.
  # We first run the first six then the second. 
  # This can be easily modified for a different execution environment.
  exp.groups <- list(a = 1:6,
                     b = 7:12)	 
  all.results <- NULL
  for (exp.group in names(exp.groups)){
    grp <- exp.groups[exp.group] 
    results <- foreach(i = as.vector(unlist(grp))) %dopar% {
                 trait.name <- trait.list[i] 
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
                
                 result <- FitGlasso(df, group)
                 trait.result <- c(trait.name, result)
                 
                 write(trait.result,
                       append = TRUE, 
                       ncolumns = length(trait.result),
                       file = "../output/group_results/grlasso_results.txt")
                 
                 trait.result
             }
      all.results <- c(all.results, results)
    }
    
    results <- NULL
    for (i in 1:length(all.results)){
      trait.results <- all.results[[i]] 
      results <- rbind(results, trait.results)
  }

  rownames(results) <- NULL
  write.csv(results, file = "../output/group_results/grlasso_results.csv") 
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD", 
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", "PLT_POST",
                "SDHT")

PerformExperiments(trait.list, TRUE)
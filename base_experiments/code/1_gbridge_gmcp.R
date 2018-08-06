#!/usr/bin/Rscript --vanilla

require(methods)
library(parallel)
library(doMC)
library(caret)
library(data.table)
library(grpreg)
source("lib_general_functions.R")

FitGmcp <- function(df, group){
  grmcp.fit <-  grpreg(X = df$x.train,
                       y = df$y.train,
                       group = group,
                       penalty = "cMCP",
                       family = "gaussian",
                       lambda.min = 0.05,
                       eps = 1e-8,
                       max.iter = 3e+8)
  best.lambda <- select(grmcp.fit, criterion="AIC")$lambda

  preds <- predict(grmcp.fit, 
                   df$x.test, 
                   lambda = best.lambda)
  rsquared <- GetPerformanceMetrics(df$y.test, preds)$rsquared
        
  return(rsquared)
}

FitGbridge <- function(df, group){
  grbridge.fit <-  gBridge(X = df$x.train,
                           y = df$y.train,
                           group = group, 
                           family = "gaussian",
                           lambda.min = 0.05,
                           eps = 1e-8,
                           max.iter = 3e+8,
                           gamma = 0.95)
  best.lambda <- select(grbridge.fit, criterion="AIC")$lambda

  preds <- predict(grbridge.fit, 
                   df$x.test, 
                   lambda = best.lambda)
  rsquared <- GetPerformanceMetrics(df$y.test, preds)$rsquared
     
  return(rsquared)
}

FitMultipleModels <- function(df, group){
	predictions <- c(gbridge = FitGbridge(df, group),
                   gmcp = FitGmcp(df, group))

  return(predictions)
}

PerformExperiments <- function(trait.list){
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
    
    predictions <- FitMultipleModels(df, group)
    write(c(trait.name, predictions), 
          ncolumn = length(c(trait.name, predictions)),
          append = TRUE, 
          file = "../output/group_base_results/gbridge_gmcp_results.txt")
    results <- rbind(results, predictions)
	}
  
  rownames(results) <- trait.list
  write.csv(results,
            file = "../output/group_base_results/gbridge_gmcp_results.csv")
}

trait.list <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD", 
                "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", "PLT_POST",
                "SDHT")

PerformExperiments(trait.list)
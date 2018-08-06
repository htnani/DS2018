GenerateTestIndices <- function(sample.number, percent.split = 0.25){
  sample.size <- floor(percent.split * sample.number)
  set.seed(1357)
  sample(seq_len(sample.number), size = sample.size)
}

GetPerformanceMetrics <- function(y.true, model.pred){
  rsquared <- 1 - (sum((y.true - model.pred)^2) / sum((y.true - mean(y.true))^2))
  mse <- sum((y.true - model.pred)^2) / length(y.true)
  rmse <- sqrt(sum((y.true - model.pred)^2) / length(y.true))

  list(predictions = model.pred, rsquared = rsquared, mse = mse, rmse = rmse)
}

GetChromosomeGroupings <- function(dataset.markers){
  group <- NULL
  for (chr in 1:12){
    chromosome.markers <- as.matrix(read.table(file = 
                                               paste0("../all_snps/",
                                                      chr, ".txt"), 
                                               header = FALSE)
                                   )[, 1]  
    dataset.chr.markers <- intersect(chromosome.markers, dataset.markers) 
    chr.group <- rep(chr, length(dataset.chr.markers))
    group <- c(group, chr.group)
  } 

  return(group)
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
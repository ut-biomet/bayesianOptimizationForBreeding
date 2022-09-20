# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# The optimization algorithm create a one `.rds` per optimization repetition. This script aggregates all the results in one `.rds` file.

#### PACKAGES ####
library(pbapply)
library(dplyr)

#### CODE ####

# 1. Inputs ----
# where to save the "aggregated results" file
outputFolder <- 'aggregatedResults/'

# location of the optimizations' results files
# inputFolder <- 'output-resultRepetition/2022-02-15_16-32-29/'
inputFolder <- 'output-resultRepetition/1024_reps_15iters_excl/'

# name of the "aggregated results" file
# aggregatedFileName <- '16x32_repetitions'
# aggregatedFileName <- '1024x1_repetitions'
aggregatedFileName <- '1024-15itersx1_repetitions'


# 2. Aggregation ----

# create new name if output file already exists
aggregatedFile <- paste0(outputFolder, aggregatedFileName, '.rds')

i <- 2
while (file.exists(aggregatedFile)) {
  aggregatedFile <- paste0(outputFolder,
                           aggregatedFileName,
                           '_', i,
                           '.rds')
  i <- i + 1
}

resFiles <- list.files(inputFolder, full.names = TRUE)
resFiles <- resFiles[grep(pattern = '^.*\\/repRes_(randOpt|bayesOpt)_.*\\.rds$',
                  resFiles)]


aggregatedResultList <- pblapply(resFiles, function(f) {
  x <- readRDS(f)
  x$optResults$mboRawRes <- NULL # take too much memory
  x$optResults$simSetup$objFun <- NULL # take too much memory
  x$optResults$simSetup$fixedParams$initPop <- NULL # take too much memory
  x$optResults$simSetup$fixedParams$trait <- NULL # take too much memory
  x$optResults$simSetup$fixedParams$phenotyper <- NULL # take too much memory

  # add id
  if (unique(x$results$opMethod) == 'Random Optimization') {
    m <- 'ro'
  } else {m <- 'bo'}
  x$id <- paste0(m, '_', x$mainSeed)
  x
}, cl = 10)


names(aggregatedResultList) <- sapply(aggregatedResultList,
                                      function(x) x$id)


# Summary of the results
n_file <- length(aggregatedResultList)
scenarios <- as.data.frame(t(sapply(aggregatedResultList, function(x){
  c(he = x$optResults$simSetup$setupInfo$he,
    nGen = x$optResults$simSetup$fixedParams$nGen,
    budget = x$optResults$simSetup$fixedParams$budget,
    optMethod = x$optResults$resultType,
    n_resEval = x$params$nRep)
})))
scenarios <- data.frame(scenarios)
scenarios$budget <- as.numeric(scenarios$budget)
scenarios$nGen <- as.numeric(scenarios$nGen)
scenarios$genB <- scenarios$budget / scenarios$nGen

scenarios %>%
  group_by(he, nGen, genB, optMethod, n_resEval) %>%
  summarise(n_repeted_opt = n())


# save aggregated results
saveRDS(aggregatedResultList,
        aggregatedFile)

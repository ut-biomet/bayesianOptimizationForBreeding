# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
#




singleSimulation <- function(
  mainSeed,
  nCpus,
  nRep,
  i,
  iHomo,
  bRep,
  phenoFreq,
  setup,
  fp){

  #### SOURCES functions ####
  startTime <- Sys.time()

  set.seed(mainSeed)


  fun <- function(x) {
    # extract optimized parameters values
    i = x[1]
    iHomo = x[2]
    bRep = x[3]
    phenoFreq = x[4]

    # get simulation parameters
    newParams <- list(i = i,
                      iHomo = iHomo,
                      bRep = bRep,
                      phenoFreq = phenoFreq,
                      budget = fp$budget,
                      nGen = fp$nGen,
                      nIndIni = fp$initPop$nInd,
                      plotCost = fp$plotCost,
                      newIndCost = fp$newIndCost)
    params <- do.call(getSimulParams, newParams)

    # simulation
    finalGVs <- simuleBreeding(nPheno = params$nPheno,
                               nSelected = params$nSelected,
                               nNew = params$nNew,
                               nGen = fp$nGen,
                               initPop = fp$initPop,
                               trait = fp$trait,
                               phenotyper = fp$phenotyper,
                               createModel = fp$createModel,
                               selectInds = fp$selectInds,
                               matInds = fp$matInds,
                               verbose = FALSE)

    results <- fp$aggrFun(finalGVs)
    results
  }

  # 5 parallel simulations -------------------------------------------------------


  print(getSimulParams(i = i,
                 iHomo = iHomo,
                 bRep = bRep,
                 phenoFreq = phenoFreq,
                 budget = fp$budget,
                 nGen = fp$nGen,
                 nIndIni = fp$initPop$nInd,
                 plotCost = fp$plotCost,
                 newIndCost = fp$newIndCost))


  rep.x = replicate(nRep, c(i, iHomo, bRep, phenoFreq), simplify = FALSE)
  res <- unlist(mclapply(rep.x, fun, mc.cores = nCpus))

  results <- data.frame(i = i,
                    iHomo = iHomo,
                    bRep = bRep,
                    phenoFreq = phenoFreq,
                    res)
  colnames(results)[5] <- paste0("BV_", fp$aggrFunName)


  # add fixed parameters to results
  for (i in names(fp)) {
    if (is.numeric(fp[[i]]) | is.character(fp[[i]])) {
      results <- cbind(fp[[i]], results)
      colnames(results)[1] <- i
    }
  }

  # add simulation setup parameters to results
  for (i in names(setup)) {
    if (is.numeric(setup[[i]]) | is.character(setup[[i]])) {
      results <- cbind(setup[[i]], results)
      colnames(results)[1] <- i
    }
  }

  # create
  output <- list(results = results,
                 fp = fp,
                 setup = setup,
                 startTime = startTime)
  output

}

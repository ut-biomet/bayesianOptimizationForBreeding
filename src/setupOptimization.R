# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of function setting up the optimization

#' Setup optimization
#'
#' This function create a list
#'
#' @param simSetup_dataFile genetic data (vcf file) of the initial population.
#' @param simSetup_nQTN number of QTN used for the simulation
#' @param simSetup_nSNP number of SNP used for the simulation
#' @param simSetup_he heritability of the phenotypic trait for the initial population
#' @param simSetup_setupSeed random seed for the simulation setup
#' @param funEval number of objective function evaluation at each optimization cycle, in order to lower the noise of the objective function. The mean of the results of these evaluations will be return to the optimizer.
#' @param nCpusFunEval number of cpu to use when evaluating in parallel the objective function at each optimization cycle
#' @param nGen total number of generation
#' @param plotCost cost for phenotyping one plot
#' @param newIndCost cost for creating one plot
#' @param totalBudget total budget for the breeding campaign
#'
#' @return list containing all the information necessary for the optimization
#'
#' @examples
setupOpt <- function(simSetup_dataFile,
                     simSetup_nQTN,
                     simSetup_nSNP,
                     simSetup_he,
                     simSetup_setupSeed = NULL,
                     funEval,
                     nCpusFunEval,
                     nGen,
                     plotCost,
                     newIndCost,
                     totalBudget) {


  # 1 Setup simulation -----------------------------------------------------------
  # Create R objects necessary for the simulation using the package breedSimulatR

  cat("Calculate data MD5 fingerprint...\n")
  dataID <- paste0(simSetup_dataFile, "_",
                   digest(simSetup_dataFile, file = TRUE))


  setupSim <- setupSimul(data = simSetup_dataFile,
                         lchrCm = 120,
                         nQTN = simSetup_nQTN,
                         nSNP = simSetup_nSNP,
                         mu = 0,
                         ve = NULL,
                         he = simSetup_he,
                         specName = "bayesOpt",
                         popName = "initialPop",
                         output = "src/setupParams_FUN.rds",
                         seed = simSetup_setupSeed,
                         verbose = TRUE)
  setupFP <- digest::digest(setupSim)
  cat(paste("Setup Simulation MD5 fingerprint:", setupFP))





  # 2 Simulation parameters ------------------------------------------------------
  source("src/selectMatFun.R") # selection and mating functions

  ## 2.1 fixed parameters: ----
  fp <- list() # list of fixed parameters (constraints)

  fp$nGen = nGen                                              # number of generations
  fp$budget = totalBudget                                     # total budget
  fp$plotCost = plotCost                                      # cost for phenotyping one plot
  fp$newIndCost = newIndCost                                  # cost for creating one new individual
  fp$initPop = setupSim$initPop                                        # initial population
  fp$initPopID = paste0(setupSim$initPop$name, "_", digest(setupSim$initPop))
  fp$trait = setupSim$trait                                            # trait of interest
  fp$traitId = paste0("trait", "_", digest(setupSim$trait))
  fp$phenotyper = setupSim$phenotyper
  fp$phenotyperID = paste0("phenoLab", "_", digest(setupSim$phenotyper))
  fp$selectInds = selectGS                                    # selection function
  fp$selectIndsId = paste0("selectGS", "_", digest(selectGS))
  fp$createModel = createModel                                # function for predicting BV from phenotypes
  fp$createModelId = paste0("createModel", "_", digest(createModel))
  fp$matInds = matInds                                        # mating function
  fp$matIndsId = paste0("matInds", "_", digest(matInds))
  fp$aggrFun = mean                                          # function aggregating the calculated BV
  fp$aggrFunName = "mean"                                          # name of the function aggregating the calculated BV
  rm(setupSim)

  ## 2.2 optimized parameters: ----
  optParams <- makeParamSet(
    makeNumericParam(id = "i",
                     lower = 0.01,
                     upper = 0.99,
                     default = 0.5),
    makeNumericParam(id = "iHomo",
                     lower = 0.01,
                     upper = 0.99,
                     default = 0.5),
    makeNumericParam(id = "bRep",
                     lower = 0.01,
                     upper = 0.99,
                     default = 0.5),
    makeIntegerParam(id = "phenoFreq",
                     lower = 1,
                     upper = fp$nGen,
                     default = 1)
  )





  # Create the objective function ---------------------------------------------------------
  source("src/simulation.R")   # simulation functions
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

  if (simSetup_funEval >= 2) {
    funPar <- function(x){
      rep.x = replicate(simSetup_funEval, x, simplify = FALSE)
      res <- unlist(mclapply(rep.x, fun, mc.cores = nCpusFunEval))
      mean(res)
    }
  } else {
    funPar <- fun
  }

  obj.fun <- makeSingleObjectiveFunction(
    name = "breedingSimulation",
    id = "breedingSimulation",
    description = "breedingSimulation",
    fn = funPar,
    vectorized = FALSE,
    par.set = optParams,
    noisy = TRUE,
    minimize = FALSE,
  )

  out <- list(simSetup = list(dataFile = simSetup_dataFile,
                              nQTN = simSetup_nQTN,
                              nSNP = simSetup_nSNP,
                              he = simSetup_he,
                              setupSeed = simSetup_setupSeed,
                              funEval = funEval,
                              dataID = dataID),
              fp = fp,
              obj.fun = obj.fun,
              fun = fun)
  out
}

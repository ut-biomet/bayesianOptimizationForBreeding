# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of the simulation functions



#' Breeding simulation function parametrized with the optimized parameters
#'
#' @param i
#' @param iHomo
#' @param bRep
#' @param phenoFreq
#' @param budget
#' @param nGen
#' @param initPop
#' @param plotCost
#' @param newIndCost
#' @param trait
#' @param phenotyper
#' @param createModel
#' @param selectMateInds
#' @param aggrFun
#' @param aggrFunName
#' @param verbose
#'
#' @return aggregated results of the breeding simulation
breedSimOpt <- function(i,
                        iHomo,
                        bRep,
                        phenoFreq,
                        seed = NA,
                        budget,
                        nGen,
                        initPop,
                        plotCost,
                        newIndCost,
                        trait,
                        phenotyper,
                        createModel,
                        selectMateInds,
                        aggrFun,
                        aggrFunName,
                        verbose) {
    if (!is.na(seed)) {
      set.seed(seed)
    }

  if (identical(Sys.getenv("BAYESOPT_TEST"), 'TRUE')) {
    results <- runif(1, -10000, -9999)
    names(results) <- paste0('BV_', aggrFunName)
    return(results)
  }

  # get simulation parameters
    newParams <- list(i = i,
                      iHomo = iHomo,
                      bRep = bRep,
                      phenoFreq = phenoFreq,
                      budget = budget,
                      nGen = nGen,
                      nIndIni = initPop$nInd,
                      plotCost = plotCost,
                      newIndCost = newIndCost)
    params <- do.call(getSimulParams, newParams)

    finalGVs <- simuleBreeding(nPheno = params$nPheno,
                               nSelected = params$nSelected,
                               nNew = params$nNew,
                               nGen = nGen,
                               initPop = initPop,
                               trait = trait,
                               phenotyper = phenotyper,
                               createModel = createModel,
                               selectMateInds = selectMateInds,
                               verbose = FALSE)

    results <- aggrFun(finalGVs)
    names(results) <- paste0('BV_', aggrFunName)
    results
}




#' Breeding simulation
#'
#' @param nPheno number of phenotyped individulas for each generations
#' @param nSelected number of selected individuals for each generations
#' @param nNew number of new individuals for each generations
#' @param nGen number of generations
#' @param initPop initial population
#' @param trait trait of interest
#' @param phenotyper breedSimulatR's phenotyper object for the trait
#' @param createModel function creating a GS model for predicting genotypics
#' values
#' @param selectInds selection function
#' @param matInds mating function
#'
#' @return GVs of the final population
simuleBreeding <- function(nPheno,
                           nSelected,
                           nNew,
                           nGen,
                           initPop,
                           trait,
                           phenotyper,
                           createModel,
                           selectMateInds,
                           verbose = TRUE) {

  # initial checks ----
  if (length(nPheno) == 1) {
    nPheno <- rep(nPheno, nGen)
  } else if (length(nPheno) != nGen) {
    stop("length(nPheno) must be '1' or 'nGen'")
  }

  if (length(nSelected) == 1) {
    nSelected <- rep(nSelected, nGen)
  } else if (length(nSelected) != nGen) {
    stop("length(nSelected) must be '1' or 'nGen'")
  }

  if (length(nNew) == 1) {
    nNew <- rep(nNew, nGen)
  } else if (length(nNew) != nGen) {
    stop("length(nNew) must be '1' or 'nGen'")
  }

  # simulation initialization ----
  currentPop <- initPop
  phenoDta <- data.frame()
  genoData <- matrix(nrow = 0, ncol = ncol(currentPop$genoMat))

  # start simulation algorithm ----
  for (gen in seq(nGen)) {
    if (verbose) cat("Generation:", gen, "/", nGen, "\n")

    ## phenotyping ----
    if (nPheno[gen] != 0) {
      rep <- nPheno[gen]/currentPop$nInd
      if (floor(rep) != rep) {
        if (verbose) {
          message(paste("rep is not a multiple of nInds.",
                      "Some individual will be phenotyped 1 plot more.",
                      "(random selection)"))
        }
        nRep <- floor(rep)
        rep <- base::rep(nRep, currentPop$nInd)
        rep[sample(currentPop$nInd, nPheno[gen] - sum(rep))] <- nRep + 1
      }
      if (verbose) cat("\tModel calculation\n")
      newPhenoDta <- phenotyper$trial(currentPop, rep = rep)$data
      phenoDta <- rbind(phenoDta, newPhenoDta)
      newGenoData <- currentPop$genoMat
      genoData <- rbind(genoData, newGenoData)
      model <- createModel(phenoDta, genoData, trait$name)
    }

    # if (verbose) cat("\tModel accuracy:")
    # truthVSpred <- data.frame(GV = trait$gv(currentPop),
    #                           GP = predictGV(currentPop, model))
    # if (verbose) cat(paste("\tr2: ", cor(truthVSpred)[1,2]^2, "\n"))

    ## selection ----
    if (verbose) cat("\tindividuals selection and mating\n")
    crosses <- selectMateInds(pop = currentPop,
                              model = model,
                              nSel = nSelected[gen],
                              nNewTot = nNew[gen],
                              basename = paste0("gen", gen))

    if (verbose) cat("\tindividuals crosses\n")
    newInd <- makeCrosses(crosses, currentPop)

    currentPop <- breedSimulatR::population$new(name = paste0("gen_", gen),
                                             inds = newInd,
                                             verbose = FALSE)
  }

  trait$gv(currentPop)
}






#' Calculate the simulation parameters from the optimizations parameters
#'
#' @param i
#' @param iHomo
#' @param bRep
#' @param phenoFreq
#' @param budget
#' @param nGen
#' @param nIndIni
#' @param plotCost
#' @param newIndCost
#'
#' @return list of four elements: `nPheno`, `nSelected`, `nNew` (see documentation of `simuleBreeding` function), `eff.i`, `eff.budget` the effective budgets used for the breeding campaign.
getSimulParams <- function(i,
                           iHomo,
                           bRep,
                           phenoFreq,
                           budget,
                           nGen,
                           nIndIni,
                           plotCost = 1,
                           newIndCost = 1.07){

  # initialisation
  nSelected <- rep(0, nGen)
  nNew <- rep(0, nGen)

  # 2 - calculation of the number of new individuals for the 1st generation
  # the first generation is homozygot, so 1 individual per cross is enough
  # the number of cross is equal to the number of selected individuals
  # using the TSP mating
  nSelected[1] <- max(round(nIndIni * iHomo), 1)
  nNew[1] <- ifelse(nSelected[1] == 2, 1, nSelected[1])




  # 1 - calculation of the total number of new individuals and phenotyping plots
  x <- calcNTotnPheno(bRep = bRep,
                      budget = budget,
                      plotCost = plotCost,
                      newIndCost = newIndCost)

  indTot <- x["indTot"]
  nPhenoTot <- x["nPheno"]

  # 2 - calc number of new individuals for each generation
  nNew <- calcNnew(r = 0,
                   nTot = indTot,
                   nGen = nGen,
                   gen1 = nNew[1])

  # 3 - calc number of selected individuals for each generation
  nSel <- pmax(round(nNew[1:nGen-1] * i), 1)
  nSelected[2:nGen] <- nSel

  # 3 - calc number of phenotyping per generation
  nPheno <- calcNpheno(phenoFreq = phenoFreq,
                       nPhenoTot = nPhenoTot,
                       nGen = nGen,)

  # checks
  nInd <- c(nIndIni, nNew[1:nGen - 1])

  if (any(nSelected > nInd)) {
    warning("nSelected > nInd detected. set nSelected to nInd")
    nSelected[nSelected > nInd] <- nInd[nSelected > nInd]
  }
  if (any(nSelected > nNew)) {
    warning("nSelected > nNew detected. set nSelected to nNew")
    nSelected[nSelected > nNew] <- nNew[nSelected > nNew]
  }
  if (any(nPheno[nPheno != 0] < nInd[nPheno != 0])) {
    warning(paste("nPheno < nInd detected.",
                  "on some generation, not all individuals can be phenotyped"))
  }



  list(nPheno = nPheno,
       nSelected = nSelected,
       nNew = nNew,
       eff.i = nSelected / nInd,
       eff.budget = sum(nNew) * newIndCost + sum(nPheno) * plotCost)
}



#' calculate the number of new individuals to generate for each generations from the optimizations parameters
#'
#' @param nTot total number of new individuals to create for the entire breeding campaign.
#' @param nGen total number of generations for the breeding campaign.
#' @param gen1 (optional) number of new individuals for the first generation.
#' @param r shape parameters. If `r = 0` all generations have the same number of new individuals, if `r \in [-1;0]` the number of new individuals per generations will linearly deacrease, and if `r \in [0;+1]`, the number of new individuals per generations will linearly increase.
#'
#' @return vector of length nGen of the number of new individuals to create at each generations.
calcNnew <- function(nTot, nGen, gen1 = NA, r = 0) {

  #checks:


  if (is.na(gen1)) {
    g <- seq(nGen)
    if (nTot < nGen) {
      warning("nTot too low in front of nGen\n nTot is set to nGen")
      nTot <- nGen
    }
  } else{
    nGen <- nGen - 1
    nTot <- nTot - gen1
    if (nTot < nGen) {
      warning("nTot too low in front of gen1.\n nTot is set to gen1 + nGen")
      nTot <- nGen
    }
    g <- seq(nGen)
  }

  # raw calculation:
  if (nGen != 1) {
    if (r != 0) {
      newInds <- (nTot / nGen) * (1 - r + ((2 * r * (g - 1)) / (nGen - 1)))
      newInds <- round(newInds)
    } else {
      newInds <- floor(rep((nTot / nGen), nGen))
      gen <- sample(nGen,  nTot %% nGen)
      newInds[gen] <- newInds[gen] + 1
      newInds
    }
  } else {
    newInds <- nTot
  }


  # minimum 1 inds:
  if (any(newInds == 0) & sum(newInds == 0) == 1) {
    newInds[newInds == 0] <- 1
    newInds[newInds == max(newInds)] <- newInds[newInds == max(newInds)] - 1
  } else if (any(newInds == 0) & sum(newInds == 0) != 1) {

    newInds <- pmax(newInds, 1)
    warning(paste("impossible to calculate nNew correctly,",
                  "add new inds to have at leaste 1 new ind.\n",
                  "sum(newInds) =", sum(newInds)))
  }

  if (!is.na(gen1)) {
    newInds <- c(gen1, newInds)
  }
  newInds
}





#' Calculate the total number of individulas and phenotyping that can be done durring the breeding campaign from the total budget and the budget repartition.
#'
#' @param bRep budget repartition.
#' @param budget total budget of the breeding campaign.
#' @param plotCost cost for phenotyping one plot.
#' @param newIndCost cost for creating one new individula.
#'
#' @return named vector containing: `indTot` and `nPheno` the total number of individulto create and the total number of phenotyping that can be done durring the breeding campaign and `remain` the remaining budget.
calcNTotnPheno <- function(bRep, budget, plotCost = 1, newIndCost = 1.07){
  indB <- budget * bRep
  indTot <- max(round(indB / newIndCost), 1)

  phenoB <- budget - indTot * newIndCost
  nPheno <- max(floor(phenoB / plotCost), 1)

  remain <- budget - nPheno*plotCost - indTot*newIndCost

  if (remain < 0) {
    warning(paste("failed to calculate nTot and nPheno correctly.\n",
                  "budget unused:", budget - remain))
  }

  c(indTot = indTot,
    nPheno = nPheno,
    remain = remain)
}






#' Calculate the number of phenotyping to do at each generations from the total number of phenotyping and the "phenotyping frequency"
calcNpheno <- function(phenoFreq,
                       nPhenoTot,
                       nGen,
                       nInds = NA){

  phenoGen <- seq(1, nGen, phenoFreq)
  nPheno <- rep(0, nGen)
  nPheno[phenoGen] <- floor(nPhenoTot / length(phenoGen))

  if (sum(nPheno) != nPhenoTot) {
    warning(paste("impossible to give the same number of nPheno for each campain.",
                  "random selection of some campain with more nPheno.\n"))
    remain <- nPhenoTot - sum(nPheno)
    phenoGen <- sample(phenoGen, remain, replace = F)
    nPheno[phenoGen] <- nPheno[phenoGen] + 1
  }

  if (nPheno[1] < 3) {
    warning(paste("nPheno[1] < 3 detected, trouble with cv.glmnet.",
                  "set nPheno[1] to 3.\nadd",
                  3 - nPheno[1], "pheno."))
    nPheno[1] <-  3
  }

  nPheno
}

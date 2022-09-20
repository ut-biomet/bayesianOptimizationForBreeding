# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Definition of function setting up the simulation

#' Setup the simulation for the optimization
#'
#' @param dataFile genetic data of the initial population. ("vcf" file)
#' @param lchrCm length of the chromosome in centi-morgan
#' @param nQTN total number of QTN
#' @param nGen Number of generation
#' @param plotCost cost for phenotyping one plot
#' @param newIndCost cost for creating one new individual
#' @param selectMateInds selection and mating function
#' @param createModel function for creating prediction model
#' @param breedSimOpt Breeding simulation function parametrized with optimized
#' and fixed parameters
#' @param aggrFun function aggregating the calculated BV
#' @param aggrFunName
#' @param totalBudget total budget of the breeding campaign
#' @param plotBudgetPerGen budget express as number of plotCost per generation
#' the total budget will be calculated by:
#' `totalBudget = nGen * plotCost * plotBudgetPerGen`
#' @param nSNP total number of SNP
#' @param mu phenotypic mean (see: breedSimulatR::phenotyper)
#' @param he heritability (see: breedSimulatR::phenotyper)
#' @param ve environmental variance (see: breedSimulatR::phenotyper)
#' @param specName specie's name
#' @param popName population name
#' @param traitName trait name
#' @param outputFolder folder where to save the setup informations. Set it to
#' NULL to not save the results. Default "simSetups".
#' @param replace If they already exists, should the setup files be replaced or
#' saved with a new name. (default: TRUE)
#' @param setupName Name of this setup. It will be recored in the setup, and
#' will be used to create the file's name of the saved file. Default: "noName"
#' @param seed [optional] random seed
#' @param verbose [optional] Display information
#'
#' @return list of setup information: specie,snps,initPop = initPop,trait and
#' phenotyper. (see: breedSimulatR's package documentation)
setupSimulation <- function(dataFile,
                            lchrCm,
                            nQTN,
                            nGen,
                            plotCost,
                            newIndCost,
                            selectMateInds,
                            createModel,
                            breedSimOpt,
                            aggrFun,
                            aggrFunName,
                            totalBudget = NULL,
                            plotBudjetPerGen = NULL,
                            nSNP = NULL,
                            mu = 0,
                            he = NULL,
                            ve = NULL,
                            specName = "bayesOpt",
                            popName = "initialPop",
                            traitName = "trait1",
                            outputFolder = 'simSetups',
                            replace = TRUE,
                            setupName = "noName",
                            varSetupParams = NULL,
                            seed = NULL,
                            verbose = TRUE) {


  # test totalBudget and plotBudjetPerGen
  if ((is.null(totalBudget) && is.null(plotBudjetPerGen))
    || (!is.null(totalBudget) && !is.null(plotBudjetPerGen))) {
    stop('`totalBudget` or (exlusif) `plotBudjetPerGen` should be specified')
  }


  # Load geno data ----

  # set random seed:
  if (!is.null(seed)) {
    prevSeed <- .Random.seed
    set.seed(seed)
  }

  genoDta <- loadData(dataFile = dataFile,
                      nSNP = nSNP,
                      verbose = verbose)

  dataFileId <- paste0(basename(dataFile), "_",
                       digest(dataFile, file = TRUE))
  filteredDataId <- digest(genoDta)

  # reset the seed
  if (!is.null(seed)) {
    set.seed(prevSeed)
  }


  # Create breedSimulatR's objects ----
  BSR_Obj <- createBreedSimObj(genoDta = genoDta,
                               specName = specName,
                               lchrCm = lchrCm,
                               popName = popName,
                               nQTN = nQTN,
                               traitName = traitName,
                               mu = mu,
                               ve = ve,
                               he = he,
                               verbose = verbose)
  remove(genoDta) # clear memory

  # Create fixed parameters list ----

  # calculate total budget
  if (is.null(totalBudget)) {
    totalBudget <- nGen * plotCost * plotBudjetPerGen
  }

  fp <- list() # list of fixed parameters (constraints)
  fp$nGen           = nGen
  fp$budget         = totalBudget
  fp$plotBudjetPerGen = totalBudget / nGen
  fp$plotCost       = plotCost
  fp$newIndCost     = newIndCost
  fp$initPop        = BSR_Obj$initPop      # initial population
  fp$trait          = BSR_Obj$trait        # trait of interest
  fp$phenotyper     = BSR_Obj$phenotyper   # phenotyper
  fp$selectMateInds = selectMateInds
  fp$createModel    = createModel
  fp$aggrFun        = aggrFun
  fp$aggrFunName    = aggrFunName

  initPopId        = paste0("initPop", "_", digest(BSR_Obj$initPop))
  traitId          = paste0("trait", "_", digest(BSR_Obj$trait))
  phenotyperId     = paste0("phenoLab", "_", digest(BSR_Obj$phenotyper))
  createModelId    = paste0("createModel", "_", digest(createModel))
  selectMateIndsId = paste0("selectMateInds", "_", digest(selectMateInds))
  aggrFunNameId    = paste0(aggrFunName, "_", digest(aggrFun))

  # Create optimized parameters: ----
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



  # Create Objective function ----
  objFun <- createObjFun(fixedParams = fp,
                         breedSimOpt = breedSimOpt)


  # Finalization ----
  setup <- list(
    setupInfo = list(dataFile = dataFile,
                     dataFileId = dataFileId,
                     nSNP = nSNP,
                     setupSeed = seed,
                     filteredDataId = filteredDataId,
                     nQTN = nQTN,
                     mu = mu,
                     ve = ve,
                     he = he,
                     initPopId = initPopId,
                     traitId = traitId,
                     phenotyperId = phenotyperId,
                     createModelId = createModelId,
                     selectMateIndsId = selectMateIndsId,
                     aggrFunNameId = aggrFunNameId),
    setupName = setupName,
    fixedParams = fp,
    objFun = objFun,
    optParams = optParams,
    varOfInterest = varSetupParams
  )

  if (verbose) {
    setupHash <- digest::digest(setup)
    cat(paste("Simulation setup MD5 fingerprint:", setupHash, "\n"))
  }


  # save setup in RDS file ----
  if (!is.null(outputFolder)) {
    resultFile <- paste0(outputFolder,
                         "/simSetup_",
                         setupName)
    resF <-  paste0(resultFile, '.rds')
    if (!replace) {
      i <- 1
      while (file.exists(resF)) {
        i <- i + 1
        resF <- paste0(resultFile, '_', i, '.rds')
      }
    }
    saveRDS(setup, file = resF)
    if (verbose) {
      fileHash <- digest::digest(resF, file = TRUE)
      cat(paste("Setup file MD5 fingerprint:", fileHash, "\n"))
    }
  }

  setup
}







loadData <- function(dataFile , nSNP, verbose) {

  # read data file ----
  if (verbose) cat("Load genotypic data ...\n")

  dta <- read.vcf(dataFile, convert.chr = FALSE)

  if (verbose) {
    dtaHash <- digest::digest(dta)
    cat(paste("Data MD5 fingerprint:", dtaHash, "\n"))
  }

  # sample some SNP ----
  if (verbose) cat("\nFilter genotypic data\n")

  if (!is.null(nSNP)) {
    if ((nSNP <= 0) | (nSNP > nrow(dta@snps))) {
      stop("nSNP should be a positive number lower or equal to the total number of SNP in the provided data.")
    }
    sampledMarkers <- sort(sample(nrow(dta@snps), nSNP))
  } else {
    sampledMarkers <- seq(nrow(dta@snps))
  }

  if (verbose) {
    markersHash <- digest::digest(sampledMarkers)
    cat(paste("Markers MD5 fingerprint:", markersHash, "\n"))
  }

  # filter geno data according to the sampled markers ----
  dta <- gaston::select.snps(dta, sampledMarkers)

  dta
}








createBreedSimObj <- function(genoDta,
                           specName,
                           lchrCm,
                           popName,
                           nQTN,
                           traitName,
                           mu,
                           ve,
                           he,
                           verbose) {


  ## specie definition ---------------------
  if (verbose) cat("\nCreate breedSimulatR's specie object\n")
  chrLength <- round(tapply(genoDta@snps$pos, genoDta@snps$chr, max) * 1.05)

  simulSpecie <- breedSimulatR::specie$new(specName = specName,
                                           nChr = length(chrLength),
                                           lchr = chrLength,
                                           lchrCm = lchrCm)

  ## SNP definition ---------------------
  if (verbose) cat("\nCreate breedSimulatR's snps object\n")
  # create SNP ids
  snpId <- sprintf(fmt = paste0(
    "snp%0",
    floor(log10(nrow(genoDta@snps))) + 1, "i"),
    seq(nrow(genoDta@snps)))

  snps <- breedSimulatR::SNPinfo$new(
    SNPcoord = data.frame(
      chr = genoDta@snps$chr,
      physPos = genoDta@snps$pos,
      SNPid = snpId),
    specie = simulSpecie)


  ## population definition ----------------
  if (verbose) cat("\nCreate breedSimulatR's population object\n")
  genoDta <- gaston::as.matrix(genoDta)
  colnames(genoDta) <- snpId
  initPop <- breedSimulatR::createPop(geno = genoDta,
                                      SNPinfo = snps,
                                      popName = popName,
                                      verbose = verbose)
  if (verbose) {
    print(initPop)
  }


  ## trait definition ---------------------
  if (verbose) cat("\nCreate breedSimulatR's trait object\n")
  if (is.null(nQTN)) {
    nQTN <- snps$nSNP()
  } else if((nQTN <= 0) | (nQTN > snps$nSNP())) {
    stop("nSNP should be a positive number lower or equal to nSNP or to the total number of SNP in the provided data.")
  }

  ## sample qtn
  qtnId <- sample(snps$SNPcoord$SNPid, nQTN)

  ## qtn effects
  qtnEffects <- rexp(nQTN, 1) * sample(c(-1, 1), nQTN, replace = T)

  trait1 <- breedSimulatR::trait$new(name = traitName,
                                     qtn = qtnId,
                                     qtnEff = qtnEffects)

  ## phenotyper definition ----------------
  if (verbose) cat("\nCreate breedSimulatR's phenotyper object\n")
  phenoLab <- breedSimulatR::phenotyper$new(name = "PhenoLab1",
                                            traits = list(trait1),
                                            plotCost = 1,
                                            mu = mu,
                                            ve = ve,
                                            he = he,
                                            pop = initPop)

  ## output ----
 list(initPop = initPop,
      trait = trait1,
      phenotyper = phenoLab)
}








createObjFun <- function(fixedParams,
                         breedSimOpt) {

  # Create a simple function ----
  # This function is doing the simulation with only one vector (x)
  # as parameter. This is necessary for using `mlrMBO`:
  fun <- function(x) {
    breedSimOpt(i = x[1],
                iHomo = x[2],
                bRep = x[3],
                phenoFreq = x[4],
                seed = x[5],
                budget = fixedParams$budget,
                nGen = fixedParams$nGen,
                initPop = fixedParams$initPop,
                plotCost = fixedParams$plotCost,
                newIndCost = fixedParams$newIndCost,
                trait = fixedParams$trait,
                phenotyper = fixedParams$phenotyper,
                createModel = fixedParams$createModel,
                selectMateInds = fixedParams$selectMateInds,
                aggrFun = fixedParams$aggrFun,
                aggrFunName = fixedParams$aggrFunName,
                verbose = FALSE)
  }
  fun
}

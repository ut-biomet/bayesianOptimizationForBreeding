# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Main R script of bayesian optimization for breeding.




#' Bayesian optimization
#'
#' @param simSetup list return by the function `setupsimulation`
#' @param optP list containing:
#'  - kernel see `help(DiceKriging::km)`: a character string specifying the
#' covariance structure to be used, to be chosen between ‘"gauss"’,
#' ‘"matern5_2"’, ‘"matern3_2"’, ‘"exp"’ or ‘"powexp"’.
#'  - propose.points number of points evaluated at each iteration
#'  - nCpusProposePoints number of CPU to use for evalating points at each
#' iteration
#'  - final.method see `help(mlrMBO::makeMBOControl)` How should the final
#' point be proposed. Possible values are:
#' “best.true.y”: Return best point ever visited according to true
#' value of target function. Can be bad if target function is
#' noisy. “last.proposed”: Return the last point proposed by the
#' model. “best.predicted”: Use the final model to predict all
#' points ever visited and use the best one. This might average-out
#' noisy function values.
#'  - save.on.disk.at.time (optional) see `help(mlrMBO::makeMBOControl)` time (in seconds)
#' which have to be passed from last save to save the state of the optimization.
#' Default: Inf (never save) Note: the results will always be saved.
#'  - acquFunction see `help(mlrMBO::setMBOControlInfill)` [‘MBOInfillCrit’]
#' acquisition function.
#'  - filter.proposed.points `help(mlrMBO::setMBOControlInfill)` boolean
#' If ‘TRUE’, proposed points whose distance to design points or other
#' current candidate points is smaller than ‘filter.proposed.points.tol’,
#' are replaced by random points.
#'  - filter.proposed.points.tol `help(mlrMBO::setMBOControlInfill)`
#' Tolerance value filtering of proposed points. We currently use a
#' maximum metric to calculate the distance between points. Default is 0.0001.
#'  - totalIter total number of iterations (including the evalualtion of the
#' initial training set)
#'  - time.budget see `help(mlrMBO::setMBOControlTermination)`
#' Running time budget in seconds. Note that the actual mbo run can take more
#' time since the condition is checked after each iteration.
#'  - initTrainSize number of points to evaluate for the initial training set
#'  - funEval (optional, if null, set to 1) number of objective function
#' evaluations for each proposed points in order to reduce the noise
#'  - nCpusFunEval  (must be specify if funEval > 1) number of CPU to use to run
#' the objective function `funEval` time
#' @param outputfolder directory where to save the results
#' @param mainseed (optional) seed for random number generation
#'
#' @return
bayesianOptimization <- function(
  simSetup,
  optP,
  mainSeed = NULL,
  outputFolder = 'results'){


  startTime <- Sys.time()
  cat("Start Bayesian Optimization...\n")
  if (!is.null(mainSeed)) {
    set.seed(mainSeed)
  }

  # check funEval
  if (is.null(optP$funEval)) {
    optP$funEval <- 1
    optP$nCpusFunEval <- 1
  }



  ## save some information ----
  optP$nCpus <- optP$nCpusFunEval * optP$nCpusProposePoints
  optP$opMethod <- "Bayesian Optimization"
  resultFile <-  paste0(outputFolder,
                        "/boRes_",
                        simSetup$setupName, '_',
                        optP$setupName, '_',
                        mainSeed)
  optP$acquFunctionName <- paste0(optP$acquFunction$name, "_",
                                  digest(optP$acquFunction))
  save.file.path <-  paste0(outputFolder,"/mboProgress_",Sys.info()["nodename"],"_",
                            format(startTime, "%Y-%m-%d_%H-%M-%S"),
                            ".Rdata")
  if (optP$nCpus > detectCores()) {
    stop("nCpusProposePoints and nCpusFunEval doesn't match your total number of cores on your computer")
  }



  # mlrMBO setup ----

  ## Learner ----
  learner = makeLearner(
    "regr.km",
    predict.type = "se",
    covtype = optP$kernel,
    nugget.estim = TRUE,
    control = list(trace = FALSE),
    scaling = TRUE
  )

  ## General Control ----
  if (is.null(optP$save.on.disk.at.time)) {
    optP$save.on.disk.at.time <- Inf
  }
  control = makeMBOControl(
    n.objectives = 1,
    propose.points = optP$propose.points,
    final.method = optP$final.method,
    final.evals = 0,
    y.name = paste0("BV_", simSetup$fixedParams$aggrFunName),
    on.surrogate.error = "warn",
    save.on.disk.at.time = optP$save.on.disk.at.time,
    save.file.path = save.file.path
  )

  ## Aquisition function ----
  control = setMBOControlInfill(
    control,
    crit = optP$acquFunction,
    filter.proposed.points = optP$filter.proposed.points,
    filter.proposed.points.tol = optP$filter.proposed.points.tol
  )

  ## Multi-points proposition ----
  control = setMBOControlMultiPoint(
    control,
    method = "cl",
    cl.lie = min
  )

  ## Termination condition ----
  control = setMBOControlTermination(
    control,
    iters = optP$totalIter,
    time.budget = optP$time.budget
  )



  # Create function which repeats the simulation and return the mean values
  if (optP$funEval != 1) {
    fun <- function(x) {
      rep.x = replicate(optP$funEval, x, simplify = FALSE)
      res <- unlist(mclapply(rep.x, simSetup$objFun,
                             mc.cores = optP$nCpusFunEval))
      mean(res)
    }
  } else {
    fun <- simSetup$objFun
  }


  fun <- makeSingleObjectiveFunction(
    name = "breedingSimulation",
    id = "breedingSimulation",
    description = "breedingSimulation",
    fn = fun,
    vectorized = FALSE,
    par.set = simSetup$optParams,
    noisy = TRUE,
    minimize = FALSE)


  # Initial design ----
  initDes = generateDesign(n = optP$initTrainSize,
                           par.set = simSetup$optParams,
                           fun = lhs::improvedLHS)




  # Bayesian optimization  ----
  parallelStartMulticore(cpus = optP$nCpusProposePoints,
                         show.info = TRUE,
                         load.balancing	= TRUE)
  run <- mbo(
    fun = fun,
    design = initDes,
    learner = learner,
    control = control,
    show.info = TRUE
  )
  parallelStop()

  # Generated outputs ----
  results <- as.data.frame(run$opt.path)
  results$seed <- NA
  # add optimization parameters to results
  for (i in names(optP)) {
    if (is.numeric(optP[[i]]) | is.character(optP[[i]])) {
      results <- cbind(optP[[i]], results)
      colnames(results)[1] <- i
    }
  }


  # add fixed parameters to results
  for (i in names(simSetup$fixedParams)) {
    if (is.numeric(simSetup$fixedParams[[i]]) | is.character(simSetup$fixedParams[[i]])) {
      results <- cbind(simSetup$fixedParams[[i]], results)
      colnames(results)[1] <- i
    }
  }


  # add simulation setup parameters to results
  for (i in names(simSetup$setupInfo)) {
    if (is.numeric(simSetup$setupInfo[[i]])
        | is.character(simSetup$setupInfo[[i]])) {
      results <- cbind(simSetup$setupInfo[[i]], results)
      colnames(results)[1] <- i
    }
  }
  # final point prediction

  # if (length(run$models) != 0) {
  #   model <- run$models[[1]]$learner.model
  # } else if (length(run$final.opt.state$models) != 0) {
    model <- run$final.opt.state$models$models[[1]]$learner.model
  # } else {
  #   model <- NULL
  #   warning("Couldn't retreive the latest Model.")
  # }

  if (!is.null(model)) {
    bestPoint_pred <- predict(
      model,
      results[run$best.ind, names(simSetup$optParams$pars)],
      type = "UK"
    )$mean
  } else {
    bestPoint_pred <- NULL
  }

  # create
  output <- list(resultType = "bayesOpt",
                 results = results,
                 optP = optP,
                 optSeed = mainSeed,
                 simSetup = simSetup,
                 bestPoint = run$x,
                 bestPoint_pred = bestPoint_pred,
                 mboRawRes = run,
                 startTime = startTime)

  # save results ----
  i <- 1
  resF <-  paste0(resultFile, '.rds')
  while (file.exists(resF)) {
    i <- i + 1
    resF <- paste0(resultFile, '_', i, '.rds')
  }
  saveRDS(output, file = resF)
  output
}








#' random exploration
#'
#' @param simSetup list return by the function `setupsimulation`
#' @param optP list containing:
#'  - totalIter total number of optimization iterations to do
#'  - propose.points number of points evaluated at each iteration
#'  - funEval (optional, if null, set to 1) number of objective function
#' evaluations for each proposed points in order to reduce the noise
#'  - nCpusFunEval (must be specify if funEval > 1) number of CPU to use to run
#' the objective function `funEval` time
#'  - nCpusProposePoints number of CPU to use for evalating points at each
#' iteration
#'  - initTrainSize number of points to evaluate for the initial training set
#' @param mainseed (optional) seed for random number generation
#' @param outputfolder directory where to save the results
#'
#' @return
randomExploration <- function(
  simSetup,
  optP,
  mainSeed = NULL,
  outputFolder = 'results') {

  startTime <- Sys.time()

  # for practical reason, this function can accept optimization parameters
  # (optP) used by the `bayesianOptimization` function so in order to avoid
  # introducing missleading information in the restults let's set the
  # parameters not used by this function to NA.
  usedOptP <- c(
    'totalIter',
    'time.budget',
    'propose.points',
    'funEval',
    'nCpusFunEval',
    'initTrainSize',
    'nCpusProposePoints',
    "setupName")
  for (p in names(optP)) {
    if (!p %in% usedOptP ) {
      optP[[p]] <- NA
    }
  }


  cat("Start Random Exploration...\n")
  if (!is.null(mainSeed)) {
    set.seed(mainSeed)
  }

  # check funEval
  if (is.null(optP$funEval)) {
    optP$funEval <- 1
    optP$nCpusFunEval <- 1
  }

  ## save some information ----
  optP$nCpus <- optP$nCpusFunEval * optP$nCpusProposePoints
  optP$opMethod <- "Random Optimization"
  resultFile <-  paste0(outputFolder,
                        "/randOptRes_",
                        simSetup$setupName, '_',
                        optP$setupName, '_',
                        mainSeed)

  if (optP$nCpus > detectCores()) {
    stop("nCpusProposePoints and nCpusFunEval",
         "doesn't match your total number of cores on your computer")
  }

  # Draw optimized parameters values ----
  drawParam <- function(n){
    params <- sapply(simSetup$optParams$pars, function(param) {
      if (param$type == "numeric") {
        round(runif(n , param$lower, param$upper), digits = 4)
      } else if (param$type == "integer") {
        sample(param$lower:param$upper, n, replace = TRUE)
      }
    })
    if (optP$funEval == 1) {
      # draw list of seeds to be able to reproduce the results
      seed <- as.integer(sample.int(.Machine$integer.max, n))
      params <- cbind(params, seed)
    } else {
      params <- cbind(params, NA)
    }
  }


  # Mimic the Initial training data sampling
  initTrain_params <- drawParam(optP$initTrainSize)

  # Draw the rest of the points if totalIter is known
  if (!is.null(optP$totalIter)) {
    sampledParams <- drawParam(optP$totalIter * optP$propose.points)
  } else {
    sampledParams <- data.frame()
  }

  # Create function which repeats the simulation and return the mean values
  if (optP$funEval != 1) {
    fun <- function(x) {
      # remove the seed to not repeat the same simulation
      x[5] <- NA
      rep.x = replicate(optP$funEval, x, simplify = FALSE)
      res <- unlist(mclapply(rep.x, simSetup$objFun,
                             mc.cores = optP$nCpusFunEval))
      mean(res)
    }
  } else {
    fun <- simSetup$objFun
  }


  # Random search ----
  # initial training data
  parallelStartMulticore(cpus = optP$nCpusProposePoints,
                         show.info = TRUE,
                         load.balancing = TRUE)
  parallelRegisterLevels(levels = "objective")
  results <- parallelMap(fun,
                         split(initTrain_params, seq(nrow(initTrain_params))),
                         simplify = TRUE,
                         level = "custom.objective")
  results <- data.frame(initTrain_params, results)
  colnames(results)[6] <- paste0("BV_", simSetup$fixedParams$aggrFunName)
  results$dob <- 0


  if (is.null(optP$time.budget)) {
    # results <- unlist(mclapply(sampledParamsList,
    #                            simSetup$obj.fun,
    #                            mc.cores = nCpusProposePoints))
    # `parallelMap` seems to use less memory.
    res <- parallelMap(fun,
                           split(sampledParams, seq(nrow(sampledParams))),
                           simplify = TRUE,
                           level = "custom.objective")
    res <- data.frame(sampledParams, res)
    colnames(res)[6] <- paste0("BV_", simSetup$fixedParams$aggrFunName)
    res$dob <- c(rep(1:optP$totalIter, each = optP$propose.points))
    results <- rbind(results, res)
  } else {
    i <- 1
    stop <- FALSE
    while (!stop) {
      # get parameter set for the current iteration
      from <- (i - 1) * optP$propose.points + 1
      to <- i * optP$propose.points
      if (to <= nrow(sampledParams)) {
        # use pre-generated parameters
        iter_params <- sampledParams[from:to,]
      } else {
        # draw parameters
        iter_params <- drawParam(optP$propose.points)
      }

      iter_results <- parallelMap(fun,
                             split(iter_params, seq(nrow(iter_params))),
                             simplify = TRUE,
                             level = "custom.objective")
      iter_results <- data.frame(iter_params, iter_results)
      colnames(iter_results)[6] <- paste0("BV_", simSetup$fixedParams$aggrFunName)
      iter_results$dob <- i
      results <- rbind(results, iter_results)

      ellapsTime <- as.numeric(difftime(Sys.time(),
                                        startTime,
                                        units = c("secs")))

      if (ellapsTime > optP$time.budget) {
        stop <- TRUE
      }
      i <- i + 1
      if (!is.null(optP$totalIter)) {
        if (i > optP$totalIter) {
          stop <- TRUE
        }
      }
    }
  }
  parallelStop()




  # add optimization parameters to results
  for (i in names(optP)) {
    if (is.numeric(optP[[i]]) || is.character(optP[[i]])) {
      results <- cbind(optP[[i]], results)
      colnames(results)[1] <- i
    }
  }

  # add fixed parameters to results
  for (i in names(simSetup$fixedParams)) {
    if (is.numeric(simSetup$fixedParams[[i]]) | is.character(simSetup$fixedParams[[i]])) {
      results <- cbind(simSetup$fixedParams[[i]], results)
      colnames(results)[1] <- i
    }
  }

  # add simulation setup parameters to results
  for (i in names(simSetup$setupInfo)) {
    if (is.numeric(simSetup$setupInfo[[i]])
        | is.character(simSetup$setupInfo[[i]])) {
      results <- cbind(simSetup$setupInfo[[i]], results)
      colnames(results)[1] <- i
    }
  }


  # get best point
  bestPoint_val <- max(results$BV_mean)
  id <- which(results$BV_mean == bestPoint_val)
  bestPoint <- list(i = results$i[id],
                    iHomo = results$iHomo[id],
                    bRep = results$bRep[id],
                    phenoFreq = results$phenoFreq[id])


  # create
  output <- list(resultType = "randOpt",
                 results = results,
                 optP = optP,
                 optSeed = mainSeed,
                 simSetup = simSetup,
                 bestPoint = bestPoint,
                 bestPoint_val = bestPoint_val,
                 startTime = startTime)

  # save results ----
  i <- 1
  resF <- paste0(resultFile, '.rds')
  while (file.exists(resF)) {
    i <- i + 1
    resF <- paste0(resultFile, '_', i, '.rds')
  }
  saveRDS(output, file = resF)

  output
}









#' Repeat simulation using the result of optimization
#'
#' @param optRes optimization results
#' @param nRep number of simulation to do
#' @param nCpus number of core to use for running these simulation in parallel
#' @param mainseed (optional) seed for random number generation
#' @param outputfolder directory where to save the results
repeatOptResults <- function(optRes,
                             nRep,
                             nCpus = 1,
                             mainSeed = NULL,
                             outputFolder = 'results') {

  startTime <- Sys.time()

  cat("Start optimization results repetition...\n")
  ## save some information ----
  resultFile <- paste0(outputFolder,
                       "/repRes_",
                       optRes$resultType, '_',
                       optRes$simSetup$setupName, '_',
                       optRes$optP$setupName, '_',
                       mainSeed)
  if (nCpus > detectCores()) {
    stop("nCpus doesn't match your total number of cores on your computer")
  }


  ## Repeat opt results ----
  results <- singleSimulation(i = optRes$bestPoint$i,
                              iHomo = optRes$bestPoint$iHomo,
                              bRep = optRes$bestPoint$bRep,
                              phenoFreq = optRes$bestPoint$phenoFreq,
                              setup = optRes$simSetup,
                              seeds = mainSeed,
                              nRep = nRep,
                              nCpus = nCpus)

  ## Generate output ----
  output <- results[c('results', 'params')]
  output$results$opMethod <- unique(optRes$results$opMethod)

  output$resultType = 'repRes'
  output$optResults = optRes
  output$mainSeed = mainSeed
  output$startTime = startTime


  ## Save results ----
  i <- 1
  resF <- paste0(resultFile, '.rds')
  while (file.exists(resF)) {
    i <- i + 1
    resF <- paste0(resultFile, '_', i, '.rds')
  }
  saveRDS(output, file = resF)

  output
}



#' Simulate one or multiple breeding campaign using the optimized parameters.
#'
#' @param i
#' @param iHomo
#' @param bRep
#' @param phenoFreq
#' @param setup list return by the function `setupsimulation`.
#' @param seeds [optional] Either a vector of random number generation seeds
#' for each replication. (The lenght of this vector should be equal to `nRep`.)
#' Or a single seed which will be used to generate seeds for each simulation)
#' @param nRep [optional] number of replication to do.
#' @param nCpus [optional] number of cpus to use for running the simulations in
#' parallel.
#'
#' @return
#' @export
#'
#' @examples
singleSimulation <- function(i,
                             iHomo,
                             bRep,
                             phenoFreq,
                             setup,
                             seeds = NULL,
                             nRep = 1,
                             nCpus = 1){

  startTime <- Sys.time()

  # create list of parameters for parallel computing
  rep_x <- matrix(rep(c(i, iHomo, bRep, phenoFreq), nRep),
                  nrow = nRep,
                  byrow = TRUE)

  if (!is.null(seeds)) {
    if (length(seeds == 1)) {
      set.seed(as.integer(seeds))
      seeds <- as.integer(sample.int(.Machine$integer.max, nRep))
    } else if (length(seeds) == nRep){
      seeds <- as.integer(seeds)
    } else {
      stop('length of seeds should be one or nRep')
    }
  } else {
    seeds <- NA
  }
  rep_x <- base::cbind(rep_x, seeds)
  rep_x <- split(rep_x, seq(nRep))

  # run the simulation in parallel
  res <- unlist(mclapply(rep_x, setup$objFun, mc.cores = nCpus))


  # aggregate results
  results <- data.frame(i = i,
                        iHomo = iHomo,
                        bRep = bRep,
                        phenoFreq = phenoFreq,
                        seed = seeds,
                        res)
  rownames(results) <- seq(nRep)

  fp <- setup$fixedParams
  colnames(results)[6] <- paste0("BV_", fp$aggrFunName)

  # add fixed parameters to results
  for (j in names(fp)) {
    if (is.numeric(fp[[j]]) | is.character(fp[[j]])) {
      results <- cbind(fp[[j]], results)
      colnames(results)[1] <- j
    }
  }

  # add simulation setup informations to results
  for (j in names(setup$setupInfo)) {
    if (is.numeric(setup$setupInfo[[j]]) | is.character(setup$setupInfo[[j]])) {
      results <- cbind(setup$setupInfo[[j]], results)
      colnames(results)[1] <- j
    }
  }

  # create
  output <- list(results = results,
                 params = list(i = i,
                               iHomo = iHomo,
                               bRep = bRep,
                               phenoFreq = phenoFreq,
                               seeds = seeds,
                               nRep = nRep,
                               nCpus = nCpus),
                 setup = setup,
                 startTime = startTime)
  output

}

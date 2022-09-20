# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Functions used for analysing the results








#' Get the number of repeated optimization for each scenario
#'
#' @param aggResults aggregated results of optimizations repetition
#'
#' @return data.frame
get_nRep_by_scenario <- function(aggResults) {
  scenarios <- t(sapply(aggResults, function(x) {
    c(he = x$optResults$simSetup$setupInfo$he,
      nGen = x$optResults$simSetup$fixedParams$nGen,
      budget = x$optResults$simSetup$fixedParams$budget,
      optMethod = x$optResults$resultType,
      n_optResultEvaluation = x$params$nRep)
  }))
  scenarios <- data.frame(scenarios)
  scenarios$budget <- as.numeric(scenarios$budget)
  scenarios$nGen <- as.numeric(scenarios$nGen)
  scenarios$genB <- scenarios$budget / scenarios$nGen

  nRepetitions <- scenarios %>%
    group_by(he, nGen, genB, optMethod, n_optResultEvaluation) %>%
    summarise(n_repeted_opt = n())
}






#' Merge all optimization progress result in one data.frame
#'
#' @param aggResults aggregated results of optimizations
#'
#' @return list of data.frame
merge_allData <- function(aggResults) {

  # 1. merge optimization data:
  allDta_optim <- do.call(bind_rows, lapply(aggResults, function(x){
    genBudget <- (as.numeric(x$optResults$simSetup$fixedParams$budget)
                  / as.numeric(x$optResults$simSetup$fixedParams$nGen))
    scenario <- paste0(
      'he_', x$optResults$simSetup$setupInfo$he,
      '_nGen_', x$optResults$simSetup$fixedParams$nGen,
      '_b_', genBudget
    )
    res <- cbind(x$optResults$results,
                 optSeed = x$optResults$optSeed,
                 scenario = scenario)
    res$selectedOptimal <- (res$i == unique(x$results$i)
                            & res$iHomo == unique(x$results$iHomo)
                            & res$phenoFreq == unique(x$results$phenoFreq)
                            & res$bRep == unique(x$results$bRep))
    res
  }))

  allDta_optim$genB <- allDta_optim$budget / allDta_optim$nGen
  allDta_optim$global_max_dob <- max(allDta_optim$dob)
  allDta_optim$global_max_bv <- max(allDta_optim$BV_mean)
  allDta_optim$global_min_bv <- min(allDta_optim$BV_mean)


  # 2. merge optimization result evaluations
  allDta_resRep <- do.call(bind_rows, lapply(aggResults, function(x){
    genBudget <- (as.numeric(x$optResults$simSetup$fixedParams$budget)
                  / as.numeric(x$optResults$simSetup$fixedParams$nGen))
    scenario <- paste0(
      'he_', x$optResults$simSetup$setupInfo$he,
      '_nGen_', x$optResults$simSetup$fixedParams$nGen,
      '_b_', genBudget
    )
    res <- cbind(x$results,
                 optSeed = x$optResults$optSeed,
                 scenario = scenario)
  }))
  allDta_resRep$genB <- allDta_resRep$budget / allDta_resRep$nGen


  # assign a "repetition id" to each arbitrary couple of bayesian optimization
  # and random optimization for all scenarios.

  # split data by scenarios and create a hash-map of the `optSeed` to the
  # repetion id `repId `
  list_optSeedByScenario <- distinct(allDta_resRep[, c('opMethod', 'scenario', 'optSeed')]) %>% group_by(scenario) %>% group_map(mutate)
  idHashMap <- rep(NA, length(unique(allDta_resRep$optSeed)))
  names(idHashMap) <- unique(allDta_resRep$optSeed)

  for (tmpDta in list_optSeedByScenario) {
    for (meth in unique(tmpDta$opMethod)) {
      ids <- seq(nrow(tmpDta[tmpDta$opMethod == meth,]))
      idHashMap[as.character(tmpDta$optSeed[tmpDta$opMethod == meth])] <- ids
    }
  }

  allDta_optim$repId <- idHashMap[as.character(allDta_optim$optSeed)]
  allDta_resRep$repId <- idHashMap[as.character(allDta_resRep$optSeed)]

  allDta <- list(optimDta = allDta_optim,
                 resRepDta = allDta_resRep)
  return(allDta)
}






#' Split optimization results data.frame of one scenario into a list of
#' optimization results data.frame according to each optimization repetition
#'
#' The function assign an arbitrary id to each optimization and each
#' optimization method:
#'  - bayesian optimization 1
#'  - random optimization 1
#'  - bayesian optimization 2
#'  - random optimization 2
#'  - etc...
#'
#' So that we can easily compare the behaviour of the two optimizations methods
#'
#' @param dta optimization results data.frame of one scenario
#'
#' @return list of data.frame
split_by_optimization <- function(dta) {

  # assign arbitrary id to each optimization by optimization method
  dta$id <- NA
  for (optM in unique(dta$opMethod)) {
    tmp <- dta[dta$opMethod == optM,]
    ids <- seq(length(unique(tmp$optSeed)))
    names(ids) <- as.character(unique(tmp$optSeed))
    # tmp$id <- as.numeric(factor(tmp$optSeed, ordered = F))
    tmp$id <- ids[as.character(tmp$optSeed)]
    dta$id[dta$opMethod == optM] <- tmp$id
  }

  # save min / max Breeding values, max nIter for the scenario
  dta$min_scenario_bv <- min(dta$BV_mean)
  dta$max_scenario_bv <- max(dta$BV_mean)
  if ('dob' %in% colnames(dta)) {
    dta$max_scenario_dob <- max(dta$dob)
  }

  dta %>% group_by(id) %>% group_map(mutate)
}

getMostRepresentativeOptim <- function(dta) {
  cumMaxDta <- dta %>%
    group_by(opMethod, repId, dob) %>%
    summarise(maxBV = max(BV_mean)) %>% # max at each iter*repetition
    group_by(opMethod, repId) %>%
    summarise(cumMax = cummax(maxBV), dob) # cummax at each iter*repetition

  meanCumMax <- cumMaxDta %>%
    group_by(opMethod, dob) %>%
    summarise(meanCumMax = mean(cumMax))


  # cumMaxDta$repId <- as.factor(cumMaxDta$repId)
  # meanCumMax$repId <- "meanCumMax"
  # colnames(meanCumMax) <- c("opMethod", "dob", "cumMax", "repId")
  # tmpDta <-  rbind(cumMaxDta[cumMaxDta$opMethod == "Bayesian Optimization",],
  #                  meanCumMax[meanCumMax$opMethod == "Bayesian Optimization",])
  # p <- ggplot(tmpDta, aes(x = dob, y = cumMax, col = repId == 'meanCumMax', group = repId))
  # p <- p + geom_line()

  minRepId <- rep(NA, length(unique(cumMaxDta$opMethod)))
  names(minRepId) <- unique(cumMaxDta$opMethod)
  for (method in unique(cumMaxDta$opMethod)) {
    x <- cumMaxDta[cumMaxDta$opMethod == method,]
    currMin <- Inf
    for (repId in unique(x$repId)) {
      x2 <- full_join(meanCumMax[meanCumMax$opMethod == method,],
                      x[x$repId == repId,],
                      by = "dob")
      mse <- mean((x2$meanCumMax - x2$cumMax)^2)
      if (mse < currMin ) {
        currMin <- mse
        minRepId[method] <- repId
      }
    }
  }


  out <- do.call(rbind, lapply(unique(meanCumMax$opMethod), function(meth){
    dta[dta$opMethod == meth & dta$repId == minRepId[meth],]
  }))
  out$id <- paste0('bo_', minRepId['Bayesian Optimization'], '-',
                   'ro_', minRepId['Random Optimization'])

  # save min / max Breeding values, max nIter for the scenario
  out$min_scenario_bv <- min(out$BV_mean)
  out$max_scenario_bv <- max(out$BV_mean)
  if ('dob' %in% colnames(dta)) {
    out$max_scenario_dob <- max(dta$dob)
  }

  out

}

# plots ----

#' draw and save a box plot of the number of maximum optimization iteration
#' for each method.
#' The data must be
#'
#' @param dta optimization results data frame of one scenario
#' @param outDir main folder where to save the results.
#'
#' @return
boxPlot_nIter <- function(dta, outDir) {
  # get number of iteration for each optimization method and each repetition
  optLenght <- dta %>%
    group_by(optSeed, opMethod) %>%
    summarise(nIter = max(dob))
  # summarise by optimization method
  stats <- optLenght %>% group_by(opMethod) %>%
    summarise(min_nIter = min(nIter),
              `1stQu_nIter` = quantile(nIter, 0.25),
              med_nIter = median(nIter),
              mean_nIter = mean(nIter),
              `3rdQu_nIter` = quantile(nIter, 0.75),
              max_nIter = max(nIter),
              nRep = n())

  p <- (ggplot(optLenght, aes(x = opMethod, y = nIter))
        + geom_boxplot()
        + annotate("text",
                   x = stats$opMethod,
                   y = stats$med_nIter,
                   label = paste('n obs:', stats$nRep),
                   vjust = -0.5)
        + annotate("text",
                   x = stats$opMethod,
                   y = stats$min_nIter,
                   label = paste('min:', stats$min_nIter),
                   hjust = -0.1)
        + annotate("text",
                   x = stats$opMethod,
                   y = stats$max_nIter,
                   label = paste('max:', stats$max_nIter),
                   hjust = -0.1)
        + coord_cartesian(ylim = c(0, unique(dta$global_max_dob) + 1))
        + labs(x = "Optimization method",
               y = "Number of optimization iterations",
               title = "Number of optimization iterations",
               subtitle = paste0(
                 'he = ', unique(dta$he),
                 '; nGen = ', unique(dta$nGen),
                 '; B/gen = ', unique(dta$budget/dta$nGen),
                 '\nnumber of Bayesian optimization: ', sum(optLenght$opMethod == "Bayesian Optimization"),
                 '\nnumber of Random optimization: ', sum(optLenght$opMethod == "Random Optimization"))
        ))
  file <- file.path(outDir, 'nIter_boxPlot', unique(dta$scenario))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}



#' Given the data of one scenario, draw a boxplot of the cumulative maximum for
#' each optimization iteration
#'
#' @param dta optimization results data frame of one scenario
#' @param outDir main folder where to save the results.
#'
#' @return
boxPlot_cumMax <- function(dta, outDir) {
  cumMaxDta <- dta %>%
    group_by(opMethod, repId, dob) %>%
    summarise(maxBV = max(BV_mean)) %>% # max at each iter*repetition
    group_by(opMethod, repId) %>%
    summarise(cumMax = cummax(maxBV), dob) # cummax at each iter*repetition
  cumMaxDta$dob <- as.factor(cumMaxDta$dob)

  p <- (ggplot(cumMaxDta, aes(x = dob, y = cumMax, col = opMethod))
        + geom_boxplot() #
        # + geom_boxplot(data = cumMaxDta_ro) #
        # + geo_boxplot() #aes(col = opMethod)
        + coord_cartesian(xlim = c(0, unique(dta$global_max_dob)),
                          ylim = c(unique(dta$global_min_bv) - 5,
                                   unique(dta$global_max_bv) + 5))
        + labs(x = "Optimization iteration",
               y = "Mean breeding value of the final population",
               col = "Optimization method",
               title = "Cumulative observed maximum over the optimization iterations",
               subtitle = paste0(
                 'he = ', unique(dta$he),
                 '; nGen = ', unique(dta$nGen),
                 '; B/gen = ', unique(dta$budget/dta$nGen))
        )
  )
  file <- file.path(outDir, 'cumMax_boxPlot', paste0('cumMax_boxPlot_', unique(dta$scenario)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}






#' Given the data of one scenario, draw all the lines of the cumulative maximum
#' for each optimization iteration
#'
#' @param dta optimization results data frame of one scenario
#' @param outDir main folder where to save the results.
#'
#' @return
lines_cumMax <- function(dta, outDir) {
  cumMaxDta <- dta %>%
    group_by(opMethod, repId, dob) %>%
    summarise(maxBV = max(BV_mean)) %>% # max at each iter*repetition
    group_by(opMethod, repId) %>%
    summarise(cumMax = cummax(maxBV), dob) # cummax at each iter*repetition
  cumMaxDta$optId <- paste(cumMaxDta$opMethod, cumMaxDta$repId)

  p <- (ggplot(cumMaxDta, aes(x = dob, y = cumMax, col = opMethod, group = optId)))
  p <- p + geom_line(alpha = 1/10)
  p <- (p + coord_cartesian(xlim = c(0, unique(dta$global_max_dob)),
                            ylim = c(unique(dta$global_min_bv) - 5,
                                     unique(dta$global_max_bv) + 5))
        + labs(x = "Optimization iteration",
               y = "Mean breeding value of the final population",
               col = "Optimization method",
               title = "Cumulative observed maximum over the optimization iterations",
               subtitle = paste0(
                 'he = ', unique(dta$he),
                 '; nGen = ', unique(dta$nGen),
                 '; B/gen = ', unique(dta$budget/dta$nGen))
        ))

  file <- file.path(outDir, 'cumMax_lines', paste0('cumMax_lines_', unique(dta$scenario)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}









#' Draw the optimization progress for BO and RO.
#'
#' @param dta optimization result for one optimization repetition
#' @param outDir main folder where to save the results.
#'
#' @return
plot_optimization_progress <- function(dta, outDir) {
  dta$prop.type[is.na(dta$prop.type)] <- 'rand_opt'

  # calculate cumulative maximum
  dta_cumMax <- dta %>% group_by(opMethod, dob) %>%
    summarise(maxBV_mean = max(BV_mean)) %>%
    group_by(opMethod) %>%
    summarise(cumMax = cummax(maxBV_mean), dob)

  # draw plot
  p <- (ggplot(dta, aes(x = dob, y = BV_mean))
        + geom_point(aes(col = opMethod, shape = prop.type))
        + geom_line(data = dta_cumMax,
                    aes(x = dob,
                        y = cumMax,
                        col = opMethod))
        + coord_cartesian(xlim = c(0, unique(dta$max_scenario_dob)),
                          ylim = c(unique(dta$min_scenario_bv) - 5,
                                   unique(dta$max_scenario_bv) + 5))
        + annotate("segment",
                   x = dta$dob[dta$selectedOptimal],
                   xend = dta$dob[dta$selectedOptimal],
                   y = dta$BV_mean[dta$selectedOptimal] + 30,
                   yend = dta$BV_mean[dta$selectedOptimal],
                   size = 0.5,
                   arrow = arrow(length = unit(0.2, "cm")))
        + annotate("text",
                   x = dta$dob[dta$selectedOptimal],
                   y = dta$BV_mean[dta$selectedOptimal] + 30,
                   label = paste('Optimal point'),
                   hjust = 0)
        + labs(x = "Optimization iteration",
               y = "Mean breeding value of the final population",
               col = "Optimizaton method",
               shape = "Sampling method",
               title = "Optimization progress",
               subtitle = paste0(
                 'Scenario: ',
                 'he = ', unique(dta$he),
                 '; nGen = ', unique(dta$nGen),
                 '; B/gen = ', unique(dta$budget/dta$nGen)#,
                 # '; optimization repetition: ',
                 # unique(dta$id))
                 )
        )
  )

  file <- file.path(outDir, unique(dta$scenario), 'optProgress',
                    paste0('optProgress_', unique(dta$scenario), '_repetition_', unique(dta$id)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}







#' Draw PCA plots for optimization progress:
#'  - Active individuals: points explored by the BO algo
#'  - Active variables: optimizes variables: i, iHomo, bRep, phenoFreq
#'  - Suplementary individuals: points explored by the RO algo
#'
#'  - color:
#'    - (1st graph) iteration of the points
#'    - (2nd graph) value of the objective function (BV_mean)
#'
#' @param dta optimization result for one optimization repetition
#' @param outDir main folder where to save the results.
#'
#' @return
plot_optimization_pca <- function(dta, outDir) {
  dta$prop.type[is.na(dta$prop.type)] <- 'rand_opt'
  optVar <- c('i', 'iHomo', 'bRep', 'phenoFreq')
  pcaDta <- rbind(dta[dta$opMethod == "Bayesian Optimization", optVar],
                  dta[dta$opMethod == "Random Optimization", optVar]
  )

  npcaInd <- sum(dta$opMethod == "Bayesian Optimization")
  pcaRes <- PCA(pcaDta,
                scale.unit = TRUE,
                ind.sup = seq(npcaInd + 1, nrow(pcaDta)),
                ncp = 2,
                graph = T)

  pcaDta <- mutate(dta[dta$opMethod == "Bayesian Optimization",],
                   pca1 = pcaRes$ind$coord[, 1],
                   pca2 = pcaRes$ind$coord[, 2])

  pcaDta <- rbind(pcaDta,
                  mutate(dta[dta$opMethod == "Random Optimization",],
                         pca1 = pcaRes$ind.sup$coord[, 1],
                         pca2 = pcaRes$ind.sup$coord[, 2]))
  dtaList <- list()
  dtaList$BO <- pcaDta[pcaDta$opMethod == "Bayesian Optimization",]
  dtaList$RO <- pcaDta[pcaDta$opMethod == "Random Optimization",]

  plotList <- list()

  for (m in names(dtaList)) {
    method <- unique(dtaList[[m]]$opMethod)
    method <- paste0(unlist(
      regmatches(method,
                 gregexpr('(?<=^|\\s)[\\w]', method, perl = TRUE))
    ), collapse = '')
    p <- (
      ggplot(dtaList[[m]],aes(x = pca1,
                              y = pca2,
                              col = dob,
                              shape = prop.type))
      + geom_vline(xintercept = 0)
      + geom_hline(yintercept = 0)
      # add variables arrows
      + annotate("segment",
                 x = 0,
                 xend = pcaRes$var$coord[,1],
                 y = 0,
                 yend = pcaRes$var$coord[,2],
                 size = 0.5,
                 arrow = arrow(length = unit(0.2, "cm")))
    # + annotate("text",
    #            x = pcaRes$var$coord[,1],
    #            y = pcaRes$var$coord[,2],
    #            label = row.names(pcaRes$var$coord),
    #            hjust = 1)
    # add optimal point
    + geom_point(aes())
    + annotate("segment",
               x = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal] - 1,
               xend = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal],
               y = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal] + 1,
                 yend = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal],
                 size = 0.5,
                 arrow = arrow(length = unit(0.2, "cm")))
      + annotate("text",
                 x = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal] - 1,
                 y = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal] + 1,
                 label = paste('Optimal point:', method),
                 hjust = 1)
      + coord_cartesian(xlim = c(min(pcaDta$pca1) - 0.5,
                                 max(pcaDta$pca1) + 0.5),
                        ylim = c(min(pcaDta$pca2) - 0.5,
                                 max(pcaDta$pca2) + 0.5))
      + labs(x = "PCA axis 1",
             y = "PCA axis 2",
             col = "Optimization iteration",
             shape = "Sampling method",
             title = "Optimization points PCA visualisation - Iteration",
             subtitle = paste0(
               'he = ', unique(dtaList[[m]]$he),
               '; nGen = ', unique(dtaList[[m]]$nGen),
               '; B/gen = ', unique(dtaList[[m]]$budget/dtaList[[m]]$nGen)))#,
               # '; optimization repetition: ',
               # unique(dtaList[[m]]$id),
               # '; ', unique(dtaList[[m]]$opMethod)))
      + scale_color_gradientn(
        colours = rev(brewer.pal(10, "Spectral"))
      )
    )
    # same range for x and y axes
    axisLim <- c(min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1],
                     ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]),
                 max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2],
                     ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
    )
    # p <- p + xlim(axisLim) + ylim(axisLim)
    suppressMessages({
      p <- p + coord_fixed(ratio = 1,
                           xlim = axisLim,
                           ylim = axisLim)
    })



    file <- file.path(outDir, unique(dtaList[[m]]$scenario), 'pca',
                      paste0('pca_nIter_', unique(dta$scenario), '_rep_', unique(dtaList[[m]]$id),
                             '_', method))
    suppressWarnings(dir.create(dirname(file), recursive = TRUE))
    ggsave(file, p)
    plotList[[m]]$n_iter <- p

    p <- (
      ggplot(dtaList[[m]],aes(x = pca1,
                            y = pca2,
                            col = BV_mean,
                            shape = prop.type))
      + geom_vline(xintercept = 0)
      + geom_hline(yintercept = 0)
      + annotate("segment",
                 x = 0,
                 xend = pcaRes$var$coord[,1],
                 y = 0,
                 yend = pcaRes$var$coord[,2],
                 size = 0.5,
                 arrow = arrow(length = unit(0.2, "cm")))
      + geom_point(aes())
      + annotate("segment",
                 x = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal] - 1,
                 xend = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal],
                 y = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal] + 1,
                 yend = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal],
                 size = 0.5,
                 arrow = arrow(length = unit(0.2, "cm")))
      + annotate("text",
                 x = dtaList[[m]]$pca1[dtaList[[m]]$selectedOptimal] - 1,
                 y = dtaList[[m]]$pca2[dtaList[[m]]$selectedOptimal] + 1,
                 label = paste('Optimal point:', method),
                 hjust = 1)
      + coord_cartesian(xlim = c(min(pcaDta$pca1) - 0.5,
                                 max(pcaDta$pca1) + 0.5),
                        ylim = c(min(pcaDta$pca2) - 0.5,
                                 max(pcaDta$pca2) + 0.5))
      + labs(x = "PCA axis 1",
             y = "PCA axis 2",
             col = "Mean breeding value ",
             shape = "Sampling method",
             title = "Optimization points PCA visualisation - Breeding Value",
             subtitle = paste0(
               'he = ', unique(dtaList[[m]]$he),
               '; nGen = ', unique(dtaList[[m]]$nGen),
               '; B/gen = ', unique(dtaList[[m]]$budget/dtaList[[m]]$nGen)))#,
               # '; optimization repetition: ',
               # unique(dtaList[[m]]$id),
               # '; ', unique(dtaList[[m]]$opMethod)))
      + scale_color_gradientn(
        colours = rev(brewer.pal(10, "Spectral"))
      )
    )
    # same range for x and y axes
    axisLim <- c(min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1],
                     ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]),
                 max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2],
                     ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
    )
    # p <- p + xlim(axisLim) + ylim(axisLim)

    suppressMessages({
      p <- p + coord_fixed(ratio = 1,
                           xlim = axisLim,
                           ylim = axisLim)
    })
    p <- p + theme(legend.box.margin = margin(0, 0, 0, 10))



    file <- file.path(outDir, unique(dtaList[[m]]$scenario), 'pca',
                      paste0('pca_BVmean_', unique(dta$scenario), '_rep_', unique(dtaList[[m]]$id),
                             '_', method))
    suppressWarnings(dir.create(dirname(file), recursive = TRUE))
    ggsave(file, p)
    plotList[[m]]$BV_mean <- p
  }


  p <- (fviz_pca_var(pcaRes)
        + labs(
          title = 'Optimization points PCA - Graph of variables',
          subtitle = paste0(
            'he = ', unique(dtaList[[m]]$he),
            '; nGen = ', unique(dtaList[[m]]$nGen),
            '; B/gen = ', unique(dtaList[[m]]$budget/dtaList[[m]]$nGen)))
            # '; optimization repetition: ',
            # unique(dtaList[[m]]$id)))
        + theme(panel.background = element_rect('white')))

  plotList$var <- p
  file <- file.path(outDir, unique(pcaDta$scenario), 'pca',
                    paste0('pca_variables_', unique(dta$scenario), '_rep_', unique(pcaDta$id),
                           '_', method))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  plotList
}





#' Plot empirical cumulative distribution of the optimization's results evaluation
#'
#' @param dta data.frame of the optimization result evaluation for one scenario
#' @param outDir main folder where to save the results.
#'
#' @return
plot_empCumDist <- function(dta, outDir) {

  # 2 by 2 comparisons
  dta_BO <- dta[dta$opMethod == 'Bayesian Optimization',]
  dta_RO <- dta[dta$opMethod == 'Random Optimization',]

  # boHigher <- rep(NA, nrow(dta_BO) * nrow(dta_RO))
  boHigher <- c()
  for (i in seq(nrow(dta_BO))) {
    boHigher <- c(boHigher, as.numeric(dta_BO[i, 'BV_mean']) > dta_RO$BV_mean)
  }
  boHigher_prop <- sum(boHigher) / length(boHigher)

  # get mean value for each optimization result
  minEval <- min(dta$BV_mean)
  maxEval <- max(dta$BV_mean)

  meanDta <- dta %>% group_by(opMethod, optSeed) %>%
    summarise(meanRestult = mean(BV_mean))

  meanDta_BO <- meanDta[meanDta$opMethod == 'Bayesian Optimization',]
  meanDta_RO <- meanDta[meanDta$opMethod == 'Random Optimization',]

  cumDistFun_BO <- as.function(ecdf(meanDta_BO$meanRestult))
  cumDistFun_RO <- as.function(ecdf(meanDta_RO$meanRestult))

  cumDistDta <- data.frame(x = seq(minEval, maxEval, length.out = 200))
  cumDistDta$BO <- cumDistFun_BO(cumDistDta$x)
  cumDistDta$RO <- cumDistFun_RO(cumDistDta$x)

  cumDistDta <- cumDistDta %>% pivot_longer(c("BO","RO"),
                                            names_to = 'opMethod',
                                            values_to = 'quantile',
                                            values_drop_na = FALSE)
  cumDistDta$opMethod[cumDistDta$opMethod == 'BO'] <- 'Bayesian Optimization'
  cumDistDta$opMethod[cumDistDta$opMethod == 'RO'] <- 'Random Optimization'

  cumDistDta$opMethod[cumDistDta$opMethod == 'Bayesian Optimization'] <- paste('Bayesian Optimization\nnObs:',sum(meanDta$opMethod == 'Bayesian Optimization'))
  cumDistDta$opMethod[cumDistDta$opMethod == 'Random Optimization'] <- paste('Random Optimization\nnObs:',sum(meanDta$opMethod == 'Random Optimization'))
  p <- (ggplot(cumDistDta, aes(x = x, y = quantile, col = opMethod))
        + geom_line()
        + labs(x = "Mean breeding value",
               y = "Quantile",
               col = "Optimisation method",
               title = "ECDF of the optimizations results evaluations",
               subtitle = paste0(
                 'he = ', unique(dta$he),
                 '; nGen = ', unique(dta$nGen),
                 '; B/gen = ', unique(dta$budget/dta$nGen),
                 # '\nBayesian Optimization: ',
                 # sum(meanDta$opMethod == 'Bayesian Optimization'), ' points',
                 # '\nRandom Optimization: ',
                 # sum(meanDta$opMethod == 'Random Optimization'), ' points',
                 '\nBO is better ',
                 round(boHigher_prop*100, 2), '% of the time'))
  )
  file <- file.path(outDir, 'ECDF_optEval', paste0('ECDF_', unique(dta$scenario)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}


plot_boxPlot_optResEval <- function(dta, outDir) {

  # plot
  dta$optMethod.nobs <- ''
  dta$optMethod.nobs[dta$opMethod == 'Bayesian Optimization'] <- paste('BO\nnObs =', sum(dta$opMethod == "Bayesian Optimization"))
  dta$optMethod.nobs[dta$opMethod == 'Random Optimization'] <- paste('RO\nnObs =', sum(dta$opMethod == "Random Optimization"))

  p <- (ggplot(dta, aes(x = optMethod.nobs, y = BV_mean))
        + geom_boxplot())

  # 2 by 2 comparisons
  dta_BO <- dta[dta$opMethod == 'Bayesian Optimization',]
  dta_RO <- dta[dta$opMethod == 'Random Optimization',]

  # boHigher <- rep(NA, nrow(dta_BO) * nrow(dta_RO))
  boHigher <- c()
  for (i in seq(nrow(dta_BO))) {
    boHigher <- c(boHigher, as.numeric(dta_BO[i, 'BV_mean']) > dta_RO$BV_mean)
  }

  propBoHigher <- sum(boHigher) / length(boHigher)

  p <- p + labs(x = "Optimization method",
                y = "Mean breeding value",
                title = 'Comparison of the optimized breeding scheme from BO and RO',
                subtitle = paste0('Scenario: ',
                               'he = ', unique(dta$he),
                               '; nGen = ', unique(dta$nGen),
                               '; B/gen = ', unique(dta$budget/dta$nGen),
                               '\nBO is better in ',
                               round(propBoHigher*100, 2), '% of the time')
  )

  file <- file.path(outDir, unique(dta$scenario), 'resRepBoxPlot',
                    paste0(unique(dta$scenario), '_boxPlot_repetition_', unique(dta$id)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}





getPropBObetter <- function(dta) {
  # 2 by 2 comparisons
  dta_BO <- dta[dta$opMethod == 'Bayesian Optimization',]
  dta_RO <- dta[dta$opMethod == 'Random Optimization',]
  out <- NULL

  if (nrow(dta_BO) > 1) {
    # t.test
    tt <- t.test(dta_BO$BV_mean, dta_RO$BV_mean)
    tt$p.value

    # # t.test bayes
    # ttbayes <- ttestBF(dta_BO$BV_mean, dta_RO$BV_mean)
    # ttbayes <- ttestBF(dta_BO$BV_mean, dta_RO$BV_mean, posterior = TRUE, iterations = 10000)
    # ttbayes@bayesFactor

    # boHigher <- rep(NA, nrow(dta_BO) * nrow(dta_RO))
    boHigher <- c()
    for (i in seq(nrow(dta_BO))) {
      boHigher <- c(boHigher, as.numeric(dta_BO[i, 'BV_mean']) > dta_RO$BV_mean)
    }
    boHigher_prop <- sum(boHigher) / length(boHigher)
    names(boHigher_prop) <- unique(dta$scenario)

    prior = 0
    out <- list(
      prop = boHigher_prop,
      nSuccess = sum(boHigher),
      nFailure = length(boHigher) - sum(boHigher),
      t.test.pVal = tt$p.value
    )

    # out$p_betaSup50 = 1 - pbeta(0.5,
    #                             prior + out$nSuccess,
    #                             prior + out$nFailure)
    # out$HDP_95 = HDP_beta(prior + out$nSuccess,
    #                       prior + out$nFailure)
    out$nSuccess.cor <- out$prop*nrow(dta)/2
    out$nFailure.cor <- (1-out$prop)*nrow(dta)/2
    out$p_betaSup50 = 1 - pbeta(0.5,
                                out$nSuccess.cor,
                                out$nFailure.cor)
    out$HDP_95 = HDP_beta(out$nSuccess.cor,
                          out$nFailure.cor)
  }
  out
}


HDP_beta <- function(shape1, shape2, a = 0.95){
  mode <- (shape1-1)/(shape1+shape2-2)

  getP2 <- function(p1) {
    q <- min(pbeta(p1, shape1, shape2) + a, 1)
    qbeta(q, shape1, shape2)
  }

  p1 <- uniroot(function(p1){
    dbeta(p1,shape1, shape2) - dbeta(getP2(p1), shape1, shape2)
  }, c(0, mode))$root

  HDP <- c(p1, getP2(p1))

  # prob <- pbeta(HDP, shape1, shape2)
  # prob <- prob[2] - prob[1]
  #
  # dbeta(HDP, shape1, shape2)
  HDP
}





plot_boxPlot_optParams <- function(dta, outDir) {

  dta <- dplyr::select(dta,
                       i, iHomo, bRep, phenoFreq,
                       opMethod, scenario)
  dta <- distinct(dta)
  dta$opMethod <- as.factor(dta$opMethod)
  dta$scenario <- as.factor(dta$scenario)

  # order scenarios
  dta$order <- NA
  dta[dta$scenario == 'he_0.3_nGen_5_b_200', 'order'] <- 1
  dta[dta$scenario == 'he_0.3_nGen_5_b_600', 'order'] <- 2
  dta[dta$scenario == 'he_0.3_nGen_10_b_200', 'order'] <- 3
  dta[dta$scenario == 'he_0.3_nGen_10_b_600', 'order'] <- 4
  dta[dta$scenario == 'he_0.7_nGen_5_b_200', 'order'] <- 5
  dta[dta$scenario == 'he_0.7_nGen_5_b_600', 'order'] <- 6
  dta[dta$scenario == 'he_0.7_nGen_10_b_200', 'order'] <- 7
  dta[dta$scenario == 'he_0.7_nGen_10_b_600', 'order'] <- 8
  dta$scenario <- reorder(dta$scenario, dta$order) # reorder factor levels

  # add number of observation
  n_rep <- dta %>%
    group_by(scenario, opMethod, order) %>%
    summarise(nRep = n()) %>%
    spread(key = opMethod, value = nRep)
  n_rep$scenario_nrep <- paste0(letters[n_rep$order], ': ', n_rep$scenario, '\nnObs: ',
                                n_rep$`Bayesian Optimization`,
                                '\n')
  # n_rep$scenario_nrep <- paste0(n_rep$order, ': ', n_rep$scenario, '\nnObs, ',
  #                               paste('BO:', n_rep$`Bayesian Optimization`,
  #                                     'RO:', n_rep$`Random Optimization`),
  #                               '\n')
  dta <- full_join(dta, n_rep)
  dta$scenario_nrep <- reorder(dta$scenario_nrep, dta$order) # reorder factor levels

  # use paper parameter names
  dta$pheno_p <- dta$phenoFreq
  dta$iInit <- dta$iHomo
  dta$order <- as.character(dta$or)


  # remove random opt
  dta <- filter(dta, opMethod == 'Bayesian Optimization')

  plotList <- lapply(c('i', 'iInit', 'bRep', 'pheno_p'), function(param){
    p <- ggplot(dta, aes(x = letters[as.numeric(dta$order)], y = dta[, param], col = scenario_nrep))
    p <- p + geom_boxplot()
    p <- p + labs(y = param,
                  title = paste0('Boxplot of the optimized parameter "', param,
                                 '"'),
                  x = 'Scenario',
                  col = 'Scenario')



    # save plot
    file <- file.path(outDir, 'optimizedParams',
                      paste0('optParam_boxplot_', param))
    suppressWarnings(dir.create(dirname(file), recursive = TRUE))
    ggsave(file, p)
    p
  })

  # # plot with both BO and RO :
  #
  # plotList <- lapply(c('i', 'iInit', 'bRep', 'pheno_p'), function(param){
  #   p <- ggplot(dta, aes(x = opMethod, y = dta[, param], col = scenario_nrep))
  #   p <- p + geom_boxplot()
  #   p <- p + labs(y = param,
  #                 title = paste0('Boxplot of the optimized parameter "', param,
  #                                '"'),
  #                 col = 'Scenario')
  #   # ,
  #   #               subtitle = paste0(
  #   #                 "Number of observations:\n",
  #   #                 paste0(n_rep$scenario,
  #   #                        ': BO ', n_rep$`Bayesian Optimization`,
  #   #                        ', RO ', n_rep$`Random Optimization`, collapse = '\n')
  #   #                 ))
  #
  #
  #   # anova
  #   mod_bo <- lm(paste0(param, ' ~ scenario'),
  #                data = dta[dta$opMethod == 'Bayesian Optimization',])
  #   summary(mod_bo)
  #   anova(mod_bo)
  #   mod_ro <- lm(paste0(param, ' ~ scenario'),
  #                data = dta[dta$opMethod == 'Random Optimization',])
  #   summary(mod_ro)
  #   anova(mod_ro)
  #
  #   # save plot
  #   file <- file.path(outDir, 'optimizedParamsBoxplot',
  #                     paste0('optParam_boxplot_', param))
  #   suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  #   ggsave(file, p)
  #   p
  # })


  # dta_BO <- dta[dta$opMethod == 'Bayesian Optimization',]
  # dta_RO <- dta[dta$opMethod == 'Random Optimization',]
  #
  #
  # p_i <- ggplot(dta, aes(x = opMethod, y = i, col = scenario))
  # p_i <- ggplot(dta, aes(x = scenario, y = i, col = opMethod, groups = opMethod))
  # # p_i <- ggplot(dta_BO, aes(x = scenario, y = iHomo))
  # p_i + geom_boxplot()
  # mod <- lm(i ~ scenario, data = dta_BO)
  # summary(mod)
  # anova(mod)
  #
  # p_iHomo <- ggplot(dta, aes(x = opMethod, y = iHomo, col = scenario))
  # p_iHomo <- ggplot(dta, aes(x = scenario, y = iHomo, col = opMethod, groups = opMethod))
  # # p_iHomo <- ggplot(dta_BO, aes(x = scenario, y = iHomo))
  # p_iHomo + geom_boxplot()
  # mod <- lm(iHomo ~ scenario, data = dta_BO)
  # summary(mod)
  # anova(mod)
  #
  # p_bRep <- ggplot(dta, aes(x = opMethod, y = bRep, col = scenario))
  # p_bRep <- ggplot(dta, aes(x = scenario, y = bRep, col = opMethod, groups = opMethod))
  # # p_bRep <- ggplot(dta_BO, aes(x = scenario, y = bRep))
  # p_bRep + geom_boxplot()
  # mod <- lm(bRep ~ scenario, data = dta_BO)
  # summary(mod)
  # anova(mod)
  #
  # p_pF <- ggplot(dta, aes(x = opMethod, y = phenoFreq, col = scenario))
  # p_pF <- ggplot(dta, aes(x = scenario, y = phenoFreq, col = opMethod, groups = opMethod))
  # # p_pF <- ggplot(dta_BO, aes(x = scenario, y = phenoFreq))
  # p_pF + geom_boxplot()
  # mod <- lm(phenoFreq ~ scenario, data = dta_BO)
  # summary(mod)
  # anova(mod)
  plotList
  }






plot_PCA_optParams <- function(dta, outDir) {

  # remove random opt
  dta <- filter(dta, opMethod == 'Bayesian Optimization')

  dta <- dplyr::select(dta,
                       i, iHomo, bRep, phenoFreq,
                       opMethod, scenario)
  dta <- distinct(dta)
  dta$opMethod <- as.factor(dta$opMethod)
  dta$scenario <- as.factor(dta$scenario)

  # order scenarios
  dta$order <- NA
  dta[dta$scenario == 'he_0.3_nGen_5_b_200', 'order'] <- 1
  dta[dta$scenario == 'he_0.3_nGen_5_b_600', 'order'] <- 2
  dta[dta$scenario == 'he_0.3_nGen_10_b_200', 'order'] <- 3
  dta[dta$scenario == 'he_0.3_nGen_10_b_600', 'order'] <- 4
  dta[dta$scenario == 'he_0.7_nGen_5_b_200', 'order'] <- 5
  dta[dta$scenario == 'he_0.7_nGen_5_b_600', 'order'] <- 6
  dta[dta$scenario == 'he_0.7_nGen_10_b_200', 'order'] <- 7
  dta[dta$scenario == 'he_0.7_nGen_10_b_600', 'order'] <- 8
  dta$scenario <- reorder(dta$scenario, dta$order) # reorder factor levels
  dta$iInit <- dta$iHomo
  dta$pheno_p <- dta$phenoFreq

  # PCA
  pcaDta <- dplyr::select(dta, i, iInit, bRep, pheno_p, scenario)
  pcaRes <- PCA(pcaDta,
                scale.unit = TRUE,
                quali.sup = 5,
                ncp = 4,
                graph = T)
  scenarioGravityCenter <- as.data.frame(pcaRes$quali.sup$coord)
  scenarioGravityCenter$scenario <- row.names(scenarioGravityCenter)

  # add PCs to the data
  dta <- mutate(dta,
                pca1 = pcaRes$ind$coord[, 1],
                pca2 = pcaRes$ind$coord[, 2])

  plotDta <- dta
  plotDta <- filter(dta, scenario %in% c('he_0.3_nGen_5_b_200',
                                         'he_0.7_nGen_10_b_600'))
  plotGravity <- scenarioGravityCenter
  plotGravity <- filter(scenarioGravityCenter, scenario %in% c('he_0.3_nGen_5_b_200',
                                             'he_0.7_nGen_10_b_600'))

  # density plot
  p <- ggplot(plotDta, aes(x = pca1,
                           y = pca2,
                           col = scenario))
  p <- p + geom_vline(xintercept = 0)
  p <- p + geom_hline(yintercept = 0)
  p <- p + geom_density_2d(aes())
  # p <- p + geom_point(alpha = 1/10) # add all points
  # add variables arrows
  p <- p + geom_point(data = plotGravity,
                      size = 4,
                      aes(x = Dim.1,
                          y = Dim.2,
                          col = scenario))
  p <- p + annotate("segment",
             x = 0,
             xend = pcaRes$var$coord[,1],
             y = 0,
             yend = pcaRes$var$coord[,2],
             size = 0.5,
             arrow = arrow(length = unit(0.2, "cm")))
  p <- p + annotate("text",
             x = pcaRes$var$coord[,1],
             y = pcaRes$var$coord[,2],
             label = row.names(pcaRes$var$coord),
             hjust = 'outward',
             vjust = 'outward')

  p <- p + labs(title = paste0("Optimized parameters's PCA visualisation"),
                subtitle = "Points correspond to the center of gravity of each scenario."
                )

  # same range for x and y axes
  axisLim <- c(min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1],
                   ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1]),
               max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2],
                   ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2])
  )
  p <- p + xlim(axisLim) + ylim(axisLim)
  p

  # save plot
  file <- file.path(outDir, 'optimizedParams',
                    'optParam_pca')
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p
}

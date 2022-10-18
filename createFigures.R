# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Script to create article's figures

#### PACKAGES ####
library(ggplot2)
library(tidyverse)
library(pbapply)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

#### OPTIONS ####
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_bw())
theme_update(text = element_text(size = 15)) # use 20 for eps
set.seed(2022)
nClust = 8

ggsave <- function(file, ...) {
  file <- paste0(file, '.jpg') # .svg
  ggplot2::ggsave(filename = file, width = 20, height = 15, units = "cm", dpi = 300)
}

# get functions
source('src/resultsAnalysis/resultsAnalysisFunctions.R')

#### CODE ####
aggregatedFiles <- c('aggregatedResults/16x32_repetitions.rds',
                     'aggregatedResults/1024-15itersx1_repetitions.rds')
aggregatedFile <- 'aggregatedResults/16x32_repetitions.rds'
# aggregatedFile <- 'aggregatedResults/1024-15itersx1_repetitions.rds'

for (aggregatedFile in aggregatedFiles) {

  # create output folder
  outDir <- paste0('figures/', tools::file_path_sans_ext(basename(aggregatedFile)))
  suppressWarnings(dir.create(outDir))

  # read aggregateresults
  aggResults <- readRDS(aggregatedFile)


  ## Number of optimization repetition for each scenario ----
  nRep_df <- get_nRep_by_scenario(aggResults)
  write.csv(nRep_df, paste0(outDir, '/n_repetitions.csv'))


  allData <- merge_allData(aggResults)


  # 1. Optimization process: ----

  # summary graph per scenario ----
  # extract all optimization progress in one data.frame
  all_optim_data <- allData$optimDta

  # split data by scenarios:
  # we get one data.frame per scenario
  list_byScenario <- all_optim_data %>%
    group_by(scenario) %>%
    group_map(mutate)

  ## Box plot, number of iteration ----
  # Draw box plot of maximum number of optimization iteration for each scenario
  # this plot is interesting when the optimization stopping criterion was
  # a specific time resources, in such case, the number of iteration is not
  # defined.
  cat("plot nIter\n")
  invisible(pblapply(list_byScenario, boxPlot_nIter, outDir = outDir, cl = nClust))

  ## Box plot, cumulative maximum ----
  # Draw box plot of the cumulative maximum at each iteration for each scenario
  cat("plot cumMax boxplot\n")
  invisible(pblapply(list_byScenario, boxPlot_cumMax, outDir = outDir, cl = nClust))

  ## Line plot, cumulative maximum ----
  # Draw line plot of the cumulative maximum at each iteration for each scenario
  cat("plot cumMax lines\n")
  invisible(pblapply(list_byScenario, lines_cumMax, outDir = outDir, cl = nClust))


  # draw optimization progress ----
  # select a subset of optimization (to avoid creating too many plots)
  cat('select subset of optimization\n')
  list_byScenario_byOpt <- pblapply(list_byScenario, function(dta) {
    list_byOpt <- split_by_optimization(dta)
    selectedId <- sort(sample(length(list_byOpt),
                              min(16, length(list_byOpt))))
    list_byOpt[selectedId]
  })
  list_byOpt <- unlist(list_byScenario_byOpt, recursive = FALSE)

  cat('plot optimization progress\n')
  invisible(pblapply(list_byOpt, plot_optimization_progress,
                     outDir = outDir, cl = nClust))

  cat('plot optimization progress PCA\n')
  invisible(pblapply(list_byOpt, plot_optimization_pca,
                     outDir = outDir, cl = nClust))

  # get most representative optim
  mostRepresentativeDta <- pblapply(list_byScenario,
                                    getMostRepresentativeOptim,
                                    cl = nClust)
  invisible(pblapply(mostRepresentativeDta, plot_optimization_progress,
                     outDir = outDir, cl = nClust))
  invisible(pblapply(mostRepresentativeDta, plot_optimization_pca,
                     outDir = outDir, cl = nClust))










  # 2. Optimization result evaluation: ----

  # extract all optimization results evaluation in one data.frame
  # all_optResultEval_data <- merge_all_optResultsEval_data(aggResults)
  all_optResultEval_data <- allData$resRepDta

  # empirical cumulative distribution ----
  # split data by scenarios:
  # we get one data.frame per scenario
  list_byScenario_optResEval <- all_optResultEval_data %>%
    group_by(scenario) %>%
    group_map(mutate)

  cat('plot ECD\n')
  invisible(pblapply(list_byScenario_optResEval, plot_empCumDist, outDir = outDir, cl = nClust))


  BObetter <- pblapply(list_byScenario_optResEval, getPropBObetter, cl = nClust)


  cat('select subset of optimization\n')
  list_byScenario_byOpt_optResEval <- pblapply(list_byScenario_optResEval, function(dta) {
    list_byOpt <- split_by_optimization(dta)
    selectedId <- sort(sample(length(list_byOpt),
                              min(16, length(list_byOpt))))
    list_byOpt[selectedId]
  })
  list_byOpt_optResEval <- unlist(list_byScenario_byOpt_optResEval, recursive = FALSE)

  cat('plot boxplot repeted results\n')

  invisible(pblapply(list_byOpt_optResEval, plot_boxPlot_optResEval, outDir = outDir, cl = nClust))


  mostRepresentativeDta_repetedSim <- mapply(function(repetedSimDta, mostRepDta){
    out <- repetedSimDta[repetedSimDta$optSeed %in% unique(mostRepDta$optSeed),]
    out$id <- paste0(
      'bo_', unique(out[out$opMethod == 'Bayesian Optimization', 'repId']), '-',
      'ro_', unique(out[out$opMethod == 'Random Optimization', 'repId']))
    out
  },
  repetedSimDta = list_byScenario_optResEval,
  mostRepDta = mostRepresentativeDta, SIMPLIFY = FALSE)

  invisible(pblapply(mostRepresentativeDta_repetedSim, plot_boxPlot_optResEval, outDir = outDir, cl = nClust))

  BObetter2 <- pblapply(mostRepresentativeDta_repetedSim, getPropBObetter, cl = nClust)

  for (i in seq(length(BObetter2))) {
    if (!is.null(BObetter2[[i]])) {

      plotDta <- data.frame(x = seq(0, 1, length.out = 200)) %>%
        dplyr::mutate(y1 = dbeta(x, BObetter2[[i]]$nSuccess.cor, BObetter2[[i]]$nFailure.cor))
      # plot_ly(
      #   type = "scatter",
      #   mode = "lines",
      #   data = plotDta,
      #   x = ~ x,
      #   y = ~ y1,
      #   name = "y1"
      # ) %>%
      #   layout(hovermode = 'compare')
      plot(plotDta$x, plotDta$y1)
    }
  }


  # 3. optimized parameters values ----
  cat('plot optimized parameters\n')
  invisible(plot_boxPlot_optParams(all_optResultEval_data, outDir))
  invisible(plot_PCA_optParams(all_optResultEval_data, outDir))

}

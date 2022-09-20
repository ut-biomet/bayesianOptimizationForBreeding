# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Definition of the selection an mating functions

#' Select and mate individuals
#'
#' @param pop breedSimulatR's population
#' @param model a list containing two elements: `inter`, the intercept, and `qtnEff` the estimated effects of all the qtn.
#' @param nSel number of selected individuals
#' @param nNewTot total number of new progenies to create
#' @param basename prefix of the name of the new individuals
#'
#' @return vector of the name of the selected individuals
selectMate <- function(pop,
                       model,
                       nSel,
                       nNewTot,
                       basename) {

  ## selection ----
  selectedInds <- selectGS(pop = pop,
                           n = nSel,
                           model = model)

  ## crossing ----
  crosses <- matInds(pop = pop,
                     selected = selectedInds,
                     nNewTot,
                     qtnEff = model$qtnEff,
                     basename)

  crosses
}


# ##### true genetic value selection ####
# selectGV <- function(pop, n, trait) {
#   gv <- trait$gv(pop)
#   names(pop$inds)[order(gv, decreasing = TRUE)][1:n]
# }

#' Select individuals based on the a given prediction model
#'
#' @param pop breedSimulatR's population
#' @param n number of selected individuals
#' @param model a list containing two elements: `inter`, the intercept, and `qtnEff` the estimated effects of all the qtn.
#'
#' @return vector of the name of the selected individuals
selectGS <- function(pop, n, model) {
  geno <- pop$genoMat[, names(model$qtnEff)]
  predGV <- model$inter + geno %*% matrix(model$qtnEff, ncol = 1)
  names(pop$inds)[order(predGV, decreasing = TRUE)][1:n]
}


#### Mating functions ####
#' "TSP" mating
#'
#' @param pop breedSimulatR's population
#' @param selected vector of the name of the selected individuals
#' @param nTot  total number of new progenies to create
#' @param basename prefix of the name of the new individuals
#' @param qtnEff estimated qtn's effects
#'
#' @return a cross table, see: `breedSimulatR::makeCrosses`
matInds <- function(pop,
                    selected,
                    nTot,
                    basename,
                    qtnEff) {


  stopifnot(length(qtnEff) == ncol(pop$genoMat[selected,]))

  # fast return
  if (length(selected) == 1) {
    crossTable <- data.frame(ind1 = selected,
                             ind2 = selected,
                             n = nTot,
                             names = basename)
    return(crossTable)
  } else if (length(selected) == 2) {
    crossTable <- data.frame(ind1 = selected[1],
                             ind2 = selected[2],
                             n = nTot,
                             names = basename)
    return(crossTable)
  }

  if (length(selected) == 3) {
    ind1 <- selected
    ind2 <- selected[c(2,3,1)]
    nMat <- 3
  } else {
    # weighted distances between individuals
    dw <- dist(pop$genoMat[selected,]  %*% diag(qtnEff))

    # find "longest" path
    tsp <- TSP::as.TSP(1/(dw + 10^-3)) # avoid 0 as denominator
    best_tour <- TSP::TOUR(1:TSP::n_of_cities(tsp), tsp = tsp)
    min_length <- TSP::tour_length(best_tour)
    for (i in 1:10) {
      # tour <- TSP::solve_TSP(tsp, method = "nn")
      tour <- TSP::solve_TSP(tsp, method = "two_opt",
                             # tour = tour,
                             two_opt_repetitions = 1000)
      if (TSP::tour_length(tour) < min_length) {
        # if the obtained length is the shorter than the current best
        min_length <- TSP::tour_length(tour) # update the shortest path
        best_tour <- tour
      }
    }
    # show the best crosses
    tmp <- labels(best_tour)
    best.path <- c(tmp, tmp[1])
    nMat <- length(best.path) - 1
    ind1 <- best.path[-length(best.path)]
    ind2 <- best.path[-1]
  }


  nRep <- nTot / nMat
  if (floor(nRep) == nRep) {
    n <- rep(nRep, nMat)
  } else {
    n <- rep(floor(nRep), nMat)
    n[sample(nMat, nTot - sum(n))] <- floor(nRep) + 1 # randomly select some cross to create 1 more offspring
  }

  crossTable <- data.frame(ind1,
                           ind2,
                           n = n,
                           names = paste(basename,
                                         seq(nMat),
                                         sep = "_"))
  crossTable
}


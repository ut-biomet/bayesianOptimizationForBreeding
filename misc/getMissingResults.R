# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Get seed and scenario params of missing results

#### OPTIONS & PACKAGE ####
options(stringsAsFactors = FALSE)
library(tidyverse)

#### CODE ####


# get list of seed scenario params of results
resFolder <- 'output-resultRepetition/1024_reps_15iters_excl/'

files <- list.files(resFolder)
files[1]
splitFile <- strsplit(files, '_')
dta <- as.data.frame(do.call(rbind, splitFile))
dta <- dta[, -c(1, 6, 7)]
colnames(dta) <- c('method', 'nGen', 'budget', 'he', 'seed')
for (colname in c('nGen', 'budget', 'he', 'seed')) {
  dta[, colname] <- regmatches(dta[, colname], regexpr("\\d+(\\.\\d)*", dta[, colname]))
  dta[, colname] <- as.numeric(dta[, colname])
}
dta <- tibble(dta)
dta <- arrange(dta, method, nGen, budget, he, seed)


# get list of planed scenario
dta2 <- read.csv('seedList_all.csv', header = T, )
nrow(dta2)/8

dta2 <- separate(dta2, name, sep = "_", into = c('nGen', 'budget', 'he'))
for (colname in c('nGen', 'budget', 'he')) {
  dta2[, colname] <- regmatches(dta2[, colname], regexpr("\\d+(\\.\\d)*", dta2[, colname]))
  dta2[, colname] <- as.numeric(dta2[, colname])
}

dta2 <- pivot_longer(dta2, cols = c(6,7), names_to = 'method', values_to = 'seed')
dta2[dta2$method == 'repBO', "method"] <- 'bayesOpt'
dta2[dta2$method == 'repRand', "method"] <- 'randOpt'
dta2 <- dta2[, c('method', 'nGen', 'budget', 'he', 'seed')]


dta2 <- arrange(dta2, method, nGen, budget, he, seed)




# get missing scenario
missing <- setdiff(dta2, dta)
View(missing)

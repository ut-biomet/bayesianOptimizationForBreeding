# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# shortcut to run the results

library(dotenv)
load_dot_env()
source("src/utils.R")

SelIntmode = 0 #selection intensity function (slope_i function)
# if SelIntmode = 0, then SelIntSlope=0 and NewIndSlope=0
BORep <- 1
nDiffStartPopulation = 2
# set.seed(1)
# (SeedsForPopulationSetup = floor(runif(nDiffStartPopulation)*1e6))
# digest::digest2int(as.character(seq(nDiffStartPopulation)))
SeedsForPopulationSetup = round(1:nDiffStartPopulation * pi*1e6) %>%
  rep(., each=BORep)


for (StartPop_k in 1:length(SeedsForPopulationSetup)){
t <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

PopSetupSeed = SeedsForPopulationSetup[StartPop_k]

rmarkdown::render('runRepeatOpt.Rmd',
                  output_file = paste0('runRepeatOpt_', t, "-initPop-",StartPop_k,'|',length(SeedsForPopulationSetup)),
                  params = list(PopSetupSeed = PopSetupSeed,
                                SelIntmode = SelIntmode,
                                StartPop_k = StartPop_k))
sendNotification("optimization finished !!!")
}

slackr::slackr_setup()
slackr::slackr_upload(filename = "runRepeatOpt.md")

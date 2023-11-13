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
nDiffStartPopulation = 2
SeedsForPopulationSetup = round(1:nDiffStartPopulation * pi*1e6)

for (startseed in SeedsForPopulationSetup){
t <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
PopSetupSeed = startseed

rmarkdown::render('runRepeatOpt.Rmd',
                  output_file = paste0('runRepeatOpt_', t),
                  params = list(PopSetupSeed = PopSetupSeed,
                                SelIntmode = SelIntmode))
sendNotification("optimization finished !!!")
}

slackr::slackr_setup()
slackr::slackr_upload(filename = "runRepeatOpt.md")

# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# shortcut to run the results

library(dotenv)
load_dot_env()
source("src/utils.R")
nDiffStartPopulation = 2
SeedsForPopulationSetup = round(1:nDiffStartPopulation * pi*1e6)
for (startseed in SeedsForPopulationSetup){
t <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
PopSetupSeed = startseed

rmarkdown::render('runRepeatOpt.Rmd',
                  output_file = paste0('runRepeatOpt_', t),
                  params = list(PopSetupSeed = PopSetupSeed))
sendNotification("optimization finished !!!")
}

slackr::slackr_setup()
slackr::slackr_upload(filename = "runRepeatOpt.md")

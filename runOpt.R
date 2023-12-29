# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# shortcut to run the results
# library(rmarkdown)
library(dplyr)
library(dotenv)
load_dot_env()
source("src/utils.R")

SelIntmode = 1 #selection intensity function (slope_i function)
Fixed_nNew = TRUE
Fixed_SelInt = FALSE

if (SelIntmode == 0) { # Control : all parameters are constant from the 2nd generation
  Fixed_nNew = TRUE
  Fixed_SelInt = TRUE
}


# if SelIntmode = 0, then SelIntSlope=0 and NewIndSlope=0
nDiffStartPopulation = 1
BORep <- 1

SeedsForPopulationSetup <- round(1:nDiffStartPopulation * pi*1e6) %>%
  rep(., each=BORep)
SeedsForBO <- digest::digest2int(as.character(1:length(SeedsForPopulationSetup))) # set 1 different seed for each BO

for (StartPop_k in 1:length(SeedsForPopulationSetup)){
  t <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  PopSetupSeed = SeedsForPopulationSetup[StartPop_k]
  BOseed <- SeedsForBO[StartPop_k]
  rmarkdown::render('runRepeatOpt.Rmd',
                    output_file = paste0('runRepeatOpt_', t, "-initPop-",floor(StartPop_k/BORep),'|',
                                         floor(length(SeedsForPopulationSetup)/BORep)),
                    params = list(PopSetupSeed = PopSetupSeed,
                                  SelIntmode = SelIntmode,
                                  Fixed_nNew = Fixed_nNew,
                                  Fixed_SelInt = Fixed_SelInt,
                                  StartPop_k = floor(StartPop_k/BORep),
                                  BOseed = BOseed))
  # sendNotification("optimization finished !!!")
}

# slackr::slackr_setup()
# slackr::slackr_upload(filename = "runRepeatOpt.md")



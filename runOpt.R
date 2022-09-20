# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# shortcut to run the results

library(dotenv)
load_dot_env()
source("src/utils.R")
t <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
rmarkdown::render('runRepeatOpt.Rmd',
                  output_file = paste0('runRepeatOpt_', t))
sendNotification("optimization finished !!!")

slackr::slackr_setup()
slackr::slackr_upload(filename = "runRepeatOpt.md")

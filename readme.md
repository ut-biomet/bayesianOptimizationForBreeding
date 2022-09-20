# Bayesian optimization of breeding schemes

## Introduction

This repository contain the `R` code related the the study "Bayesian Optimization for breeding" submitted to *Frontiers in Plant Science* by DIOT Julien and IWATA Hiroyoshi. This study presents the usage of Bayesian optimization for optimizing Breeding schemes and breeding simulation thanks to the `R` packages [`breedSimulatR`](https://github.com/ut-biomet/breedSimulatR) and [`mlrMBO`](https://mlrmbo.mlr-org.com/).

This readme is here to help reproduce the results or adapt the method on different use case.

### Repository structure

- [`runRepeatOpt.Rmd`](./runRepeatOpt.Rmd): **Main code running repeated Bayesian and Random optimization.**
- [`src`](./src): Folder containing all the `R` scripts defining the functions used in [`runRepeatOpt.Rmd`](./runRepeatOpt.Rmd).
- [`data`](./data): Location of the raw genotype data of the initial breeding population.
- [`simSetups`](./simSetups):  Location of the saved "simulation setups" files. Those files might be helpful to run another optimization using the same simulation setup (eg. same breeding simulation function, same heritability, same initial population...)
- [`optSetups`](./optSetups): Location of the saved "optimization setups" files. Those files might be helpful to run another optimization using the same setup.
- [`output-optimization`](./output-optimization): Location of the optimization results.
- [`output-resultRepetition`](./output-resultRepetition): Location of the Breeding schemes replication done after an optimization.
- [`aggregatedResults`](./aggregatedResults): Folder containing the Bayesian and Random optimization results aggregated in singles `.rds` files.
- [`aggregateResults.R`](./aggregateResults.R): Script aggregating the results of the Bayesian and Random optimizations in singles `.rds` files. 
- [`createFigures.R`](./createFigures.R): Script generating the figures presented in the paper.
- [`figures`](./figures): Location of the figures generatied by [`createFigures.R`](./createFigures.R).
- [`misc`](./misc): Miscellaneous folder maid for files without specific interest for this repo but which might not worth deleting:
  - [`misc/errorMsg.txt`](./misc/errorMsg.txt)an example of error messages that can appear in when doing the bayesian optimization
  - [`misc/getMissingResults.R`](./misc/getMissingResults.R)a script that can get the id of the optimization runs which have generated the errors (based on the missing results file).
  - [`misc/seedList_all.csv`](./misc/seedList_all.csv): list of the ids of all the results.
  - [`misc/singleSimulation.R`](./singleSimulation.R): Script to Simulate 1 breeding scheme according to breeding scheme parameters return by the optimization.

<!-- - [`bayesianOptimizationForBreeding.Rproj`](./bayesianOptimizationForBreeding.Rproj): -->
<!-- - [`LICENSE`](./LICENSE): -->
<!-- - [`readme.md`](./readme.md): -->
<!-- - [`readme.Rmd`](./readme.Rmd): -->
<!-- - [`renv`](./renv): -->
<!-- - [`renv.lock`](./renv.lock): -->
<!-- - [`runOpt.R`](./runOpt.R): -->
<!-- - [`runRepeatOpt_2022-04-26_14-37-31.html`](./runRepeatOpt_2022-04-26_14-37-31.html):
- [`runRepeatOpt_2022-04-26_14-37-31.md`](./runRepeatOpt_2022-04-26_14-37-31.md):
- [`runRepeatOpt_2022-04-26_14-37-36.html`](./runRepeatOpt_2022-04-26_14-37-36.html):
- [`runRepeatOpt_2022-04-26_14-37-36.md`](./runRepeatOpt_2022-04-26_14-37-36.md):
- [`runRepeatOpt_dec2021.md`](./runRepeatOpt_dec2021.md): -->
<!-- - [`utils`](./utils): -->

## Install dependencies

All computation are done with the [`R`](https://www.r-project.org/) (v4.1.1) language and uses the [`renv`](https://rstudio.github.io/renv/index.html) package to manage the library dependencies.

To install the necessary libraries, open a R console in the project folder and run:

```r
renv::restore()
```

## Optimization workflow

The main optimization workflow is detailed in the file [`runRepeatOpt.Rmd`](./runRepeatOpt.Rmd). It is a R-markdown file which can be run from R-Studio directly or with the shell command:

```sh
R -e "rmarkdown::render('runRepeatOpt.Rmd')"
```

> Note: Be aware that the calculation can be long.

The workflow can be summarized as follow:

1. Create the "Simulation setup". This setup contain all the necessary information to run a breeding simulation (eg. genotype of the initial population, marker effects, heritability, budget...).
2. Create the "Optimization setup". This setup contain all the necessary information to run an optimization (eg. number of iterations, kernel, acquisition function...).
3. Run the optimizations and repeat the breeding simulation with the optimized parameters:
   - Bayesian optimization
   - Repeat the breeding simulation using the "bayesian optimized" parameters
   - Random optimization
   - Repeat the breeding simulation using the "random optimized" parameters

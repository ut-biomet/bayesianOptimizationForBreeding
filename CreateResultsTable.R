require(dbplyr)
require(stringr)


#' Create table from all the Result files in a directory
#' the result files names have to fit a certain format '_ParameterName-Value_'
#' Or you have to give the list of names
#'
#' @param repositoryPath
#' @param paramSep separator string
#' @param valueSep
#' @param fileExtension
#' @param ParamNames list of the parameter names contained in the FileNames
CreateResFileTable <- function(repositoryPath,
                        paramSep = '_',
                        valueSep = '-',
                        fileExtension = '.rds',
                        ParamNames = NULL){

  # browser()
  # Explore Repository
  paths <- list.files(repositoryPath, recursive = T)

  # keep only file with the right file extension  #### ^.*/dossier/.*\.rds$
  paths <- paths[grepl(pattern = paste0(fileExtension, "$"), x = paths)]

  # Create table
  ResFileTable <- data.frame(Path = unlist(lapply(paths, function(path) stringr::str_extract(path, '^.*/'))),
                             FileNames = unlist(lapply(paths, function(path) gsub('^.*/','',path))))

  # Extract FileNames information
  ParamInfo <- lapply(ResFileTable$FileNames, function(fileName)
    unlist(strsplit(x = gsub(pattern = fileExtension, "", x = fileName), split = paramSep))) %>%
    do.call(rbind,.)

  # Get Parameter Names
  if (!is.null(ParamNames)){
    if (length(ParamNames) != ncol(ParamInfo)){
      errorCondition("ParamNames not the right size (", length(ParamNames),"instead of", ncol(ParamInfo),")")
    }
  } else {
    ParamNames <- unlist(lapply(ParamInfo[1,], function(param)
      gsub(pattern = "[-.0-9]|FALSE|TRUE","",x = param, ignore.case = TRUE) %>%
        ifelse(. == '', "NoName",.)))
  }
  colnames(ParamInfo) <- ParamNames

  # Remove Text from data
  ParamInfo <- lapply(colnames(ParamInfo), function(colname)
    unlist(lapply(ParamInfo[,colname], function(cell) unlist(strsplit(cell,split = '-')) %>%
      tail(1)))) %>%
    do.call(cbind,.) %>%
    as.data.frame
  colnames(ParamInfo) = ParamNames

  # change variable types
  for (colname in ParamNames) {
    if (all(grepl(pattern = 'TRUE|FALSE', x = ParamInfo[, colname]))) {
      ParamInfo[, colname] <- as.logical(ParamInfo[, colname])
    } else if (any(grepl(paste0(c("[", letters, ']'), collapse = ""), x = ParamInfo[, colname]))) {
      ParamInfo[, colname] <- as.character(ParamInfo[, colname])
    } else {
      ParamInfo[, colname] <- as.numeric(ParamInfo[, colname])
    }
  }
  ParamInfo <- as.data.frame(ParamInfo)
  #Gather folder and FileNames table
  ResFileTable <- cbind(ParamInfo, ResFileTable)
  return(ResFileTable)
}


#' Keep only the caracteristics in the Table
#' check if the sample sizes are identical
#' @param ResFileTable output of the CreateResFileTable function
#' @param ParamNames Names of the column to keep
#' @param UnifSampleSize Does the sample sizes have to be identical (TRUE/FALSE)
#' @param VarsOfInterest Names of the column that need to have the same sample sizes
CleanResFileTable <- function(ResFileTable,
                              ParamNames = NULL,
                              BalanceResCheck = TRUE,
                              VarsOfInterest = NULL) {
  # browser()
  # Remove unused columns
  if (is.null(ParamNames)){
    #remove the column if all values are identical
    ResFileTable <- unlist(lapply(colnames(ResFileTable), function(colname)
      length(table(ResFileTable[,colname])) > 1)) %>%
      ResFileTable[,.]
    ParamNames = colnames(ResFileTable)
  } else {
    if (!all(ParamNames %in% colnames(ResFileTable))) {
      errorCondition("ParamNames contains parameter that does not exist in ResFileTable.")
    }
    # keep path information in ResFileTable if not mentioned
    ParamNames <- unique(c(ParamNames, "Path","FileNames") )
    ResFileTable <- ResFileTable[,ParamNames]
  }

  if (BalanceResCheck) {
    if (is.null(VarsOfInterest)) {
      VarsOfInterest = ParamNames
    }
    # Checking if the sample sizes are the same for all parameter combinations
    BalanceResCheck = (length(unique(table(ResFileTable[,VarsOfInterest]))) == 1)
    cat("BalanceResCheck : ", BalanceResCheck, "\n")
  }
  return(ResFileTable)
}


#' Get optimised parameter form a list of result files
#' @param FilePaths.list
#' @param GetInitPopHash
#' @param KeepOptValue Boolean TRUE: keep the value given by the Bayesian Optimization
#'                             FALSE : get the mean value from the output-repetition folder
#' @param IsRepFile adapt the function if the repetition folder is used
GetOptimizedParams <- function(FilePaths.list,
                               GetInitPopHash = FALSE,
                               KeepOptValue = FALSE,
                               IsRepFile = FALSE) {
  lapply(FilePaths.list, function(FilePath){
    Resfile = readRDS(FilePath)
    if (IsRepFile) {optimizedParams <- Resfile$params[!(names(Resfile$params) %in% c("seeds", "nRep", "nCpus"))]
    } else {optimizedParams <-  Resfile$bestPoint}

    # Get BreedingImprovement Value
    if (IsRepFile) {
      BVmean <- Resfile$results$BV_mean[1]
      if (!KeepOptValue) BVmean <- mean(Resfile$results$BV_mean)
    } else {
      BVmean = Resfile$bestPoint_pred
      if (!KeepOptValue) {
        RepFile <- gsub("output-optimization/", "output-resultRepetition/", FilePath) %>%
          gsub(pattern = "boRes_.*$",replacement = "",x=.) %>%
        paste0(.,list.files(.)) %>%
          .[grepl('bayesOpt',.)]
        BVmean <- mean(readRDS(RepFile)$results$BV_mean)
      }
    }

    #Set unoptimized values to 0
    if (is.null(optimizedParams$NewIndSlope)) optimizedParams = c(optimizedParams, NewIndSlope = 0)
    if (is.null(optimizedParams$SelIntSlope)) optimizedParams = c(optimizedParams, SelIntSlope = 0)
    if (IsRepFile) Resfile <- Resfile$optResults
    if (GetInitPopHash) {
      initPop = Resfile$simSetup$fixedParams$initPop
      optimizedParams$InitPopHash <-  digest::digest2int(digest::digest(initPop))
      Pheno <- Resfile$simSetup$fixedParams$phenotyper$trial(initPop)
      optimizedParams$InitPopBV <- max(Pheno$data$trait1)
    }
    return(unlist(c(BVmean = BVmean,optimizedParams)))
    }) %>%
    do.call(rbind,.) %>%
    as.data.frame
}


#' Compare the improvement
#' @param FilePaths.list
CompareImprovTable <- function(FilePaths.list) {
  lapply(FilePaths.list, function(FilePath){
    Resfile = readRDS(FilePath)
    return(data.frame(SelIntmode = rep(Resfile$simSetup$fixedParams$SelIntmode,length(Resfile$results$BV_mean)),
                      InitPopHash = digest::digest2int(digest::digest(Resfile$simSetup$fixedParams$initPop)),
                      FixednNew = is.null(Resfile$simSetup$optParams$pars$NewIndSlope) %>% rep(length(Resfile$results$BV_mean)),
                      BOiter = Resfile$results$dob,
                      BV_mean = Resfile$results$BV_mean))
  }) %>%
    do.call(rbind,.) %>%
    as.data.frame
}








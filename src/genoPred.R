
##### Genetic prediction selection ####
# GS model (REML)
#' Create a genomic prediction model based on a Ridge regression
#'
#' @param pheno phenotypic values
#' @param geno genotypic values
#' @param traitName trait name
#'
#' @return a list containing two elements: `inter`, the intercept, and `qtnEff` the estimated effects of all the qtn.
createModel <- function(pheno, geno, traitName) {

  # check
  if (any(!pheno$ind %in% rownames(geno))) {
    # missing genotype
    warning("Missing genotype for some individuals")
  }

  # filter individuals
  geno <- geno[pheno$ind,]

  # check
  if (nrow(pheno) != nrow(geno)) {
    stop("nrow(pheno) != nrow(geno)")
  }

  ## rrBLUP:mixed.solve is too long (~20min for 35k markers)
  # mod <- mixed.solve(y = pheno[, trait$name],
  #                    Z = geno,
  #                    K = diag(1, ncol(geno)),
  #                    method = "REML")
  # names(mod$u) <- colnames(geno)
  # mod

  mod2 <- cv.glmnet(x = geno,
                    y = pheno[, traitName],
                    alpha = 0, # ridge regression
                    nfolds = 10)
  mod2
  beta <- coef(mod2, s = "lambda.min")

  list(inter = beta[1,],
       qtnEff = beta[-1,])
}

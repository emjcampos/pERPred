#' @title R2_pERP
#' @description This function will regress the true pERPs on the estimated pERPs.
#'
#' @param true_pERPs The dataframe containing the true pERPs.
#' @param pERPs The dataframe containing the estimated pERPs.
#'
#' @return a dataframe containing the number of pERPs and the associated R2_pERP value
#' @export
R2_pERP <- function(true_pERPs, pERPs) {
  datapoints <- length(true_pERPs[,1])
  num_pERPs <- ncol(pERPs)

  X <- pERPs
  svdICA <- svd(X)
  ProjICA <- svdICA$u %*% t(svdICA$u)

  residual <- true_pERPs - ProjICA %*% true_pERPs

  residual.var <- norm(residual, type = "F")^2
  original.var <- norm(true_pERPs, type = "F")^2
  rsquared <- 1 - residual.var/original.var
  data.frame("num_pERPs" = num_pERPs,
             "R2" = rsquared)
}

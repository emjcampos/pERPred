# code adapted from factoextra::get_eig
# factoextra not used due to ggplot2 dependency
get_eigenvalue = function (X)
{
  if (inherits(X, c("PCA", "CA", "MCA", "FAMD", "MFA", "HMFA",
                    "sPCA", "sCA", "sMCA", "sMFA", "sHMFA")))
    eig <- X$eig
  else {
    if (inherits(X, "prcomp") | inherits(X, "princomp"))
      eig <- (X$sdev)^2
    else if (inherits(X, c("pca", "coa", "acm")) & inherits(X,
                                                            "dudi"))
      eig <- X$eig
    else if (inherits(X, "ca"))
      eig <- X$sv^2
    else if (inherits(X, "mjca"))
      eig <- X$inertia.e
    else if (inherits(X, "correspondence"))
      eig <- X$cor^2
    else if (inherits(X, "expoOutput"))
      eig <- X$ExPosition.Data$eigs
    else stop("An object of class : ", class(X), " can't be handled by the function get_eigenvalue()")
    variance <- eig * 100/sum(eig)
    cumvar <- cumsum(variance)
    eig <- data.frame(eigenvalue = eig, variance = variance,
                      cumvariance = cumvar)
  }
  colnames(eig) <- c("eigenvalue", "variance.percent", "cumulative.variance.percent")
  rownames(eig) <- paste0("Dim.", 1:nrow(eig))
  eig
}

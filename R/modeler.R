#' @title Modeler
#' @description This function will regress the individual record on the estimated pERPs, we will need to demean the record before regression.
#' @param df A dataframe containing a single record from the observed data.
#' @param pERPs The estimated pERPs from the pERP-RED algorithm.
#'
#' @return The linear model from regressing the observed record on the estimated pERPs.
#' @export

modeler <- function(df, pERPs) {
  Signal <- NULL
  rm(list = "Signal")
  dat <- cbind(Signal = df$Signal, pERPs) %>%
    mutate(demean = scale(Signal, center = TRUE, scale = FALSE))
  if(sum(is.na(dat$Signal)) > 0){
    model <- NA
  } else {
    model <- lm(demean ~ . - 1 - Signal, data = dat)
  }
}

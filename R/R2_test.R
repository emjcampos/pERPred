# this pair of functions will also regress one record from the observed data on
# the estimated components, but it is faster!

#' @title sum_of_squares_calc
#' @description This function will regress one record from the observed data on the estimated pERPs at a time and calculate the variance of the residuals and the variance of the observed signal to be used in the global R2 calculation.
#'
#' @param df The observed data, typically the test set.
#' @param pERPs The estimated pERPs from the pERP-RED algorithm.
#'
#' @return  a dataframe containing the number of pERPs and the associated R2_test value
#'
#' @importFrom stats lm var
sum_of_squares_calc <- function(df, pERPs) {
  dat <- cbind(Signal = df$Signal, pERPs)
  if(sum(is.na(dat$Signal)) > 0){
    data.frame("RSS_g" = 0, "TSS_g" = 0)
  } else {
    model <- lm(Signal ~ .-1, data = dat)
    data.frame("RSS_g" = var(model$residuals), "TSS_g" = var(dat$Signal))
  }
}

#' @title R2_test
#' @description This function will regress one record from the observed data on the estimated pERPs at a time and calculate a global R2 value.
#'
#' @param df The observed data, typically the test set.
#' @param pERPs The estimated pERPs from the pERP-RED algorithm.
#'
#' @return  a dataframe containing the number of pERPs and the associated R2_test value
#' @export
R2_test <- function(df, pERPs) {
  Subject   <- NULL
  Time      <- NULL
  Task      <- NULL
  Signal    <- NULL
  Electrode <- NULL
  RSS_g     <- NULL
  TSS_g     <- NULL
  rm(list = c("Subject", "Time", "Task", "Signal", "Electrode", "RSS_g", "TSS_g"))

  rows    <- nrow(pERPs)
  columns <- ncol(pERPs)
  df %>%
    gather(Electrode, Signal, -c(Task, Subject, Time)) %>%
    group_by(Task, Subject, Electrode) %>%
    do(sum_of_squares_calc(., pERPs)) %>%
    ungroup() %>%
    summarise(R2 = 1 - sum(RSS_g) / sum(TSS_g)) %>%
    mutate(pERPs = columns)
}

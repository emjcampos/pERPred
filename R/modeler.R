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



#' @title pERP_scorer
#' @description This function calculates all of the individual score values.
#'
#' @param df A dataframe containing all of the observed ERP records.
#' @param pERPs The estimated pERPs from the pERP-RED algorithm.
#'
#' @return The dataframe containing each individual score for each of the records.
#' @export
#'
#' @importFrom broom tidy
#' @importFrom tidyr spread gather unite

pERP_scorer <- function(df, pERPs) {
  Electrode <- NULL
  Signal <- NULL
  Subject <- NULL
  Time <- NULL
  Task <- NULL
  model <- NULL
  term <- NULL
  rm(list = c("Electrode", "Signal", "Subject", "Time", "Task", "model", "term"))

  df %>%
    gather(Electrode, Signal, -c(Task, Subject, Time)) %>%
    group_by(Task, Subject, Electrode) %>%
    do(model = modeler(., pERPs)) %>%
    tidy(model) %>%
    mutate(term = gsub("`", "", term)) %>%
    ungroup()
}

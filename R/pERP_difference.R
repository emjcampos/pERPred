
#' @title Calculate t-tests for pERPs
#' @description This function will take in two comparison groups.
#'
#' @param scores The dataframe with all of the individual weights.
#' @param electrode A string of the electrode for testing.
#' @param task1 A string of the first task for comparison.
#' @param group1 A string of the first treatment/diagnosis group to test. If "Overall", all of the individuals will be averaged.
#' @param task2 A string of the second task for comparison, set to the same value as task1 by default.
#' @param group2 A string of the second treatment/diagnosis group to test. If "Overall", all of the individuals will be averaged.
#'
#' @return A dataframe of t-test results
#' @export
#' @importFrom stats pt
#'

pERP_difference <- function(scores, electrode, task1, group1 = "Overall",
                            task2 = task1, group2 = group1) {

  term <- NULL
  Average <- NULL
  APSD <- NULL
  SE <- NULL
  Average1 <- NULL
  Average2 <- NULL
  SE1 <- NULL
  SE2 <- NULL
  n1 <- NULL
  n2 <- NULL

  group_scores <- function(task, group, electrode) {

    Task <- NULL
    Electrode <- NULL
    term <- NULL
    estimate <- NULL
    Average <- NULL
    SE <- NULL
    group_member <- NULL

    if (group == "Overall") {
      scores %>%
        filter(grepl(task, Task),
               Electrode == electrode) %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average/SE,
                  n = n()) %>%
        mutate("Group" = task) %>%
        ungroup()
    } else {
      scores %>%
        filter(grepl(task, Task),
               Electrode == electrode,
               group_member == group) %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average/SE,
                  n = n()) %>%
        mutate("Group" = task) %>%
        ungroup()
    }
  }

  group_1_all <- group_scores(task1, group1, electrode)
  group_2_all <- group_scores(task2, group2, electrode)

  cbind(
    select(group_1_all, term, Average, APSD, SE, t, n),
    select(group_2_all, Average, APSD, SE, t, n)
  ) %>%
    setNames(c("pERP", "Average1", "APSD1", "SE1", "t1", "n1",
               "Average2", "APSD2", "SE2", "t2", "n2")) %>%
    mutate("Mean Difference" = Average1 - Average2,
           "SE Difference" = sqrt(SE1^2 + SE2^2),
           t = (Average1 - Average2)/sqrt(SE1^2 + SE2^2),
           significant = ifelse(2*pt(abs(t),
                                     n1 + n2 - 2,
                                     lower.tail = FALSE) < 0.054,
                                "*",
                                " ")) %>%
    setNames(c("pERP",
               paste(task1, group1, "Mean"), paste(task1, group1, "APSD"),
               paste(task1, group1, "SE"), paste(task1, group1, "t"),
               paste(task1, group1, "N"),
               paste(task2, group2, "Mean"), paste(task2, group2, "APSD"),
               paste(task2, group2, "SE"),
               paste(task2, group2, "t"), paste(task2, group2, "N"),
               "Mean Difference", "SE Difference", "t", "signif"))
}

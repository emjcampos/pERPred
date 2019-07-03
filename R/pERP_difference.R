
#' @title Calculate t-tests for pERPs
#' @description This function will take in two comparison groups.
#'
#' @param scores The dataframe with all of the individual weights.
#' @param electrode A string of the electrode for testing.
#' @param task1 A string of the first task for comparison.
#' @param group1 A string of the first treatment/diagnosis group to test. If "Overall", all of the individuals will be averaged.
#' @param task2 A string of the second task for comparison, set to the same value as task1 by default.
#' @param group2 A string of the second treatment/diagnosis group to test. If "Overall", all of the individuals will be averaged.
#' @param double_difference A logical to contrast task1 - task2. Use in cases such as correct - incorrect.
#'
#' @return A dataframe of t-test results
#' @export
#' @importFrom stats pt
#'

pERP_difference <- function(scores, electrode, task1, group1 = "Overall", task2 = task1, group2 = group1, double_difference = FALSE) {

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
  rm(list = c("term", "Average", "APSD", "SE", "Average1", "Average2", "SE1",
              "SE2", "n1", "n2"))

  group_scores <- function(task, group, electrode) {

    Task <- NULL
    Electrode <- NULL
    term <- NULL
    estimate <- NULL
    Average <- NULL
    SE <- NULL
    group_member <- NULL
    rm(list = c("Task", "Electrode", "term", "estimate", "Average", "SE",
                "group_member"))

    if (any(group == "Overall")) {
      scores %>%
        filter(Task == task,
               Electrode == electrode) %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average/SE,
                  n = n()) %>%
        ungroup()
    } else {
      scores %>%
        filter(Task == task,
               Electrode == electrode,
               group_member %in% group) %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average/SE,
                  n = n()) %>%
        ungroup()
    }
  }

  group_scores_double_diff <- function(tasks, group, electrode) {

    Task <- NULL
    Subject <- NULL
    Electrode <- NULL
    term <- NULL
    estimate <- NULL
    Average <- NULL
    SE <- NULL
    group_member <- NULL
    rm(list = c("Task", "Electrode", "term", "estimate", "Average", "SE",
                "group_member"))

    if (any(group == "Overall")) {
      scores %>%
        filter(Task %in% tasks,
               Electrode == electrode) %>%
        select(Task, Subject, group_member, Electrode, term, estimate) %>%
        spread(Task, estimate) %>%
        mutate(Difference = !!sym(tasks[1]) - !!sym(tasks[2])) %>%
        gather(Task, estimate,
               -c(Subject, group_member, Electrode, term)) %>%
        filter(Task == "Difference") %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average/SE,
                  n = n()) %>%
        ungroup()
    } else {
      scores %>%
        filter(Task %in% tasks,
               Electrode == electrode,
               group_member %in% group) %>%
        select(Task, Subject, group_member, Electrode, term, estimate) %>%
        spread(Task, estimate) %>%
        mutate(Difference = !!sym(tasks[1]) - !!sym(tasks[2])) %>%
        gather(Task, estimate,
               -c(Subject, group_member, Electrode, term)) %>%
        filter(Task == "Difference") %>%
        group_by(Electrode, term) %>%
        summarise(Average = mean(estimate),
                  APSD = sd(estimate),
                  SE = sd(estimate)/sqrt(n()),
                  t = Average / SE,
                  n = n()) %>%
        ungroup()
    }
  }

  if(double_difference) {
    group_1_all <- group_scores_double_diff(c(task1, task2), group1, electrode)
    group_2_all <- group_scores_double_diff(c(task1, task2), group2, electrode)

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
                 paste(task1, "-", task2, paste(group1, collapse = ', '), "Mean"),
                 paste(task1, "-", task2, paste(group1, collapse = ', '), "APSD"),
                 paste(task1, "-", task2, paste(group1, collapse = ', '), "SE"),
                 paste(task1, "-", task2, paste(group1, collapse = ', '), "t"),
                 paste(task1, "-", task2, paste(group1, collapse = ', '), "N"),
                 paste(task1, "-", task2, paste(group2, collapse = ', '), "Mean"),
                 paste(task1, "-", task2, paste(group2, collapse = ', '), "APSD"),
                 paste(task1, "-", task2, paste(group2, collapse = ', '), "SE"),
                 paste(task1, "-", task2, paste(group2, collapse = ', '), "t"),
                 paste(task1, "-", task2, paste(group2, collapse = ', '), "N"),
                 "Mean Difference", "SE Difference", "t", "signif"))
  } else {
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
                 paste(task1, paste(group1, collapse = ', '), "Mean"),
                 paste(task1, paste(group1, collapse = ', '), "APSD"),
                 paste(task1, paste(group1, collapse = ', '), "SE"),
                 paste(task1, paste(group1, collapse = ', '), "t"),
                 paste(task1, paste(group1, collapse = ', '), "N"),
                 paste(task2, paste(group2, collapse = ', '), "Mean"),
                 paste(task2, paste(group2, collapse = ', '), "APSD"),
                 paste(task2, paste(group2, collapse = ', '), "SE"),
                 paste(task2, paste(group2, collapse = ', '), "t"),
                 paste(task2, paste(group2, collapse = ', '), "N"),
                 "Mean Difference", "SE Difference", "t", "signif"))
  }
}

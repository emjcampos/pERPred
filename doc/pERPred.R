## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.path = "figures/",
  fig.align = "center", 
  fig.pos = "hold",
  fig.width = 6, 
  fig.height = 5,
  message = FALSE
)
library(here)
library(tidyverse)
library(latex2exp)
theme_set(theme_bw())

load(here("vignettes/data/true_pERPs.rda"))
load(here("vignettes/data/simulation_data.rda"))
load(here("vignettes/data/estimation.rda"))

## ----true_pERPs, echo = FALSE--------------------------------------------
true_pERPs %>%
  as.data.frame() %>%
  rename(
    "True pERP 1" = V1,
    "True pERP 2" = V2,
    "True pERP 3" = V3,
    "True pERP 4" = V4,
    "True pERP 5" = V5
  ) %>%
  mutate(Time = unique(simulated_data$Time)) %>%
  gather(Component, Signal, -Time) %>%
  mutate(Component = factor(Component, levels = c(
    "True pERP 1",
    "True pERP 2",
    "True pERP 3",
    "True pERP 4",
    "True pERP 5"
  ))) %>% 
  ggplot(aes(x = Time, y = Signal)) +
  geom_line() +
  facet_wrap( ~ Component, ncol = 2)

## ----head of simulated data----------------------------------------------
head(simulated_data[, 1:7])

## ----load pERPred package------------------------------------------------
library(pERPred)

## ----estimate pERPs, eval = FALSE----------------------------------------
#  subject_list <- unique(simulated_data$Subject)
#  electrode_list <- names(simulated_data)[-c(1:3)]
#  task_list <- unique(simulated_data$Task)
#  train_subjects <-
#    sort(sample(subject_list, round(2 * length(subject_list) / 3)))
#  
#  pERPs <- map(3:7,
#               ~ pERPred(simulated_data[simulated_data$Subject %in%
#                                          train_subjects, ],
#                         num_pERPs = .x,
#                         percent_variation = 80))
#  
#  # library(future)
#  # plan(multiprocess)
#  # pERPs <- future_map(3:7,
#  #                     ~ pERPred(simulated_data[simulated_data$Subject %in%
#  #                                                train_subjects, ],
#  #                               num_pERPs = .x,
#  #                               percent_variation = 80),
#  #                     .progress = TRUE)

## ----r2 test, eval = FALSE-----------------------------------------------
#  map_dfr(pERPs, ~R2_test(simulated_data, .x))

## ----load r2 test, echo = FALSE------------------------------------------
R2

## ----re-estimate pERPs, eval = FALSE-------------------------------------
#  pERPs5 <- pERPred(simulated_data,
#                   num_pERPs = 5,
#                   percent_variation = 80)

## ----plot_pERPs----------------------------------------------------------
pERPs5 %>%
  mutate(Time = unique(simulated_data$Time)) %>%
  gather(pERP, Amplitude, -Time) %>%
  ggplot() +
  geom_line(aes(x = Time, y = Amplitude)) +
  facet_wrap(~ pERP, ncol = 2)

## ----individual scores, eval = FALSE-------------------------------------
#  individual_scores <- pERP_scorer(simulated_data, pERPs5)

## ----head individual scores----------------------------------------------
head(individual_scores)

## ----individual_reconstruction-------------------------------------------
record <- simulated_data %>% 
  filter(Subject == "Subject_001", 
         Task == "Task_1") %>% 
  select(Task, Time, "CZ")
scores <- individual_scores %>% 
  filter(Subject == "Subject_001", 
         Electrode == "CZ", 
         Task == "Task_1") %>% 
  pull(estimate)
projected <- data.frame("Projected" = as.matrix(pERPs5) %*% 
                          as.matrix(scores)) %>% 
  as.data.frame() %>% 
  mutate(Time = record$Time) %>% 
  mutate(Projected = Projected + mean(record$CZ))

ggplot() + 
  geom_line(data = record, 
            aes(x = Time, y = CZ, color = "Observed")) + 
  geom_line(data = projected, 
            aes(x = Time, y = Projected, color = "Reconstructed")) +
  theme(legend.title = element_blank()) + 
  labs(y = TeX("$\\mu V$"))

## ----headmaps------------------------------------------------------------
average_scores <- individual_scores %>% 
  group_by(Task, Electrode, term) %>% 
  summarise(average = mean(estimate, na.rm = TRUE)) %>% 
  split(.$Task)

# keep scale the same across all tasks
coefficient_headmap("Task_1", average_scores)

# keep scale the same across only tasks 1 and 2
coefficient_headmap("Task_1", average_scores[c("Task_1", "Task_2")])

## ----pERP difference-----------------------------------------------------
pERP_difference(scores = individual_scores,
                electrode = "CZ",
                task1 = "Task_1",
                task2 = "Task_2")

## ----pERP difference with groups-----------------------------------------
group <- data.frame(
  Subject = distinct(simulated_data, Subject),
  group_member = sample(
    c("Control", "Treatment"),
    100,
    replace = TRUE,
    prob = c(.5, .5)
  )
)

individual_scores_groups <- full_join(individual_scores, group, by = "Subject")

pERP_difference(scores = individual_scores_groups,
                electrode = "CZ",
                task1 = "Task_1",
                group1 = "Control", 
                group2 = "Treatment")


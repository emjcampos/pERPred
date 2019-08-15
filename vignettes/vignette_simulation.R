library(here)
library(KRLS2)
library(ica)
library(LaplacesDemon)
library(magrittr)
library(purrr)
library(tidyverse)
source('~/Box Sync/CanonicalERP/code/00_generate-canonical-components.R')
source('~/Box Sync/CanonicalERP/code/01_simulate-data.R')
theme_set(theme_bw())

set.seed(1234)

num_tasks      <- 9
num_electrodes <- 40
num_subjects   <- 100
num_components <- 5
c_variance     <- .25
datapoints     <- 384
start_time     <- 0
end_time       <- 1
subject_cov_multipliers <- seq(.1, .5, length.out = num_components)

noiselevel     <- sqrt(0.5) # 2/3 signal, 1/3 noise
simulated_info <- simulate_data()
simulated_data <- simulated_info$full_data %>%
  setNames(c("Task", "Subject", "Time", "AF3", "AF4", "AFZ", "C3", "C4", "CPZ",
             "CZ", "F10", "F3", "F4", "F7", "F8", "F9", "FCZ", "FP1", "FP2",
             "FPZ", "FT10", "FT7", "FT8", "FT9", "FZ", "IZ", "O1", "O2", "OZ",
             "P10", "P3", "P4", "P7", "P8", "P9", "POZ", "PZ", "T7", "T8",
             "TP10", "TP7", "TP8", "TP9")) %>%
  mutate(Subject = paste0("Subject_",
                          str_pad(Subject, 3, side = "left", pad = "0")),
         Task = paste0("Task_", Task))
true_pERPs <- simulated_info$independent_components

subject_list <- unique(simulated_data$Subject)
electrode_list <- names(simulated_data)[-c(1:3)]
task_list <- unique(simulated_data$Task)
train_subjects <-
  sort(sample(subject_list, round(2 * length(subject_list) / 3)))

pERPs <- map(3:7,
             ~ pERPred(simulated_data[simulated_data$Subject %in%
                                        train_subjects, ],
                       num_pERPs = .x,
                       percent_variation = 80))

R2 <- map_dfr(pERPs, ~R2_test(simulated_data, .x))

pERPs5 <- pERPred(simulated_data,
                  num_pERPs = 5,
                  percent_variation = 80)

individual_scores <- pERP_scorer(simulated_data, pERPs5)


simulated_data <- mutate_if(simulated_data,
                            is.numeric,
                            round, 7)

save(true_pERPs,
     file = here("vignettes/data/true_pERPs.rda"))
save(simulated_data,
     file = here("vignettes/data/simulation_data.rda"))
save(individual_scores, pERPs5, R2,
     file = here("vignettes/data/estimation.rda"))

# save the electrode locations as internal data for the package
electrodeLocs <- read_csv(
  "~/Box Sync/CanonicalERP/data/electrodeLocs81.csv",
  col_types = cols()
) %>%
  rename(Electrode = electrode)
usethis::use_data(electrodeLocs, internal = TRUE, overwrite = TRUE)

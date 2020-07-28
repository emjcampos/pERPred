## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#
#' @title Principle ERP Reduction
#' @description This function will contain the estimation algorithm of the principle ERPs (pERPs).
#'
#' @param df The dataframe with a column for Task, Subject, and Time. The remaining columns are the Electrodes.
#' @param num_pERPs The number of pERPs to estimate.
#' @param percent_variation_electrode The percent variation to use in the electrode PCA step.
#' @param percent_variation_subject The percent variation to use in the subject PCA step.
#'
#' @return pERPs The prinicple ERPs (pERPs) are the bases functions estimated by the pERP-RED algorithm.
#' @export
#'
#' @importFrom stats sd princomp complete.cases setNames
#' @import dplyr
#' @import fastICA
#' @importFrom tidyr spread gather unite
#' @importFrom purrr map_chr map_dfr
#' @importFrom tibble rownames_to_column
#' @importFrom factoextra get_eigenvalue
#' @importFrom glue glue

pERPred <- function(df, num_pERPs = 20, percent_variation_electrode = 80, percent_variation_subject = 80) {
  # to avoid issues with non-standard evaluation in tidyeval, set "global
  # variables" to NULL and remove them. this won't cause an issue with the rest
  # of the code.
  Subject        <- NULL
  Time           <- NULL
  Task           <- NULL
  Region         <- NULL
  Signal         <- NULL
  Subject_Region <- NULL
  Component      <- NULL
  SubjectRegion  <- NULL
  V1             <- NULL
  rm(list = c("Subject", "Time", "Task", "Region", "Signal", "Subject_Region",
              "Component", "SubjectRegion", "V1"))

  # Setup ------------------------------------------------------------------

  datapoints     <- length(unique(df$Time))
  start_time     <- min(df$Time)
  end_time       <- max(df$Time)
  subject_list   <- unique(df$Subject)
  num_subjects   <- length(subject_list)
  task_list      <- unique(df$Task)
  num_tasks      <- length(task_list)
  electrode_list <- names(df)[!names(df) %in% c("Task", "Subject", "Time")]
  num_electrodes <- length(electrode_list)

  # Normalize each record --------------------------------------------------

  scale_this <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }

  scaled_data <- df %>%
    group_by(Subject, Task) %>%
    mutate_at(vars(-Task, -Subject, -Time), scale_this) %>%
    ungroup()


  # Reduce electrode dimension ---------------------------------------------

  # set up with empty lists
  electrode_pcas             <- vector("list")
  electrode_loads            <- vector("list")
  electrode_var_tables       <- vector("list")
  electrode_num_comps_chosen <- vector("integer")

  for (person in subject_list) {

    try(if (num_electrodes > (datapoints * num_tasks)) {
      stop(glue::glue(
        "Uh oh! PCA requires that the number of subject-regions is smaller than the number of tasks x timepoints! You have {num_electrodes} subject-regions and {datapoints * num_tasks} tasks x timepoints, whomp whomp."
      ),
      call. = FALSE)
    })

    # store the pca results
    electrode_pcas[[person]] <- scaled_data %>%
      # fix subject
      filter(Subject == person) %>%
      # remove Task Subject and Time variables
      select(-c(Task, Subject, Time)) %>%
      # remove missing electrodes
      select_if( ~ !any(is.na(.))) %>%
      princomp()

    # store the loadings from the pca results
    electrode_loads[[person]] <-
      electrode_pcas[[person]]$loadings %>%
      as.data.frame.matrix()

    # store the eigenvalues
    electrode_var_tables[[person]] <-
      get_eigenvalue(electrode_pcas[[person]])

    # store the number of components chosen based on the percent of variation explained
    electrode_num_comps_chosen[[person]] <-
      min(
        which(
          electrode_var_tables[[person]]$cumulative.variance.percent >= percent_variation_electrode
        )
      )

  }

  reduce_electrodes <- function(thissubject) {
    # this function will calculate the rotated data based on the principal
    # electrodes (these are the scores)
    subject_data <- scaled_data %>%
      filter(Subject == thissubject) %>%
      select_if( ~ !any(is.na(.))) %>%
      select(-c(Task, Subject, Time))

    principal_electrodes <- vector(
      "list",
      electrode_num_comps_chosen[thissubject]
    )
    for (i in 1:electrode_num_comps_chosen[thissubject]) {
      principal_electrodes[[i]] <- as.matrix(subject_data) %*%
        as.matrix(electrode_loads[[thissubject]][, i])
    }

    names <- map_chr(1:electrode_num_comps_chosen[thissubject],
                     ~ ifelse(nchar(.) < 2, paste0("PE0", .), paste0("PE", .)))
    dat <- as.data.frame(do.call(cbind, principal_electrodes)) %>%
      mutate(
        Task = rep(task_list, each = datapoints),
        Subject = rep(thissubject, datapoints * num_tasks),
        Time = rep(seq(start_time, end_time, length.out = datapoints),
                   times = num_tasks)
      ) %>%
      select(Task, Subject, Time, everything())
    colnames(dat) <- c("Task", "Subject", "Time", names)
    dat

  }

  electrodes_reduced <-
    map_dfr(subject_list, reduce_electrodes) %>%
    gather(key = "Region", value = "Signal",-c(Task, Time, Subject)) %>%
    select(Task, Region, Subject, Time, Signal) %>%
    arrange(Task, Region, Subject, Time) %>%
    unite(Subject_Region, c("Subject", "Region")) %>%
    filter(complete.cases(.)) %>%
    spread(Subject_Region, Signal)


  # Re-normalize -----------------------------------------------------------

  normalize_subject_region <- function(df, task) {
    df %>%
      filter(Task == task) %>%
      select(-c(Task, Time)) %>%
      scale() %>%
      as.data.frame()
  }

  scaled_subject_regions <- map_dfr(task_list,
                                    ~ normalize_subject_region(electrodes_reduced, .x))


  # Reduce Subject-Regions -------------------------------------------------
  try(if (ncol(scaled_subject_regions) > nrow(scaled_subject_regions)) {
    stop(glue::glue(
      "Uh oh! PCA requires that the number of subject-regions is smaller than the number of tasks x timepoints! You have {ncol(scaled_subject_regions)} subject-regions and {nrow(scaled_subject_regions)} tasks x timepoints, whomp whomp."
    ),
    call. = FALSE)
  })

  subject_region_pca <- scaled_subject_regions %>%
    princomp()

  subject_region_vartable <- get_eigenvalue(subject_region_pca)
  num_subject_regions_chosen <- min(
    which(
      subject_region_vartable$cumulative.variance.percent >= percent_variation_subject
    )
  )

  principal_subject_region <- as.matrix(scaled_subject_regions) %*%
    as.matrix(subject_region_pca$loadings[, 1:num_subject_regions_chosen]) %>%
    as.data.frame() %>%
    rename_all( ~ map_chr(
      1:num_subject_regions_chosen,
      ~ ifelse(nchar(.) < 2,
               paste0("PSR0", .),
               paste0("PSR", .))
    ))

  reduced_data <- principal_subject_region %>%
    as.data.frame() %>%
    mutate(
      Task = rep(task_list, each = datapoints),
      Time = rep(seq(start_time, end_time, length.out = datapoints),
                 times = num_tasks)
    ) %>%
    select(Task, Time, everything()) %>%
    gather("SubjectRegion", "Signal",-c(Task, Time)) %>%
    unite(Component, Task, SubjectRegion) %>%
    spread(Component, Signal) %>%
    select(-Time)

  ica_results <- as.data.frame(fastICA(
    reduced_data,
    n.comp = num_pERPs,
    method = "C",
    maxit = 100
  )$S) %>%
    setNames(map_chr(1:num_pERPs,
                     ~ ifelse(nchar(.x) < 2,
                              paste0("pERP 0", .),
                              paste0("pERP ", .))))

  # ordering the estimated components by the largest peak
  est_comp_order <- ica_results %>%
    summarise_all(funs(which.max(abs(.)))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Component") %>%
    arrange(V1)

  ica_results_ordered <- ica_results %>%
    as.data.frame() %>%
    select(est_comp_order$Component) %>%
    setNames(map_chr(1:num_pERPs,
                     ~ ifelse(
                       nchar(.x) < 2,
                       paste0("pERP 0", .x),
                       paste0("pERP ", .x)
                     )))

  ica_results_ordered
}

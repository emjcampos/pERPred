#' @title Shiny App pERP Explorer
#' @description In our paper, we often referenced a Shiny app in order to explore all of the results from the pERP-space analysis. Now you can, too! This function will take all of your data and results and launch a shiny application on your machine.
#'
#' @param data A dataframe containing all of the ERP data used to calculate the pERPs. One column for Task, Subject, and Time, and the rest are electrodes with the ERP data in the rows.
#' @param groups A dataframe containing Subject, which should have the same form as the Subject column in the data, and Group, a group membership indicator.
#' @param pERPs A dataframe containing the pERPs. This is the resulting dataframe from running `pERPred`.
#' @param scores A dataframe containing each of the individual scores from regressing the observed data on the pERPs. This is the resulting dataframe from running `pERP_scorer`.
#'
#' @return shiny application
#' @export
#'
#' @import shiny
#' @import shinythemes
#' @import ggplot2
#' @importFrom kableExtra kable_styling add_header_above
#' @importFrom knitr kable
#' @importFrom gridExtra grid.arrange
#' @importFrom stringr str_remove
#' @importFrom purrr map
#' @importFrom plyr round_any

shiny_pERPs <- function(data, groups, pERPs, scores) {
  group <- NULL
  Component <- NULL
  Signal <- NULL
  Time <- NULL
  Subject <- NULL
  Task <- NULL
  Electrode <- NULL
  estimate <- NULL
  V1 <- NULL
  Projected <- NULL
  term <- NULL
  rm(list = c("group", "Component", "Signal", "Time", "Subject", "Task",
              "Electrode", "estimate", "V1", "Projected", "term"))

  theme_set(theme_bw())

  shinyApp(

    ui <- fluidPage(
      withMathJax(),

      # the application title
      titlePanel("pERP Exploration"),

      tabsetPanel(
        # panel: pERPs ----
        tabPanel(
          "pERPs", br(),
          plotOutput("pERPs_plot", height = "700px")
        ),

        # panel: individuals ----
        tabPanel(
          "Individual Reconstruction", br(),
          sidebarPanel(
            p("Each of the observed ERPs can be reconstructed using the pERPs. The observed ERP is regressed on the pERPs to obtain a set of scores and using those scores, the pERPs are projected onto the observed ERP. This typically results in a smoother ERP."),
            hr(),
            selectInput("subject", "Subject:",
                        choices = sort(unique(groups$Subject))),
            selectInput("task", "Task:",
                        choices = sort(unique(data$Task))),
            selectInput("electrode", "Electrode:",
                        choices = sort(colnames(data[,-c(1:3)])))
          ),
          mainPanel(
            plotOutput("erp")
          )
        ),

        # panel: headmaps ----
        tabPanel(
          "Headmaps", br(),
          wellPanel(
            fluidRow(
              column(5,
                     p("Here, the scores from each task/condition were averaged across all subjects and plotted for each electrode to visualize where activity occurs on the scalp."),
                     helpText("Select 'Standardized' if you would like to see the average scores divided by their standard devations. In the 'Difference' headmaps, this would indicate a significant difference.")
              ),
              column(7,
                     selectInput("headmap_task", "Task:",
                                 choices = unique(data$Task)),
                     selectInput("compare_task", "Comparison Task:",
                                 choices = c("None", unique(data$Task))),
                     checkboxInput("standard", "Standardized")
              )
            )
          ),
          plotOutput("headmap")
        ),
        # panel: groups ----
        tabPanel(
          "Condition Differences", br(),
          sidebarLayout(
            sidebarPanel(
              width = 4,
              h3("Differences"),
              p("Here, you can compare groups. Start by selecting an electrode of interest, then choosing the task/group combination that you would like to compare. The plot on the top is the projected signal using all of the pERPs but the plot on the bottom is only the pERPs you select. Use this to explore how combinations of pERPs change how the ERP looks. The table computes the t-test for the means."),
              helpText("Select 'Double Difference' to contrast Task 1 - Task 2. Use this to investigate condition differences such as correct - incorrect for different treatment/diagnosis groups."),
              selectInput(
                "diff_electrode", "Electrode:",
                choices = sort(colnames(data[, -which(names(data) %in%
                                                        c("Time","Subject",
                                                          "Task"))])),
                selected = "CZ"),
              h4("Group 1"),
              selectInput("diff_task1", "Task:",
                          choices = unique(data$Task)),
              uiOutput("condition1"),
              selectInput("diff_group1", "Group:",
                          choices = unique(groups$group_member)),
              h4("Group 2"),
              selectInput("diff_task2", "Task:",
                          choices = unique(data$Task)),
              uiOutput("condition2"),
              selectInput("diff_group2", "Group:",
                          choices = unique(groups$group_member)),
              checkboxInput("double_diff", "Double Difference"),
              list(
                tags$div(align = "left",
                         class = "multicol",
                         checkboxGroupInput(
                           "checkGroup",
                           label = "Select pERPs:",
                           choices = names(pERPs),
                           selected = "pERP 01",
                           inline = FALSE),
                         checkboxInput("bar", "All/none")
                )
              )
            ),
            mainPanel(plotOutput("perp_diff", height = "550px"),
                      tableOutput("score_table"))
          )
        )
      )
    ),

    server <- function(input, output, session) {
      pERP <- NULL
      Mean1 <- NULL
      Mean2 <- NULL
      Group <- NULL
      rm(list = c("pERP", "Mean1", "Mean2", "Group"))

      start_time <- min(data$Time)
      end_time <- max(data$Time)
      datapoints <- length(unique(data$Time))
      plot_breaks <- seq(round_any(start_time, 100, f = floor),
                         round_any(end_time, 100, f = ceiling),
                         by = 100)

      # output$pERPs_plot ----
      output$pERPs_plot <- renderPlot({
        pERPs %>%
          as.data.frame() %>%
          mutate(Time = seq(start_time, end_time, length.out = datapoints)) %>%
          gather(Component, Signal, -Time) %>%
          ggplot(aes(x = Time, y = Signal)) +
          geom_line() +
          facet_wrap(~ Component, ncol = 3) +
          scale_x_continuous(breaks = plot_breaks) +
          labs(title = "Estimated pERPs",
               x = "Time (ms)")
      })

      # output$erp ----
      output$erp = renderPlot({
        df <- data %>%
          filter(Subject == input$subject,
                 Task == input$task) %>%
          select(Task, Time, input$electrode) %>%
          rename(Signal = input$electrode)
        scores_sub <- scores %>%
          filter(Subject == str_remove(input$subject, '"'),
                 Task == input$task,
                 Electrode == input$electrode) %>%
          pull(estimate)
        if (length(scores_sub) == 0) {
          ggplot() +
            geom_text(aes(x = end_time/2, y = 7,
                          label = "Missing Electrode")) +
            geom_line(aes(x = end_time, y = 0, color = "Observed")) +
            geom_line(aes(x = start_time, y = 0, color = "Reconstructed")) +
            scale_x_continuous(breaks = plot_breaks) +
            ylim(-10, 10) +
            theme(legend.title = element_blank()) +
            labs(x = "Time (ms)",
                 y = "Signal")
        } else{
          projected <- as.matrix(pERPs) %*% as.matrix(scores_sub) %>%
            as.data.frame() %>%
            mutate(Time = df$Time) %>%
            rename(Projected = V1) %>%
            mutate(Projected = Projected + mean(df$Signal))
          ggplot() +
            geom_line(data = df,
                      aes(x = Time, y = Signal, color = "Observed")) +
            geom_line(data = projected,
                      aes(x = Time, y = Projected, color = "Reconstructed")) +
            scale_x_continuous(breaks = plot_breaks) +
            theme(legend.title = element_blank())
        }
      })
      # output$headmap ----
      output$headmap <- renderPlot({
        if(input$compare_task == "None") {
          if(input$standard) {
            average_scores <- scores %>%
              filter(Task == input$headmap_task) %>%
              group_by(Task, Electrode, term) %>%
              summarise(average = mean(estimate, na.rm = TRUE) /
                          (sd(estimate, na.rm = TRUE)/sqrt(n()))) %>%
              split(.$Task)
            coefficient_headmap(input$headmap_task, average_scores, input$standard)
          } else {
            average_scores <- scores %>%
              filter(Task == input$headmap_task) %>%
              group_by(Task, Electrode, term) %>%
              summarise(average = mean(estimate, na.rm = TRUE)) %>%
              split(.$Task)
            coefficient_headmap(input$headmap_task, average_scores)
          }
        } else {
          if(input$standard) {
            average_scores <- scores %>%
              filter(Task %in% c(input$headmap_task, input$compare_task)) %>%
              select(Task, Subject, Electrode, term, estimate) %>%
              spread(Task, estimate) %>%
              mutate(Difference = !!sym(input$headmap_task) -
                       !!sym(input$compare_task)) %>%
              gather(Task, estimate, -c(Subject, Electrode, term)) %>%
              group_by(Task, Electrode, term) %>%
              summarise(average = mean(estimate, na.rm = TRUE) /
                          (sd(estimate, na.rm = TRUE)/sqrt(n()))) %>%
              split(.$Task)

            plots <- map(
              c(input$headmap_task, input$compare_task),
              ~ coefficient_headmap(.x, average_scores, input$standard)
            )
            plots[[3]] <- coefficient_headmap("Difference", average_scores["Difference"], input$standard)
            grid.arrange(grobs = plots, ncol = 3)

          } else {
            average_scores <- scores %>%
              filter(Task %in% c(input$headmap_task, input$compare_task)) %>%
              select(Task, Subject, Electrode, term, estimate) %>%
              spread(Task, estimate) %>%
              mutate(Difference = !!sym(input$headmap_task) -
                       !!sym(input$compare_task)) %>%
              gather(Task, estimate, -c(Subject, Electrode, term)) %>%
              group_by(Task, Electrode, term) %>%
              summarise(average = mean(estimate, na.rm = TRUE)) %>%
              split(.$Task)

            plots <- map(
              c(input$headmap_task, input$compare_task),
              ~ coefficient_headmap(.x, average_scores)
            )
            plots[3] <- coefficient_headmap("Difference", average_scores["Difference"])
            grid.arrange(grobs = plots, ncol = 3)
          }

        }
      })
      # update checkGroup ----
      observe({
        updateCheckboxGroupInput(
          session, "checkGroup",
          choices = names(pERPs),
          selected = if (input$bar) names(pERPs)
        )
      })

      # output$perp_diff ----
      output$perp_diff <- renderPlot({
        compare_scores <- pERP_difference(
          scores,
          input$diff_electrode,
          input$diff_task1,
          input$diff_group1,
          input$diff_task2,
          input$diff_group2,
          input$double_diff
        ) %>%
          setNames(c("pERP", "Mean1", "APSD1", "SE1", "t1", "n1", "Mean2", "APSD2", "SE2", "t2", "n2", "difference", "SE", "t", "signif"))
        group_1_selected <- compare_scores %>%
          filter(pERP %in% input$checkGroup) %>%
          select(Mean1) %>%
          as.matrix()
        group_2_selected <- compare_scores %>%
          filter(pERP %in% input$checkGroup) %>%
          select(Mean2) %>%
          as.matrix()
        if(input$double_diff) {
          group_1_all <- data.frame(
            as.matrix(pERPs) %*%
              as.matrix(select(compare_scores, Mean1))
          ) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, "-", input$diff_task2,
                                 input$diff_group1)) %>%
            rename(Signal = Mean1)
          group_2_all <- data.frame(
            as.matrix(pERPs) %*%
              as.matrix(select(compare_scores, Mean2))
          ) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, "-", input$diff_task2,
                                 input$diff_group2)) %>%
            rename(Signal = Mean2)

          selected_pERPs <- select(pERPs, input$checkGroup) %>% as.matrix()

          group_1_selected <- data.frame(
            selected_pERPs %*% group_1_selected
          ) %>%
            rename(Signal = Mean1) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, "-", input$diff_task2,
                                 input$diff_group1))
          group_2_selected <- data.frame(
            selected_pERPs %*% group_2_selected
          ) %>%
            rename(Signal = Mean2) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, "-", input$diff_task2,
                                 input$diff_group2))
        } else {
          group_1_all <- data.frame(
            as.matrix(pERPs) %*%
              as.matrix(select(compare_scores, Mean1))
          ) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, input$diff_group1)) %>%
            rename(Signal = Mean1)
          group_2_all <- data.frame(
            as.matrix(pERPs) %*%
              as.matrix(select(compare_scores, Mean2))
          ) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task2, input$diff_group2)) %>%
            rename(Signal = Mean2)

          selected_pERPs <- select(pERPs, input$checkGroup) %>% as.matrix()

          group_1_selected <- data.frame(
            selected_pERPs %*% group_1_selected
          ) %>%
            rename(Signal = Mean1) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task1, input$diff_group1))
          group_2_selected <- data.frame(
            selected_pERPs %*% group_2_selected
          ) %>%
            rename(Signal = Mean2) %>%
            mutate(Time = unique(data$Time),
                   Group = paste(input$diff_task2, input$diff_group2))
        }

        limits <- c("min" = min(group_1_all$Signal,
                                group_2_all$Signal,
                                group_1_selected$Signal,
                                group_2_selected$Signal),
                    "max" = max(group_1_all$Signal,
                                group_2_all$Signal,
                                group_1_selected$Signal,
                                group_2_selected$Signal))

        p1 <- ggplot(mapping = aes(x = Time, y = Signal, color = Group)) +
          geom_line(data = group_1_all) +
          geom_line(data = group_2_all) +
          theme(legend.position = "bottom") +
          ylim(limits["min"], limits["max"]) +
          scale_x_continuous(breaks = plot_breaks) +
          labs(x = "Time",
               y = "Projected ERP",
               title = "Using all pERPs")
        p2 <- ggplot(mapping = aes(x = Time, y = Signal, color = Group)) +
          geom_line(data = group_1_selected) +
          geom_line(data = group_2_selected) +
          theme(legend.position = "bottom") +
          ylim(limits["min"], limits["max"]) +
          scale_x_continuous(breaks = plot_breaks) +
          labs(x = "Time",
               y = "Projected ERP",
               title = "Using selected pERPs")
        grid.arrange(grobs = list(p1, p2), ncol = 1)

      })

      # output$score_table ----
      output$score_table <- function() {
        pERP_difference(
          scores,
          input$diff_electrode,
          input$diff_task1,
          input$diff_group1,
          input$diff_task2,
          input$diff_group2,
          input$double_diff
        ) %>%
          kable(booktabs = TRUE, linesep = "", escape = TRUE,
                longtable = TRUE, digits = 2,
                col.names = c("pERP",
                              "Mean", "APSD", "SE", "t", "N",
                              "Mean", "APSD", "SE", "t", "N",
                              "Mean", "SE", "t", "Significant")) %>%
          kable_styling(full_width = FALSE,
                        latex_options = c("repeat_header")) %>%
          add_header_above(c(" " = 1, "Group 1" = 5, "Group 2" = 5,
                             "Contrast" = 4))
      }
    }
  )
}

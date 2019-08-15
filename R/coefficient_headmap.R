
#' @title Scalp distribution of the scores on each of the pERPs
#' @description This function will create a heatmap of scalp distribution of the scores on each of the pERPs.
#'
#' @param task The string for one of the tasks used to estimate the pERPs.
#' @param scores A dataframe of scores calculated for the task. In the case of multiple tasks to be compared, use a named list of dataframes where the names are the tasks.
#'
#' @return A ggplot of with the headmaps of the scores averaged within the task for each pERP
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#'

coefficient_headmap <- function(task, scores) {
  average <- NULL
  Electrode <- NULL
  x <- NULL
  y <- NULL
  rm(list = c("average", "Electrode", "x", "y"))

  theme_topo <- function(base_size = 12) {
    theme_bw(base_size = base_size) %+replace%
      theme(
        rect = element_blank(),
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()
      )
  }

  circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
    r = diameter / 2
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }

  headShape <- circleFun(c(0, 0), round(max(electrodeLocs$x)), npoints = 100)
  nose <- data.frame(x = c(-0.075, 0, .075), y = c(.495, .575, .495))

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))

  limits <- map_dfr(scores, ~.x) %>%
    ungroup () %>%
    summarise("max" = max(average),
              "min" = min(average))


  scores_long <- scores[[task]] %>%
    left_join(electrodeLocs, by = "Electrode")

  ggplot(headShape, aes(x, y)) +
    geom_path(size = 1.5) +
    geom_point(data = scores_long,
               aes(x, y, colour = average),
               size = 7) +
    geom_text(data = scores_long, aes(x, y, label = Electrode), size = 3) +
    scale_colour_gradientn(colours = jet.colors(10),
                           guide = "colourbar",
                           limits = c(limits$min, limits$max)) +
    geom_line(data = nose, aes(x, y, z = NULL), size = 1.5) +
    theme_topo() +
    facet_wrap(~ term, ncol = ceiling(length(unique(scores_long$term))/2)) +
    coord_equal() +
    labs(title = task)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefficient_headmap.R
\name{coefficient_headmap}
\alias{coefficient_headmap}
\title{Scalp distribution of the scores on each of the pERPs}
\usage{
coefficient_headmap(task, scores)
}
\arguments{
\item{task}{The string for one of the tasks used to estimate the pERPs.}

\item{scores}{A dataframe of scores calculated for the task. In the case of multiple tasks to be compared, use a named list of dataframes where the names are the tasks.}
}
\value{
A ggplot of with the headmaps of the scores averaged within the task for each pERP
}
\description{
This function will create a heatmap of scalp distribution of the scores on each of the pERPs.
}

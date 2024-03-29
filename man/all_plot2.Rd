% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_plot2.R
\name{all_plot2}
\alias{all_plot2}
\title{Plots all effect estimates against p values with each specific variable in the models}
\usage{
all_plot2(
  data,
  xlabels = c(0, 0.001, 0.01, 0.05, 0.2, 0.5, 1),
  xlim = c(0, 1),
  xlab = "P value",
  ylim = NULL,
  ylab = NULL,
  yscale_log = FALSE,
  title = NULL
)
}
\arguments{
\item{data}{\emph{Object} from \code{all_cox}, \code{all_glm}, \code{all_speedglm}, or \code{all_glm}, including all effect estimate values.}

\item{xlabels}{\emph{numeric vector} x-axis tick labels. Default is
\code{"c(0, 0.001, 0.01, 0.05, 0.2, 0.5, 1)"}}

\item{xlim}{\emph{vector} of 2 numeric values for x-axis limits. Default is \code{"c(0, 1)"}.}

\item{xlab}{\emph{Character} string for x-axis name. Default is \code{"P value"}.}

\item{ylim}{\emph{vector} of 2 numeric values for y-axis limits.}

\item{ylab}{\emph{Character} string for y-axis name. Default depends on original model types.}

\item{yscale_log}{\emph{TRUE or FALSE} re-scale y-axis to "log10". Default is \code{"FALSE"}.}

\item{title}{\emph{Character} title. Default is \code{"NULL"}.}
}
\value{
A \pkg{ggplot2} object: scatter plot
}
\description{
\code{all_plot2()}  generates a panel of scatter plots with effect estimates of all possible models
again p values. Each plot includes effect estimates from all models including a specific variable.
}
\examples{
vlist <- c("Age", "Sex", "Married", "BMI", "Income")
results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
all_plot2(results)
}

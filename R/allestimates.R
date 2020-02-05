#'
#' Effect estimates from models with all possible combinations of variables
#'
#'  To assess treatment effects in clinical trials and risk factors in bio-medical
#'  and epidemiological research, we use
#'  \emph{regression coefficients}, \emph{odds ratios} or \emph{hazard ratios} as
#'  \emph{effect estimates}. \code{allestimates} allows users to quickly obtain
#'  effect estimates from models with all possible combinations of a list of variables
#'  specified by users. \code{all_lm} for linear regression, \code{all_glm} for
#'  logistic regression, \code{all_speedglm}
#'  using \code{speedlm} as a faster alternative of \code{all_glm}, and
#'  \code{all_cox} for Cox Proportional Hazards Models. Users can further
#'  use those values in a returned list of results.
#'  \code{all_plot} draws scatter plots with all effect
#'  estimate values against p values, as \code{Stata confall} command
#'  (Wang Z (2007) <doi:10.1177/1536867X0700700203>).
#'  Those plots divide estimates into four categories:
#'
#'  \itemize{
#'    \item{positive and significant:} {left-top quarter}
#'    \item{negative and significant:} {left-bottom quarter}
#'    \item{positive and non-significant:} {right-top quarter}
#'    \item{negative and non-significant:} {right-bottom quarter}
#' }
#'
#'  \code{all_plot2} draws multiple plots. Each of those plots
#'  indicates whether a specific variable is included or
#'  not included in models.
#'  Those effect estimates help users better understand
#'  confounding effects, uncertainty of their estimates, as well as
#'  inappropriately including variables in the models. This is a tool for
#'  calculating and exploring effect estimates from all possible models.
#'  Interpretation of the results should be in the context of other
#'  analyses and biological knowledge.
#'
#' @docType package
#' @importFrom  stats AIC as.formula binomial glm lm na.omit nobs
#' @importFrom broom tidy
#' @import utils
#' @import ggplot2
#' @importFrom stringr str_count str_detect
#' @importFrom dplyr mutate filter
#' @import tidyr
#' @name allestimates
#' @examples
#'
#' ? all_speedglm
#' ? all_glm
#' ? all_cox
#' ? all_lm
#' ? all_plot
#' ? all_plot2
utils::globalVariables(c(
  "p_value", "variables", "p", "aic", "n",
  "estimate", "crude", "value"))
NULL

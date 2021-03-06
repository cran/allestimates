#'
#' Example data: Health outcomes of 2372 adults with and without diabetes
#'
#' A data frame with 2372 rows and 14 variables with diabetes status
#' \code{diabetes} and mortality status \code{endpoint}. For the purpose
#' of demonstrate, assume that we are interested in the association
#' between \code{diabetes} and \code{endpoint}. Other variables are
#' considered as possible confounding factors. The purposes of this dataset
#' is to illustrate those functions in \pkg{chest} and \pkg{allestimates} packages only.
#' Therefore, we assume it is a cohort design for Cox Proportional Hazard regression,
#' and a case-control design for logistic regression.
#'
#' @format A data frame with 2372 rows and 14 variables:
#' \describe{
#'   \item{Diabetes}{diabetes status 1: with diabetes 0: without diabetes }
#'   \item{Endpoint}{mortality status 1: reached end point, and 0: survived }
#'   \item{Age}{Age, in years}
#'   \item{Sex}{sex, 1: male, 2: Female}
#'   \item{BMI}{Body mass index}
#'   \item{Married}{marital status 1: married, 0: not}
#'   \item{Smoke}{smoking status 1: smoker, 0: non-smoker}
#'   \item{CVD}{cardiovascular disease 1: yes 0: no}
#'   \item{Cancer}{cancer 1: yes, 0: no}
#'   \item{Education}{education 1: high, 0: low}
#'   \item{Income}{income 1: high, 0: low}
#'   \item{t0}{time (age) at the start of the follow-up}
#'   \item{t1}{time (age) at the end of the follow-up}
#'   \item{mid}{matched set id, for conditional logistic regression }
#'   }
"diab_df"

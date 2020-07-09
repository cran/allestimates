#' Estimates all possible effect estimates using Cox Proportional Hazards regression models
#'
#' Estimates hazard ratios using Proportional Hazards Regression models
#' (\code{"coxph"} from \pkg{survival} package) from models with all
#' possible combinations of a list of variables.
#'
#' @export
#' @import survival
#' @param crude An object of \emph{formula} for initial model, generally crude model.
#' However, any other variables can also be included here as the initial model.
#' The left-hand side of ~ is the outcome of interest, and the variable on the
#' right-hand side of ~ is the exposure of the interest (either a treatment or a risk factor)
#' @param xlist A \emph{vector} of a list of variable names.
#' @param data \emph{Data frame}.
#' @param na_omit Remove all missing values. Default is \code{"na_omit = TRUE"}.
#' @param ... Further optional arguments.
#' @return A list of all effect estimates.
#' @seealso
#' \code{surival}
#' @examples
#' vlist <- c("Age", "Sex", "Married", "BMI", "Education", "Income")
#' results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
#' results
#' @name all_cox
all_cox <- function(crude, xlist, data, na_omit = TRUE, ...) {
  data <- data.frame(c(data[all.vars(as.formula(crude))], data[xlist]))
  if (na_omit) {
    data <- na.omit(data)
  }
  mod_0 <- survival::coxph(as.formula(crude), ..., data = data)
  m0 <- broom::tidy(mod_0, exponentiate = TRUE, conf.int = TRUE)
  p <- m0$p.value[1]
  estimate <- m0$estimate[1]
  conf_low <- m0$conf.low[1]
  conf_high <- m0$conf.high[1]
  n <- broom::glance(mod_0)$n
  df_0 <- data.frame(
    variables = "Crude",
    estimate, conf_low,
    conf_high, p, aic = AIC(mod_0), n
  )
  comb_lst <- unlist(lapply(
    seq_along(xlist),
    function(x) combn(xlist, x, simplify = FALSE)
  ),
  recursive = FALSE
  )
  variables <- sapply(comb_lst, toString)
  md_lst <- lapply(comb_lst, function(x) paste(crude, "+", paste(x, collapse = "+")))
  models <- lapply(md_lst, function(x) coxph(as.formula(x), ..., data = data))
  p <- unlist(lapply(models, function(x) tidy(x)$p.value[1]))
  estimate <- unlist(lapply(models, function(x) tidy(x, exponentiate = TRUE)$estimate[1]))
  conf_low <- unlist(lapply(models, function(x) tidy(x, exponentiate = TRUE)$conf.low[1]))
  conf_high <- unlist(lapply(models, function(x) tidy(x, exponentiate = TRUE)$conf.high[1]))
  n <- unlist(lapply(models, function(x) broom::glance(x)$n))
  aic <- unlist(lapply(models, function(x) stats::AIC(x)))
  df_coef <- data.frame(
    variables, estimate,
    conf_low, conf_high, p, aic, n
  )
  message("estimate: Hazard ratio")
  message(paste("Crude model:", crude))
  estimate <- rbind(df_0, df_coef)
  fun <- "all_cox"
  family <- "coxph"
  lst_ret <- list(estimate, xlist, fun, crude, family)
  names(lst_ret) <- c("estimate", "xlist", "fun", "crude", "family")
  lst_ret
}

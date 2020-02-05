#' Effect estimates from all possible models using \code{speedglm}
#'
#' This is a faster alternative to \code{all_glm}. \code{all_speedglm}
#' estimates odds ratios or rate ratios using generalized linear
#' models (\code{speedglm}) with all possible combinations of a list of variables
#' (potential confounding factors) specified in \code{xlist} argument.
#'
#' @export
#' @import speedglm
#' @param crude An object of \emph{formula} for initial model, generally crude model.
#' However, any other variables can also be included here as the initial model.
#' @param xlist A \emph{vector} of characters with variable names to be included in
#' as potential confounding factors.
#' @param data \emph{Data frame}.
#' @param na_omit Remove all missing values. Default is \code{"na_omit = TRUE"}.
#' @param family Description of the error distribution. Default is \code{binomial()}.
#' @param ... Further optional arguments.
#' @return A list of all effect estimates.
#' @seealso \pkg{speedglm}
#' @examples
#' vlist <- c("Age", "Sex", "Married", "Cancer", "CVD", "Education", "Income")
#' results <- all_speedglm(crude = "Endpoint ~ Diabetes", xlist = vlist, data = diab_df)
#' results$estimate
#' @name all_speedglm
all_speedglm <- function(crude, xlist, data,
                         family = binomial(), na_omit = TRUE, ...) {
  data <- data.frame(c(data[all.vars(as.formula(crude))], data[xlist]))
  if (na_omit) {
    data <- na.omit(data)
  }
  if (is.character(family)) {
    message("family =", family, "is incorrect (see 'help(all_speedglm)'.")
  }
  mod_0 <- speedglm::speedglm(as.formula(crude),
    family = family,
    data = data, ...
  )
  p <- as.numeric(as.character(tidy(mod_0)$p.value[2]))
  estimate <- exp(as.numeric(as.character(tidy(mod_0)$estimate[2])))
  conf_low <- exp(as.numeric(as.character(tidy(mod_0, conf.int = TRUE)$conf.low[2])))
  conf_high <- exp(as.numeric(as.character(tidy(mod_0, conf.int = TRUE)$conf.high[2])))
  aic <- AIC(mod_0)
  n <- stats::nobs(mod_0)
  df_0 <- data.frame(
    variables = "Crude",
    estimate,
    conf_low, conf_high, p, aic, n
  )
  comb_lst <- unlist(lapply(
    seq_along(xlist),
    function(x) combn(xlist, x, simplify = FALSE)
  ),
  recursive = FALSE
  )
  comb_str <- sapply(comb_lst, toString)
  md_lst <- lapply(
    comb_lst,
    function(x) paste(crude, "+", paste(x, collapse = "+"))
  )
  models <- lapply(
    md_lst,
    function(x) {
      speedglm::speedglm(as.formula(x),
        family = binomial(),
        data = data, ...
      )
    }
  )
  OR <- (lapply(
    models,
    function(x) as.character(tidy(x)$estimate[2])
  ))
  lb <- unlist(lapply(
    models,
    function(x) as.character(tidy(x, conf.int = TRUE)$conf.low[2])
  ))
  ub <- unlist(lapply(
    models,
    function(x) as.character(tidy(x, conf.int = TRUE)$conf.high[2])
  ))
  p <- unlist(lapply(models, function(x) as.character(tidy(x)$p.value[2])))
  aic <- unlist(lapply(models, function(x) AIC(x)))
  n <- unlist(lapply(models, function(x) stats::nobs(x)))
  estimate <- exp(as.numeric(OR))
  conf_low <- exp(as.numeric(lb))
  conf_high <- exp(as.numeric(ub))
  p <- as.numeric(p)
  df_coef <- data.frame(
    variables = comb_str, estimate,
    conf_low, conf_high, p, aic, n
  )
  message("estimate: Odds Ratio or Rate Ratio")
  message(paste("Crude model:", crude))
  estimate <- rbind(df_0, df_coef)
  fun <- "all_speedglm"
  family <- family$family
  lst_ret <- list(estimate, xlist, fun, crude, family)
  names(lst_ret) <- c("estimate", "xlist", "fun", "crude", "family")
  lst_ret
}

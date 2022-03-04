#' Estimates all possible effect estimates using \code{glm}
#'
#' \code{all_glm} estimates odds ratios or rate ratios using
#' generalized linear models (\code{glm}) with all
#' possible combinations of a list of variables (potential confounding factors).
#'
#' @export
#' @param crude An object of \emph{formula} for initial model, generally crude model.
#' However, any other variables can also be included here as the initial model.
#' @param xlist A \emph{vector} of a list of variable names (potential confounding factors).
#' @param data \emph{Data frame}.
#' @param family \emph{family} Description of the error distribution. Default is \code{"binomial"}.
#' @param na_omit Remove all missing values. Default is \code{"na_omit = TRUE"}.
#' @param ... Further optional arguments.
#' @return A list of all effect estimates.
#' @seealso \pkg{stats}
#' @examples
#' \dontrun{
#' diab_df$Overweight <- as.numeric(diab_df$BMI >= 25)
#' vlist <- c("Age", "Sex", "Income")
#' all_glm(crude = "Diabetes ~ Overweight", xlist = vlist, data = diab_df)
#' }
#' @name all_glm
#'
all_glm <- function(crude, xlist, data,
                    family = "binomial", na_omit = TRUE, ...) {
  data <- data.frame(c(data[all.vars(as.formula(crude))], data[xlist]))
  if (na_omit) {
    data <- na.omit(data)
  }
  mod_0 <- glm(as.formula(crude),
    family = family,
    data = data, ...
  )
  p <- tidy(mod_0)$p.value[2]
  estimate <- tidy(mod_0, exponentiate = TRUE)$estimate[2]
  conf_low <- tidy(mod_0, exponentiate = TRUE, conf.int = TRUE)$conf.low[2]
  conf_high <- tidy(mod_0, exponentiate = TRUE, conf.int = TRUE)$conf.high[2]
  aic <- AIC(mod_0)
  n <- stats::nobs(mod_0)
  df_0 <- data.frame(
    variables = "Crude",
    estimate, conf_low, conf_high, p, aic, n
  )
  exposure <- unlist(strsplit(crude, "~"))[2]
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
      glm(as.formula(x),
        family = family,
        data = data, ...
      )
    }
  )
  estimate <- unlist(lapply(
    models,
    function(x) {
      tidy(x,
        exponentiate = TRUE,
        conf.int = TRUE
      )$estimate[2]
    }
  ))
  conf_low <- unlist(lapply(
    models,
    function(x) {
      tidy(x,
        exponentiate = TRUE,
        conf.int = TRUE
      )$conf.low[2]
    }
  ))
  conf_high <- unlist(lapply(
    models,
    function(x) {
      tidy(x,
        exponentiate = TRUE,
        conf.int = TRUE
      )$conf.high[2]
    }
  ))
  p <- unlist(lapply(models, function(x) tidy(x)$p.value[2]))
  aic <- unlist(lapply(models, function(x) AIC(x)))
  n <- unlist(lapply(models, function(x) stats::nobs(x)))
  df_coef <- data.frame(
    variables = comb_str,
    estimate, conf_low, conf_high, p, aic, n
  )
  message("estimate: Odds Ratio or Rate Ratio")
  message(paste("Crude model:", crude))
  estimate <- rbind(df_0, df_coef)
  fun <- "all_glm"
  family <- mod_0$family$family
  lst_ret <- list(estimate, xlist, fun, crude, family)
  names(lst_ret) <- c("estimate", "xlist", "fun", "crude", "family")
  lst_ret
}

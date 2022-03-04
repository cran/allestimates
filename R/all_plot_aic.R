#' Draws scatter plot with all effect estimates against AIC
#'
#' \code{all_plot_aic()} generates a scatter plot with all effect estimates against AIC.
#'
#' @export
#' @param data \emph{Object} from \code{all_cox}, \code{all_glm}, \code{all_speedglm}, or \code{all_glm}, including all effect estimate values.
#' @param xlab \emph{Character} string for x-axis name. Default is \code{"AIC"}
#' @param ylab \emph{Character} string for y-axis name. Default depends on original model types.
#' @param title \emph{Character} for plot title. Default is \code{"NULL"}.
#' @return A \pkg{ggplot2} object: scatter plot
#' @examples
#' vlist <- c("Age", "Sex", "BMI")
#' results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
#' all_plot_aic(results)
#' @name all_plot_aic
all_plot_aic <- function(data,
                         xlab = "AIC",
                         ylab = NULL,
                         title = NULL) {
  result_df <- data$estimate
  if (is.null(ylab)) {
    if (data$fun == "all_cox") {
      ylab <- "Hazard ratio"
    } else if (data$fun == "all_lm") {
      ylab <- "Coefficient"
    } else if (data$family == "poisson") {
      ylab <- "Rate ratio"
    } else if (data$family == "binomial") {
      ylab <- "Odds ratio"
    }
    else {
      ylab <- "Effect estimates"
    }
  }
  hline <- ifelse(data$fun == "all_lm", 0, 1)
  df_scatter <- result_df %>%
    dplyr::mutate(p_value = p^(log(0.5) / log(0.05)))
  ggplot(data = df_scatter, aes(x = aic, y = estimate)) +
    geom_point(shape = 1) +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    geom_hline(yintercept = hline, linetype = "dashed")
}

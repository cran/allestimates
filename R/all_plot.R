#' Plot all effect estimates against p values
#'
#' \code{all_plot()}  generates a scatter plot with effect estimates of all possible models
#' again p values.
#' @export
#' @param data \emph{Object} from \code{all_cox}, \code{all_glm}, \code{all_speedglm}, or \code{all_glm}, including all effect estimate values.
#' @param xlabels \emph{Numeric vector} x-axis tick labels. Default is
#' \code{"c(0, 0.001, 0.01, 0.05, 0.2, 0.5, 1)"}.
#' @param xlim \emph{Vector} of 2 numeric values for x-axis limits. Default is \code{"c(0, 1)"}.
#' @param xlab \emph{Character} string for x-axis name. Default is \code{"P value"}.
#' @param ylim \emph{Vector} of 2 numeric values for y-axis limits.
#' @param ylab \emph{Character} string for y-axis name. Default depends on original model types.
#' @param yscale_log \emph{TRUE or FALSE} to re-scale y-axis to "log10". Default is \code{"FALSE"}.
#' @param title \emph{Character} for plot title. Default is \code{"NULL"}.
#' @return A \pkg{ggplot2} object: scatter plot
#' @examples
#' vlist <- c("Age", "Sex", "Married", "BMI", "Education", "Income")
#' results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
#' all_plot(results)
#' @name all_plot
all_plot <- function(data,
                     xlabels = c(0, 0.001, 0.01, 0.05, 0.2, 0.5, 1),
                     xlim = c(0, 1),
                     xlab = "P value",
                     ylim = NULL,
                     ylab = NULL,
                     yscale_log = FALSE,
                     title = NULL) {
  xbreaks <- sapply(xlabels, function(x) x^(log(0.5) / (log(0.05))))
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
    } else {
      ylab <- "Effect estimates"
    }
  }
  hline <- ifelse(data$fun == "all_lm", 0, 1)
  df_scatter <- result_df %>%
    dplyr::mutate(p_value = p^(log(0.5) / log(0.05)))
  if (is.null(ylim) & data$fun != "all_lm") {
    maxv <- max(max(1 / result_df$estimate), max(result_df$estimate))
    ylim <- c(1 / maxv, maxv)
  } else if (is.null(ylim) & data$fun == "all_lm") {
    ylim <- c(-(max(abs(result_df$estimate))), max(abs(result_df$estimate)))
  }
  if (yscale_log) {
    my_yscale <- scale_y_log10(limits = ylim)
  } else {
    my_yscale <- scale_y_continuous(limits = ylim)
  }
  ggplot(data = df_scatter, aes(x = p_value, y = estimate)) +
    geom_point(shape = 1) +
    scale_x_continuous(breaks = xbreaks, labels = xlabels, limits = xlim) +
    my_yscale +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    geom_vline(xintercept = 0.5, linetype = 2) +
    geom_hline(yintercept = hline, linetype = 2)
}

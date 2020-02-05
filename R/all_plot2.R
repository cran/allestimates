#' Plots all effect estimates against p values with each specific variable in the models
#'
#' \code{all_plot2()}  generates a panel of scatter plots with effect estimates of all possible models
#' again p values. Each plot includes effect estimates from all models including a specific variable.
#'
#' @export
#' @param data \emph{Object} from \code{all_cox}, \code{all_glm}, \code{all_speedglm}, or \code{all_glm}, including all effect estimate values.
#' @param xlabels \emph{numeric vector} x-axis tick labels. Default is
#' \code{"c(0, 0.001, 0.01, 0.05, 0.2, 0.5, 1)"}
#' @param xlim \emph{vector} of 2 numeric values for x-axis limits. Default is \code{"c(0, 1)"}.
#' @param xlab \emph{Character} string for x-axis name. Default is \code{"P value"}.
#' @param ylim \emph{vector} of 2 numeric values for y-axis limits.
#' @param ylab \emph{Character} string for y-axis name. Default depends on original model types.
#' @param yscale_log \emph{TRUE or FALSE} re-scale y-axis to "log10". Default is \code{"FALSE"}.
#' @param title \emph{Character} title. Default is \code{"NULL"}.
#' @return A \pkg{ggplot2} object: scatter plot
#' @examples
#' vlist <- c("Age", "Sex", "Married", "BMI", "Education", "Income")
#' results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
#' all_plot2(results)
#' @name all_plot2
all_plot2 <- function(data,
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
  mark_df <- sapply(data$xlist, function(x) x <- str_detect(result_df$variables, x)) * 1
  df_scatter <- cbind(df_scatter, mark_df)
  pivot_longer(df_scatter,
    cols = c(data$xlist)
  ) %>%
    mutate(value = factor(value, labels = c("Not included", "Included"))) %>%
    ggplot(aes(
      x = p_value, y = estimate,
      color = value,
      shape = value
    )) +
    geom_point() +
    my_yscale +
    scale_x_continuous(breaks = xbreaks, labels = xlabels, limits = xlim) +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_shape_manual(values = c(0, 1)) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    geom_hline(yintercept = hline, linetype = 2) +
    facet_wrap(. ~ name)
}

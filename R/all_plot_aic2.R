#' Draws multiple scatter plots of all effect estimates against AIC
#'
#' \code{all_plot_aic2()} draws multiple scatter plots of
#' all effect estimates against AIC. Each plot indicates if a specific
#' variable is included in the models.
#'
#' @export
#' @param data \emph{Object} from \code{all_cox}, \code{all_glm}, \code{all_speedglm}, or \code{all_glm}, including all effect estimate values.
#' @param xlab \emph{Character} string for x-axis name. Default is \code{"AIC"}.
#' @param ylab \emph{Character} string for y-axis name. Default depends on original model types.
#' @param title \emph{Character} for plot title. Default is \code{"NULL"}.
#' @return A  \pkg{ggplot2}  object: scatter plot.
#' @examples
#' \dontrun{
#' vlist <- c("Age", "Sex", "BMI", "Education", "Income")
#' results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
#' all_plot_aic(data = results)
#' }
#' @name all_plot_aic2
all_plot_aic2 <- function(data,
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
    } else {
      ylab <- "Effect estimates"
    }
  }
  hline <- ifelse(data$fun == "all_lm", 0, 1)
  mark_df <- sapply(
    data$xlist,
    function(x) {
      x <- str_detect(result_df$variables, x)
    }
  ) * 1
  df_scatter <- cbind(result_df, mark_df)
  pivot_longer(df_scatter,
    cols = c(data$xlist)
  ) %>%
    mutate(value = factor(value, labels = c("Not included", "Included"))) %>%
    ggplot2::ggplot(aes(
      x = aic, y = estimate,
      color = value, shape = value
    )) +
    ggplot2::geom_point() +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_shape_manual(values = c(0, 1)) +
    geom_hline(yintercept = hline, linetype = 2) +
    facet_wrap(. ~ name)
}

## ----set-options, echo=FALSE, cache=FALSE---------------------------------------------------------
options(width = 100)

## ---- include = FALSE-----------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  install.packages("allestimates")

## ----setup----------------------------------------------------------------------------------------
library(allestimates)

## ---- coxhp_result--------------------------------------------------------------------------------
vlist <- c("Age", "Sex", "Smoke", "Married", "BMI", "Income")
results <- all_cox(crude = "Surv(t0, t1, Endpoint) ~ Diabetes", xlist = vlist, data = diab_df)
results

## ---- cox_plot, fig.height=2.5, fig.width = 4-----------------------------------------------------
all_plot(results)

## ---- cox_plot2, fig.height=4, fig.width = 7------------------------------------------------------
all_plot2(results)

## ---- speedglm_plot, fig.height=2.5, fig.width = 4------------------------------------------------
diab_df$Overweight = as.numeric(diab_df$BMI >= 25)
vlist <- c("Age", "Sex", "Education","Married", "BMI", "Income")
results <- all_speedglm(crude = "Diabetes ~ Overweight", xlist = vlist, data = diab_df)
all_plot(results)

## ---- speedglm_plot2, fig.height=4, fig.width = 7-------------------------------------------------
all_plot2(results)

## ---- speedglm_plot_aic2, fig.height=4, fig.width = 7---------------------------------------------
all_plot_aic2(results)

## ---- lm_plot, fig.height=2.5, fig.width = 4------------------------------------------------------
vlist <- c("Age", "Sex", "Education", "Income", "Diabetes", "Smoke")
results <- all_lm(crude = "BMI ~ Married", xlist = vlist, data = diab_df)
all_plot(results)

## ---- lm_plot2, fig.height=4, fig.width = 7-------------------------------------------------------
all_plot2(results)


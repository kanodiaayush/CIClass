#' ---
#' title: "Causal Inference Course Final Project"
#' author: "Ayush Kanodia and Mitchell Linegar"
#' toc: true
#' toc_depth: 2
#' always_allow_html: true
#' output:
#'   pdf_document2:
#'     number_sections: true
#'     df_print: paged
#'     toc: true
#'     toc_depth: 2
#' ---
#### SETUP ####
#+ echo=FALSE
# set global options

local_dir <- "~/Dropbox/Athey/sherlock_oak" # points to  /oak/stanford/groups/athey/Stones2Milestones
data_dir <- sprintf("%s/basic_rec_system/obs_study", local_dir)

outcome_family <- outcome_family # based on whether your outcome is binary or not; input to glm call
outcome_family <- "gaussian" # if outcome not binary
treatment_family <- "gaussian"
outcome_type <- "class"
n_sims <- 20
prop_to_keep <- 1.0 # if you want to only run on a random sample of the data, if want to run on full data set to 1.0

# lambda <- c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5) #, 0.7, 1, 5, 10, 50, 100, 1000)
# lambda <- c(0.0001, 0.01)
lambda <- c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1, 5, 10, 50, 100, 1000)
prop_drop_rf <- c(0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, .5, .7, .8, .9)

propensity_bound <- c(0.01, 0.99)

library(here)
# devtools::install_github("hrbrmstr/hrbrthemes")
# library(hrbrthemes)
library(ggplot2)
theme_set(theme_classic())
library(data.table)
library(tidyverse)
library(broom)
library(grf)
library(sandwich)
devtools::install_github("swager/amlinear") # install amlinear package
library(amlinear)
library(stargazer)
library(policytree)

library(dplyr)       # Data manipulation (0.8.0.1)
library(fBasics)     # Summary statistics (3042.89)
library(corrplot)    # Correlations (0.84)
library(psych)       # Correlation p-values (1.8.12)
library(grf)         # Generalized random forests (0.10.2)
library(rpart)       # Classification and regression trees, or CART (4.1-13)
library(rpart.plot)  # Plotting trees (3.0.6)
library(treeClust)   # Predicting leaf position for causal trees (1.1-7)
library(car)         # linear hypothesis testing for causal tree (3.0-2)
library(remotes)    # Install packages from github (2.0.1)
library(readr)       # Reading csv files (1.3.1)
library(tidyr)       # Database operations (0.8.3)
library(tibble)      # Modern alternative to data frames (2.1.1)
library(knitr)       # RMarkdown (1.21)
library(kableExtra)  # Prettier RMarkdown (1.0.1)
library(ggplot2)     # general plotting tool (3.1.0)
library(haven)       # read stata files (2.0.0)
library(aod)         # hypothesis testing (1.3.1)
library(evtree)      # evolutionary learning of globally optimal trees (1.0-7)
library(estimatr)    # simple interface for OLS estimation w/ robust std errors ()

# For causal trees (Athey and Imbens, 2016)  version 0.0
remotes::install_github('susanathey/causalTree') # Uncomment this to install the causalTree package
library(causalTree)
remotes::install_github('grf-labs/sufrep') # Uncomment this to install the sufrep package
library(sufrep)

library(bookdown)
library(CBPS)

#### KNITR SETUP ####
#+ setup, include=FALSE
knitr::opts_chunk$set(
  echo = FALSE,
  cache = TRUE, 
  warning = FALSE,
  message = FALSE,
  cache.lazy = FALSE,
  dev = "cairo_pdf",
  fig.align = 'center',
  fig.pos = 'H')


pretty_print_inline <- function(x, decimals = 1) {
  if(is.integer(x)) {
    format(x,
           big.mark = ",")
  }
  else if(is.double(x)) {
    format(round(x, decimals),
           digits = decimals,
           nsmall = decimals,
           big.mark = ",",
           decimal.mark = ".")
  }
  else {
    x
  }
}
knitr::knit_hooks$set(inline = function(x, decimals = 1) {
  pretty_print_inline(x)
})

options(
  digits = 3,
  contrasts = rep("contr.treatment", 2),
  tinytex.verbose = TRUE
)
#### TUTORIAL FUNCTIONS ####
difference_in_means <- function(dataset) {
  # Filter treatment / control observations, pulls outcome variable as a vector
  y1 <- dataset %>% dplyr::filter(W == 1) %>% dplyr::pull(Y) # Outcome in treatment grp
  y0 <- dataset %>% dplyr::filter(W == 0) %>% dplyr::pull(Y) # Outcome in control group
  
  n1 <- sum(dataset[,"W"])     # Number of obs in treatment
  print(n1)
  n0 <- sum(1 - dataset[,"W"]) # Number of obs in control
  print(n0)
  
  # Difference in means is ATE
  tauhat <- mean(y1) - mean(y0)
  
  # 95% Confidence intervals
  se_hat <- sqrt( var(y0)/(n0-1) + var(y1)/(n1-1) )
  lower_ci <- tauhat - 1.96 * se_hat
  upper_ci <- tauhat + 1.96 * se_hat
  
  return(c(ATE = tauhat, lower_ci = lower_ci, upper_ci = upper_ci))
}

ate_condmean_ols <- function(dataset) {
  df_mod_centered = data.frame(scale(dataset, center = TRUE, scale = FALSE))
  
  # Running OLS with full interactions is like running OLS separately on
  # the treated and controls. If the design matrix has been pre-centered,
  # then the W-coefficient corresponds to the ATE.
  lm.interact = lm(Y ~ . * W, data = df_mod_centered)
  tau.hat = as.numeric(coef(lm.interact)["W"])
  se.hat = as.numeric(sqrt(vcovHC(lm.interact)["W", "W"]))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

ipw <- function(dataset, p) {
  W <- dataset$W
  Y <- dataset$Y
  G <- ((W - p) * Y) / (p * (1 - p))
  tau.hat <- mean(G)
  se.hat <- sqrt(var(G) / (length(G) - 1))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}

prop_score_ols <- function(dataset, p) {
  # Pulling relevant columns
  W <- dataset$W
  Y <- dataset$Y
  # Computing weights
  weights <- (W / p) + ((1 - W) / (1 - p))
  # OLS
  lm.fit <- lm(Y ~ W, data = dataset, weights = weights)
  tau.hat = as.numeric(coef(lm.fit)["W"])
  se.hat = as.numeric(sqrt(vcovHC(lm.fit)["W", "W"]))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}


aipw_ols <- function(dataset, p) {
  
  ols.fit = lm(Y ~ W * ., data = dataset)
  
  dataset.treatall = dataset
  dataset.treatall$W = 1
  treated_pred = predict(ols.fit, dataset.treatall)
  
  dataset.treatnone = dataset
  dataset.treatnone$W = 0
  control_pred = predict(ols.fit, dataset.treatnone)
  
  actual_pred = predict(ols.fit, dataset)
  
  G <- treated_pred - control_pred +
    ((dataset$W - p) * (dataset$Y - actual_pred)) / (p * (1 - p))
  tau.hat <- mean(G)
  se.hat <- sqrt(var(G) / (length(G) - 1))
  c(ATE=tau.hat, lower_ci = tau.hat - 1.96 * se.hat, upper_ci = tau.hat + 1.96 * se.hat)
}



#### OTHER FUNCTIONS ####
plot_prob <- function(prob, pred, model_name = "", data_name = ""){
  model = deparse(quote(Wmod)) %>% substr(0,1)
  stopifnot(model %in% c("W", "Y"))
  data_text = ifelse(data_name=="","", sprintf(", with %s data", data_name))
  ggplot(data.frame(prob, pred), aes(x = prob, y = pred)) + 
    geom_point(alpha = .01) +
    # geom_smooth() + 
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 4), se = TRUE) + 
    geom_abline(intercept = 0, slope = 1) + 
    labs(title = sprintf("Predicted %s Propensity vs Actual %s%s", model_name, model, data_text),
         x = sprintf("Predicted %s", model),
         y = sprintf("Observed %s", model)) + 
    xlim(0,1) + 
    ylim(0,1)
}

convert_to_prob <- function(x){
  1/(1 + exp(-x))
}

loglike <- function(pred, data){
  mean(data * log(pred) + (1 - data) * log(1 - pred))
}

setnames_W_Y <- function(df){
  W_Y_cols <- names(df)[grepl("W|Y", names(df), ignore.case = FALSE)]
  setnames(df, W_Y_cols, c("W", "Y"))
}

calculate_propensities <- function(df, continuous_treatment = FALSE){
  df_mod = copy(df)
  setDT(df_mod)
  Xmod = df_mod[,covariate_names, with=FALSE]
  Ymod = df_mod$Y
  Wmod = df_mod$W
  
  
  # takes a while to run
  if(continuous_treatment){
    p_W_CBPS <- CBPS(Wmod ~ as.matrix(Xmod), method = "exact")$fitted.values
    pW_logistic <- NULL
  } else {
    
    pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = treatment_family)
    pW_logistic <- predict(pW_logistic.fit, type = "response")
    # pW_logistic.fit.tidy <-  pW_logistic.fit %>% tidy()
    # hist(pW_logistic)
    

    # pW_logistic.fit.tidy <-  pW_logistic.fit %>% tidy()
    # hist(pW_logistic)
    
    p_W_CBPS <- NULL 
  }
  
  # no matter what do RF propensity scores
  # FIXME: does this actually make sense??
  pW_rf.fit = regression_forest(Xmod, Wmod, num.trees = 500)
  pW_rf <- predict(pW_rf.fit, newdata = Xmod)$predictions
  
  df_mod[, `:=`(p_W = pW_logistic,
                p_W_rf = pW_rf,
                p_W_CBPS = p_W_CBPS)]
  df_mod
}

causal_forest_from_df <- function(df){
  df_mod = copy(df)
  setDT(df_mod)
  Xmod = df_mod[,covariate_names, with=FALSE]
  Ymod = df_mod$Y
  Wmod = df_mod$W
  causal_forest(Xmod, Ymod, Wmod, num.trees = 500) 
}

run_all_models_on_df <- function(df_propensity, continuous_treatment = FALSE){

  df <- df_propensity[,c(covariate_names, "W", "Y"), with=FALSE]
  tauhat_rct <- difference_in_means(df)
  
  # linear only if binary treatment
  # CBPS only if continuous treatment
  # CF switches from average_treatment_effect to average_partial_effect if continuous treatment
  if (continuous_treatment){
    tauhat_logistic_ipw <- NULL
    tauhat_logistic_pscore_ols <- NULL
    tauhat_lin_logistic_aipw <- NULL
    
    tauhat_CBPS_ipw <- ipw(df, df_propensity$p_W_CBPS)
    tauhat_CBPS_pscore_ols <- prop_score_ols(df, df_propensity$p_W_CBPS)
    tauhat_lin_CBPS_aipw <- aipw_ols(df, df_propensity$p_W_CBPS)
    
    cf = causal_forest_from_df(df)
    
    ape_rf_aipw = grf::average_partial_effect(cf)
    
    tauhat_rf_aipw_ape = c(ATE=ape_rf_aipw["estimate"],
                           lower_ci=ape_rf_aipw["estimate"] - 1.96 * ape_rf_aipw["std.err"],
                           upper_ci=ape_rf_aipw["estimate"] + 1.96 * ape_rf_aipw["std.err"])
    
    tauhat_rf_aipw_ate <- NULL
  } else {
    tauhat_logistic_ipw <- ipw(df, df_propensity$p_W)
    tauhat_logistic_pscore_ols <- prop_score_ols(df, df_propensity$p_W)
    tauhat_lin_logistic_aipw <- aipw_ols(df, df_propensity$p_W)
    
    tauhat_CBPS_ipw <- NULL
    tauhat_CBPS_pscore_ols <- NULL
    tauhat_lin_CBPS_aipw <- NULL
    
    cf = causal_forest_from_df(df)
    
    ate_rf_aipw = average_treatment_effect(cf, target.sample="control")
    
    tauhat_rf_aipw_ate = c(ATE=ate_rf_aipw["estimate"],
                           lower_ci=ate_rf_aipw["estimate"] - 1.96 * ate_rf_aipw["std.err"],
                           upper_ci=ate_rf_aipw["estimate"] + 1.96 * ate_rf_aipw["std.err"])
    
    tauhat_rf_aipw_ape <- NULL
  }
  

  

  
  tauhat_rf_ipw = ipw(df, df_propensity$p_W_rf)
  

  
  tauhat_ols_rf_aipw = aipw_ols(df, df_propensity$p_W_rf)
  
  all_estimators = rbind(
    RCT_gold_standard = tauhat_rct,
    # naive_observational = tauhat_confounded,
    # linear_regression = tauhat_ols,
    logistic_propensity_weighted_regression = tauhat_logistic_pscore_ols,
    IPW_logistic = tauhat_logistic_ipw,
    AIPW_linear_plus_logistic = tauhat_lin_logistic_aipw,
    CBPS_propensity_weighted_regression = tauhat_CBPS_pscore_ols,
    IPW_CBPS = tauhat_CBPS_ipw,
    AIPW_linear_plus_CBPS = tauhat_lin_CBPS_aipw,
    IPW_forest = tauhat_rf_ipw,
    AIPW_ate_causal_forest = tauhat_rf_aipw_ate,
    AIPW_ape_causal_forest = tauhat_rf_aipw_ape, 
    AIPW_linear_plus_forest = tauhat_ols_rf_aipw
    # IPW_lasso = tauhat_lasso_ipw#,
    # FIXME: add this above
    # AIPW_lasso = tauhat_lasso_aipw,
    # FIXME: add this above
    # approx_residual_balance = tauhat_lasso_balance
  )
  all_estimators = data.frame(all_estimators)
  all_estimators$ci_length = all_estimators$upper_ci - all_estimators$lower_ci
  all_estimators
  # round(as.matrix(all_estimators), 3) %>% stargazer(summary = FALSE, header = FALSE)
}

#### LOAD DATA ####


freadom <- fread(sprintf("%s/treatment_effects.csv", data_dir))
data_dir <- "Personal"
freadom <- fread(sprintf("%s/treatment_effects.csv", data_dir))

# cleaning
# some W not really continuous treatment - just different levels
freadom[,W_reading_time] %>% hist()
freadom$W_wordcount %>% hist()

# for now only look at 'medium' and 'high' levels of W_reading_time
freadom <- freadom[W_reading_time==12.5 | W_reading_time==7.5]
# look only at people who actually opened the SOD
freadom <- freadom[W_utility != 0]
# create binary treatment variables
freadom[, W_reading_time_high := ifelse(W_reading_time==max(W_reading_time), 1, 0)]
freadom[, W_wordcount_high := ifelse(W_wordcount>median(W_wordcount), 1, 0)]
freadom[, W_wordcount_AB := ifelse(W_wordcount > median(W_wordcount), "A", "B")]
# make Utility (completion) binary
freadom[, `:=`(W_utility_high = ifelse(W_utility==3, 1, 0),
               Y_utility_high = ifelse(Y_utility==3, 1, 0))]

covariate_names <- colnames(freadom)[!grepl("W|Y", colnames(freadom), ignore.case = FALSE)]

# construct various data sets we'll need
df_time_finish <- freadom[,c(covariate_names, "W_reading_time_high", "Y_utility"), with=FALSE]
df_time_finish %>% setnames_W_Y()

df_wcount_finish <- freadom[,c(covariate_names, "W_wordcount_high", "Y_utility"), with=FALSE]
df_wcount_finish %>% setnames_W_Y()

df_time_lag <- freadom[,c(covariate_names, "W_reading_time_high", "Y_next_view_lag"), with=FALSE]
df_time_lag %>% setnames_W_Y()

df_wcount_lag <- freadom[,c(covariate_names, "W_wordcount_high", "Y_next_view_lag"), with=FALSE]
df_wcount_lag %>% setnames_W_Y()

df_finish_lag <- freadom[,c(covariate_names, "W_utility_high", "Y_next_view_lag"), with=FALSE]
df_finish_lag %>% setnames_W_Y()


#### ANALYSIS PER DATASET ####

#+ propensity plots
# Calculate and Plot Overlap
df_time_finish_overlap <- calculate_propensities(df_time_finish)
df_time_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOP Reading Time Propensity, Finish SOD Outcome")

df_wcount_finish_overlap <- calculate_propensities(df_wcount_finish)
df_wcount_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOP Word Count Propensity, Finish SOD Outcome")

df_time_lag_overlap <- calculate_propensities(df_time_lag)
df_time_lag_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOP Reading Time Propensity, Time to Next Session Outcome")

df_wcount_lag_overlap <- calculate_propensities(df_wcount_lag)
df_wcount_lag_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOP Word Count Propensity, Time to Next Session Outcome")


df_finish_lag_overlap <- calculate_propensities(df_finish_lag)
df_finish_lag_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Predicted Probability of Completing SOP")

# Subset if necessary
df_time_finish_overlap <- df_time_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_wcount_finish_overlap <- df_wcount_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_time_lag_overlap <- df_time_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_wcount_lag_overlap <- df_wcount_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_finish_lag_overlap <- df_finish_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]



#### READING TIME FINISH SOD ANALYSIS ####
#' \newpage
#' # Rate at Which Users Finish Stories of the Day by Length {#finish_sod_by_length}
#+ finish_sod_time_ate, results='asis'

df_time_finish_models <- run_all_models_on_df(df_time_finish_overlap)
round(as.matrix(df_time_finish_models), 3) %>% stargazer(summary = FALSE, header = FALSE)

#### WORD COUNT FINISH SOD ANALYSIS ####
#' \newpage
#' # Rate at Which Users Finish Stories of the Day by Word Count {#finish_sod_by_wcount}
#+ finish_sod_wcount_ate, results='asis'

df_wcount_finish_models <- run_all_models_on_df(df_wcount_finish_overlap)
round(as.matrix(df_wcount_finish_models), 3) %>% stargazer(summary = FALSE, header = FALSE)



#### READING TIME LAG SOD ANALYSIS ####
#' \newpage
#' # Rate at Which Users Finish Stories of the Day by Length {#lag_sod_by_length}
#+ lag_sod_time_ate, results='asis'

df_time_lag_models <- run_all_models_on_df(df_time_lag_overlap)
round(as.matrix(df_time_lag_models), 3) %>% stargazer(summary = FALSE, header = FALSE)

#### WORD COUNT LAG SOD ANALYSIS ####
#' \newpage
#' # Rate at Which Users Finish Stories of the Day by Word Count {#lag_sod_by_wcount}
#+ lag_sod_wcount_ate, results='asis'

df_wcount_lag_models <- run_all_models_on_df(df_wcount_lag_overlap)
round(as.matrix(df_wcount_lag_models), 3) %>% stargazer(summary = FALSE, header = FALSE)

#### INSTRUMENTAL VARIABLES AND CAUSAL FORESTS ####
#' \newpage
#' # Instrumental Variables and Causal Forests
#' In this section we estimate ATE by using instrumental variables. HTO BE DONE??


#### FINISH SOD LAG NEXT SESSION ANALYSIS ####
#' # Effect of Finishing SOD on Time to Next Session {#finish_sod_lag}
#+ finish_sod_lag_avg, results='asis'
df_finish_lag_models <- run_all_models_on_df(df_finish_lag_overlap)
round(as.matrix(df_finish_lag_models), 3) %>% stargazer(summary = FALSE, header = FALSE)

#### FINISH SOD LAG NEXT SESSION CATE ANALYSIS ####
#' We now estimate the CATE, and use it to construct quartiles. We then report the ATE as estimated with AIPW from our causal forest estimate across quartiles. 
#+ results='asis'
cf <- causal_forest_from_df(df_finish_lag_overlap)

oob_pred <- predict(cf, estimate.variance=TRUE)
df_finish_lag_overlap$cate <- oob_pred$predictions

# Manually creating subgroups
num_tiles <- 4  # ntiles = CATE is above / below the median
df_finish_lag_overlap$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))


# ATE estimates
df_finish_lag_overlap_qtile <- df_finish_lag_overlap[,.(
  avg_cf_cate = mean(cate)
),.(ntile)][order(ntile)]

df_finish_lag_overlap_qtile

estimated_aipw_ate <- lapply(
  seq(num_tiles), function(w) {
    ate <- average_treatment_effect(cf, subset = df_finish_lag_overlap$ntile == w)
  })
estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
estimated_aipw_ate$ntile <- factor(seq(num_tiles))
setDT(estimated_aipw_ate)
setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
# join CATE estimates together
df_finish_lag_overlap_qtile <- df_finish_lag_overlap_qtile[estimated_aipw_ate, on = .(ntile)]
df_finish_lag_overlap_qtile %>% stargazer(summary = FALSE, header = FALSE)

X <- as.matrix(freadom[, c(covariate_names), with=FALSE])
Y <- as.matrix(freadom[, c("Y_utility"), with=FALSE])
W <- as.matrix(freadom[, c("W_wordcount_AB"), with=FALSE])
multi.forest <- multi_causal_forest(X = X, Y = Y, W = W)
Gamma.matrix <- double_robust_scores(multi.forest)

train <- sample(1:n, 200)
opt.tree <- policy_tree(X[train, ], Gamma.matrix[train, ], depth = 2)
opt.tree

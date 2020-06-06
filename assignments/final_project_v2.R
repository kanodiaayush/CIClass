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
#' fontsize: 12pt
#' linestretch: 1.5
#' ---
#### SETUP ####
#+ echo=FALSE
# set global options

local_dir <- "~/Dropbox/Athey/sherlock_oak" # points to  /oak/stanford/groups/athey/Stones2Milestones

data_dir <- "."
# data_dir <- sprintf("%s/basic_rec_system/obs_study", local_dir)


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

# multiplier of treeplot font size
treeplot_font_size_mult <- 0.8

#### PACKAGES ####
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

devtools::install_github("grf-labs/policytree")
library(policytree)



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

forest_from_df <- function(df, FUN = causal_forest, Wname = "W", Yname = "Y"){
  df_mod = copy(df)
  setDT(df_mod)
  Xmod = df_mod[,covariate_names, with=FALSE] %>% as.matrix()
  Ymod = df_mod[,Yname, with=FALSE] %>% as.matrix()
  Wmod = df_mod[,Wname, with=FALSE] %>% as.matrix()
  FUN(Xmod, Ymod, Wmod, num.trees = 500) 
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
    
    cf = forest_from_df(df)
    
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
    
    cf = forest_from_df(df)
    
    ate_rf_aipw = average_treatment_effect(cf, target.sample="control")
    
    tauhat_rf_aipw_ate = c(ATE=ate_rf_aipw["estimate"],
                           lower_ci=ate_rf_aipw["estimate"] - 1.96 * ate_rf_aipw["std.err"],
                           upper_ci=ate_rf_aipw["estimate"] + 1.96 * ate_rf_aipw["std.err"])
    
    tauhat_rf_aipw_ape <- NULL
  }
  

  

  
  tauhat_rf_ipw = ipw(df, df_propensity$p_W_rf)
  

  
  tauhat_ols_rf_aipw = aipw_ols(df, df_propensity$p_W_rf)
  
  all_estimators = rbind(
    Difference_in_Means = tauhat_rct,
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
  # round(as.matrix(all_estimators), 3) %>% knitr::kable(format = "latex")
}

#### LOAD DATA ####
#+ echo=TRUE
# trip-level data 
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
# freadom[, W_wordcount_AB := ifelse(W_wordcount > median(W_wordcount), "A", "B")]
# make Utility (completion) binary
freadom[, `:=`(W_utility_high = ifelse(W_utility==3, 1, 0),
               Y_utility_high = ifelse(Y_utility==3, 1, 0))]

# make words read quartiles
freadom[, W_wordcount_qtile := ntile(W_wordcount, 4)]

# create words read count
freadom[, Y_wordsread := case_when(Y_utility==0 ~ 0,
                                 Y_utility==1 ~ W_wordcount * .25,
                                 Y_utility==2 ~ W_wordcount * .5,
                                 Y_utility==3 ~ W_wordcount)]

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

df_wcount_qtile_wordsread <- freadom[,c(covariate_names, "W_wordcount_qtile", "Y_wordsread"), with=FALSE]
df_wcount_qtile_wordsread %>% setnames_W_Y()

# aggregate level data
freadom_agg <- fread(sprintf("%s/treatment_effects_agg.csv", data_dir))

# don't need to do any filtering here because everything is an aggregate

# create binary treatment variables
freadom_agg[, W_reading_time_high := ifelse(W_reading_time>median(W_reading_time), 1, 0)]
freadom_agg[, W_wordcount_high := ifelse(W_wordcount>median(W_wordcount), 1, 0)]

# make Utility (completion) binary
freadom_agg[, `:=`(W_utility_high = ifelse(W_utility>1.5, 1, 0),
               Y_utility_high = ifelse(Y_utility>1.5, 1, 0))]

# make words read quartiles
freadom_agg[, W_wordcount_qtile := ntile(W_wordcount, 4)]

# create words read count
freadom_agg[, Y_wordsread := Y_utility/3 * W_wordcount]

# construct aggregate data
agg_df_time_finish <- freadom_agg[,c(covariate_names, "W_reading_time_high", "Y_utility"), with=FALSE]
agg_df_time_finish %>% setnames_W_Y()

agg_df_wcount_finish <- freadom_agg[,c(covariate_names, "W_wordcount_high", "Y_utility"), with=FALSE]
agg_df_wcount_finish %>% setnames_W_Y()

agg_df_time_lag <- freadom_agg[,c(covariate_names, "W_reading_time_high", "Y_next_view_lag"), with=FALSE]
agg_df_time_lag %>% setnames_W_Y()

agg_df_wcount_lag <- freadom_agg[,c(covariate_names, "W_wordcount_high", "Y_next_view_lag"), with=FALSE]
agg_df_wcount_lag %>% setnames_W_Y()

agg_df_finish_lag <- freadom_agg[,c(covariate_names, "W_utility_high", "Y_next_view_lag"), with=FALSE]
agg_df_finish_lag %>% setnames_W_Y()

agg_df_wcount_qtile_wordsread <- freadom_agg[,c(covariate_names, "W_wordcount_qtile", "Y_wordsread"), with=FALSE]
agg_df_wcount_qtile_wordsread %>% setnames_W_Y()

#### ANALYSIS PER DATASET ####

#+ propensity plots, echo=TRUE
# Calculate and Plot Overlap
df_time_finish_overlap <- calculate_propensities(df_time_finish)
df_time_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOD Reading Time Propensity")

df_wcount_finish_overlap <- calculate_propensities(df_wcount_finish)
df_wcount_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOD Word Count Propensity")

df_finish_lag_overlap <- calculate_propensities(df_finish_lag)
df_finish_lag_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Predicted Probability of Completing SOD Propensity")

df_wcount_qtile_wordsread_overlap <- calculate_propensities(df_wcount_qtile_wordsread)
df_wcount_qtile_wordsread_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Assigned Word Count propensity")

# aggregate data

agg_df_time_finish_overlap <- calculate_propensities(agg_df_time_finish)
agg_df_time_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOD Reading Time Propensity, 
          Averaged Over Each User's Trips")

agg_df_wcount_finish_overlap <- calculate_propensities(agg_df_wcount_finish)
agg_df_wcount_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Above Median SOD Word Count Propensity, 
          Averaged Over Each User's Trips")

agg_df_finish_lag_overlap <- calculate_propensities(agg_df_finish_lag)
agg_df_finish_lag_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Predicted Probability of Completing SOD Propensity, 
          Averaged Over Each User's Trips")

agg_df_wcount_qtile_wordsread_overlap <- calculate_propensities(agg_df_wcount_qtile_wordsread)
agg_df_wcount_qtile_wordsread_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram() + 
  ggtitle("Assigned Word Count Propensity, 
          Averaged Over Each User's Trips")

# Subset if necessary
df_time_finish_overlap <- df_time_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_wcount_finish_overlap <- df_wcount_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_time_lag_overlap <- df_time_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_wcount_lag_overlap <- df_wcount_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_finish_lag_overlap <- df_finish_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
df_wcount_qtile_wordsread_overlap <- df_wcount_qtile_wordsread_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]


agg_df_time_finish_overlap <- agg_df_time_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
agg_df_wcount_finish_overlap <- agg_df_wcount_finish_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
agg_df_time_lag_overlap <- agg_df_time_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
agg_df_wcount_lag_overlap <- agg_df_wcount_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
agg_df_finish_lag_overlap <- agg_df_finish_lag_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]
agg_df_wcount_qtile_wordsread_overlap <- agg_df_wcount_qtile_wordsread_overlap[p_W_rf %between% propensity_bound & p_W %between% propensity_bound]


#' # User-Session Level Analysis
#### READING TIME FINISH SOD ANALYSIS ####
#' \newpage
#' ## Effect of Higher than Average Estimated Reading Time on SOD Completion {#finish_sod_by_length}
#+ finish_sod_time_ate, results='asis', echo=TRUE

df_time_finish_models <- run_all_models_on_df(df_time_finish_overlap)
round(as.matrix(df_time_finish_models), 3) %>% knitr::kable(format = "latex")

#### WORD COUNT FINISH SOD ANALYSIS ####
#' \newpage
#' ## Effect of Higher than Average Estimated Word Count on SOD Completion {#finish_sod_by_wcount}
#+ finish_sod_wcount_ate, results='asis', echo=TRUE

df_wcount_finish_models <- run_all_models_on_df(df_wcount_finish_overlap)
round(as.matrix(df_wcount_finish_models), 3) %>% knitr::kable(format = "latex")



#### READING TIME LAG SOD ANALYSIS ####
#' \newpage
#' ## Effect of Higher than Average Estimated Reading Time SOD on Time to Next Session {#lag_sod_by_length}
#+ lag_sod_time_ate, results='asis', echo=TRUE

df_time_lag_models <- run_all_models_on_df(df_time_lag_overlap)
round(as.matrix(df_time_lag_models), 3) %>% knitr::kable(format = "latex")

#### WORD COUNT LAG SOD ANALYSIS ####
#' \newpage
#' ## Effect of Higher than Average Word Count SOD on Time to Next Session {#lag_sod_by_wcount}
#+ lag_sod_wcount_ate, results='asis', echo=TRUE

df_wcount_lag_models <- run_all_models_on_df(df_wcount_lag_overlap)
round(as.matrix(df_wcount_lag_models), 3) %>% knitr::kable(format = "latex")


#### WORD COUNT QUARTILE SOD WORDS READ ANALYSIS ####
#' \newpage
#' # Number of Words Read  {#lag_sod_by_wcount_qtile}
#+ lag_sod_wcount_ate_qtl, results='asis', echo=TRUE

df_wcount_lag_models <- run_all_models_on_df(df_wcount_lag_overlap)
round(as.matrix(df_wcount_lag_models), 3) %>% knitr::kable(format = "latex")



#### WORD COUNT SOD LAG NEXT SESSION CATE ANALYSIS ####
#' \newpage
#' # Word Count CATE on Time to Next Session  
#' We now estimate the CATE, and use it to construct quartiles of user-sessions. We then report the ATE as estimated with AIPW from our causal forest estimate across quartiles. 
#+ results='asis', echo=TRUE
cf <- forest_from_df(df_wcount_lag_overlap)

oob_pred <- predict(cf, estimate.variance=TRUE)
df_wcount_lag_overlap$cate <- oob_pred$predictions

# Manually creating subgroups
num_tiles <- 4  # ntiles = CATE is above / below the median
df_wcount_lag_overlap$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))


# ATE estimates
df_wcount_lag_overlap_qtile <- df_wcount_lag_overlap[,.(
  avg_cf_cate = mean(cate)
),.(ntile)][order(ntile)]

df_wcount_lag_overlap_qtile

estimated_aipw_ate <- lapply(
  seq(num_tiles), function(w) {
    ate <- average_treatment_effect(cf, subset = df_wcount_lag_overlap$ntile == w)
  })
estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
estimated_aipw_ate$ntile <- factor(seq(num_tiles))
setDT(estimated_aipw_ate)
setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
# join CATE estimates together
df_wcount_lag_overlap_qtile <- df_wcount_lag_overlap_qtile[estimated_aipw_ate, on = .(ntile)]
df_wcount_lag_overlap_qtile %>% knitr::kable(format = "latex")


#' # User-level Analysis (Averaged over User-Sessions)

#### READING TIME FINISH SOD AGGREGATE ANALYSIS ####
#' \newpage
#' ## User Average Effect of Higher than Average Estimated Reading Time on SOD Completion {#finish_sod_by_length_agg}
#+ finish_sod_time_ate_agg, results='asis', echo=TRUE

agg_df_time_finish_models <- run_all_models_on_df(agg_df_time_finish_overlap)
round(as.matrix(agg_df_time_finish_models), 3) %>% knitr::kable(format = "latex")

#### WORD COUNT FINISH SOD AGGREGATE ANALYSIS ####
#' \newpage
#' ## User Average Effect of Higher than Average Estimated Word Count on SOD Completion {#finish_sod_by_wcount_agg}
#+ finish_sod_wcount_ate_agg, results='asis', echo=TRUE

agg_df_wcount_finish_models <- run_all_models_on_df(agg_df_wcount_finish_overlap)
round(as.matrix(agg_df_wcount_finish_models), 3) %>% knitr::kable(format = "latex")



#### READING TIME LAG SOD AGGREGATE ANALYSIS ####
#' \newpage
#' ## User Average Effect of Higher than Average Estimated Reading Time SOD on Time to Next Session {#lag_sod_by_length_agg}
#+ lag_sod_time_ate_agg, results='asis', echo=TRUE

agg_df_time_lag_models <- run_all_models_on_df(agg_df_time_lag_overlap)
round(as.matrix(agg_df_time_lag_models), 3) %>% knitr::kable(format = "latex")

#### WORD COUNT LAG SOD AGGREGATE ANALYSIS ####
#' \newpage
#' ## User Average Effect of Higher than Average Word Count SOD on Time to Next Session {#lag_sod_by_wcount_agg}
#+ lag_sod_wcount_ate_agg, results='asis', echo=TRUE

agg_df_wcount_lag_models <- run_all_models_on_df(agg_df_wcount_lag_overlap)
round(as.matrix(agg_df_wcount_lag_models), 3) %>% knitr::kable(format = "latex")


#### FINISH SOD LAG NEXT SESSION CATE AGGREGATE ANALYSIS ####
#' We now estimate the CATE, and use it to construct quartiles. We then report the ATE as estimated with AIPW from our causal forest estimate across quartiles. 
#+ results='asis', echo=TRUE
cf <- forest_from_df(agg_df_wcount_lag_overlap)

oob_pred <- predict(cf, estimate.variance=TRUE)
agg_df_wcount_lag_overlap$cate <- oob_pred$predictions

# Manually creating subgroups
num_tiles <- 4  # ntiles = CATE is above / below the median
agg_df_wcount_lag_overlap$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))


# ATE estimates
agg_df_wcount_lag_overlap_qtile <- agg_df_wcount_lag_overlap[,.(
  avg_cf_cate = mean(cate)
),.(ntile)][order(ntile)]

agg_df_wcount_lag_overlap_qtile

estimated_aipw_ate <- lapply(
  seq(num_tiles), function(w) {
    ate <- average_treatment_effect(cf, subset = agg_df_wcount_lag_overlap$ntile == w)
  })
estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
estimated_aipw_ate$ntile <- factor(seq(num_tiles))
setDT(estimated_aipw_ate)
setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
# join CATE estimates together
agg_df_wcount_lag_overlap_qtile <- agg_df_wcount_lag_overlap_qtile[estimated_aipw_ate, on = .(ntile)]
agg_df_wcount_lag_overlap_qtile %>% knitr::kable(format = "latex")



#### POLICY TREE ####
#' # Optimal Policy Trees
#+ echo=TRUE
X <- as.matrix(freadom[, c(covariate_names), with=FALSE])
# Y <- as.matrix(freadom[, c("Y_utility"), with=FALSE])
# W <- as.matrix(freadom[, c("W_wordcount_high"), with=FALSE])
# multi.forest <- multi_causal_forest(X = X, Y = Y, W = W)
#' ## Policy Tree for Effect of SOD Word Count Quartile on Words Read per Session
#+ echo=TRUE
multi.forest <- forest_from_df(freadom, multi_causal_forest, "W_wordcount_qtile", "Y_wordsread")
Gamma.matrix <- double_robust_scores(multi.forest)

opt.tree <- policy_tree(X, Gamma.matrix, depth = 2)
plot(opt.tree, tweak = treeplot_font_size_mult, title = "Total Words Read")

#' ## Policy Tree for Effect of SOD Word Count Quartile on SOD Completion
#+ echo=TRUE
multi.forest <- forest_from_df(freadom, multi_causal_forest, "W_wordcount_qtile", "Y_utility")
Gamma.matrix <- double_robust_scores(multi.forest)

opt.tree <- policy_tree(X, Gamma.matrix, depth = 2)
plot(opt.tree, tweak = treeplot_font_size_mult, title = "Total Words Read")

#' ## Policy Tree for User-Level Average Effect of SOD Word Count Quartile on Words Read per Session
#+ echo=TRUE
X <- as.matrix(freadom_agg[, c(covariate_names), with=FALSE])
multi.forest <- forest_from_df(freadom_agg, multi_causal_forest, "W_wordcount_qtile", "Y_wordsread")
Gamma.matrix <- double_robust_scores(multi.forest)

opt.tree <- policy_tree(X, Gamma.matrix, depth = 2)
plot(opt.tree, tweak = treeplot_font_size_mult)

#' ## Policy Tree for User-Level Average Effect of SOD Word Count Quartile on SOD Completion
#+ echo=TRUE
X <- as.matrix(freadom_agg[, c(covariate_names), with=FALSE])
multi.forest <- forest_from_df(freadom_agg, multi_causal_forest, "W_wordcount_qtile", "Y_utility")
Gamma.matrix <- double_robust_scores(multi.forest)

opt.tree <- policy_tree(X, Gamma.matrix, depth = 2)
plot(opt.tree, tweak = treeplot_font_size_mult)


#### WRITEUP ####
#' \newpage
#' # Writeup Introduction  
#' 
#' Improvements to childhood literacy have been linked to numerous positive outcomes, including economic and social benefits. In this paper, we use data from a mobile application aimed at improving childhood reading outcomes. The app is primarily used by schoolgoing children from junior kindergarten until grade 3 to read and interact with stories.  
#' 
#' This work takes advantage of the application’s “Story of the Day” (hereafter SOD) feature. Stories of the Day are featured prominently on the app, and users read the SOD on approximately 33% of days they use the app. Our analysis focuses entirely on these stories. Several SODs are available to be assigned to users on each day, and vary primarily by estimated reading time and word count. The assignment of story of the day is generic and not personalised, and different stories are shown each day. As a result, if we consider a user opening the app as an exogenous random decision on a given day, since the stories shown to students are different each day, this gives us exogenous treatments for length of stories shown to students in terms of reading time and number of words in story. We measure the effect of this treatment on reading outcomes. Our identifying motivation is that longer stories reduce the probability of a child reading a story.  
#' 
#' Even if this effect is true on average, this does not mean that longer stories have negative effects on all users. As such, we examine CATEs across a variety of groups in Section 5.
#' 
#' Finally, in Section 6 we attempt to maximize aggregate reading time by identifying the optimal policy, and summarize possible gains from targeted assignment of Stories of the Day by length.
#' 
#' # Data Description
#' ## User Covariates:
#' In our dataset, we have a bunch of covariates describing users. These include age, grade, statistics about usage such as books read, experience points gained on the app while using it, total time spent on the app, covariates about how well a child answered questions related to their readings, and their reading interests.  
#' 
#' ## Treatment Definitions  
#' 
#' We use the following treatments for analysis:  
#' 
#' * Suggested Reading Time for a given story: The app includes suggested reading times for a story in one of five choices. We examine only stories where the estimated reading time was either 7.5 or 12.5 minutes, as estimated reading time is not continuous. These two values of estimated reading time account for 95% of all user-trips.   
#' * Number of words for a given story: We parsed the stories shown on the app to get the number of words in each story, and we use this as a treatment variable. We divide observations into those where the number of words is below and above the median, giving us a treatment and a control. We also do another analysis with multiple treatments where the treatment is characterised by the quantile in which its number of words falls; giving us 4 treatment conditions.  
#' 
#' We recall that as the assignment of the Story of the Day is at-random, so is the assignment of word count and estimated reading time. Note: We examine only users-trips where the user started reading the Story of the Day, as otherwise the user would have no estimate of the length or time required to read the story.
#' 
#' ## Outcomes  
#' 
#' We examine three outcomes in our analysis: whether the user finished reading their assigned Story of the Day, the length of time until their next session, and the estimated number of words they read. In addition to conducting our analyses at the user-session level, we also estimate user-level aggregate treatment effects. For each user we aggregate over all sessions, averaging session-level binary treatment status and outcome. 
#' 
#' # Average Treatment Effect  
#' 
#' First, we measure the average treatment effect of our various treatments on our various outcomes. For this, as we learnt in the class, we use various methods.
#' 
#' ## Models  
#' 
#' For each set of models, we compare results from the following methods:  
#' * Simple Difference in Means   
#' * Logistic Propensity Weighted Linear Regression  
#' * IPW Estimator using Logistic Propensities  
#' * IPW Estimator using Regression Forest Propensities  
#' * AIPW Estimator using Logistic Propensities and Linear Regression  
#' * AIPW Estimator using Causal Forests for both outcome and propensity models  
#' 
#' ## Results for each Treatment  
#' 
#' ### Suggested Reading Time  
#' 
#' At the user-session level, estimated reading time has no significant effect on a user’s likelihood of completing their story of the day; confidence intervals include zero in almost all models in Table 1.1. 
#' On the other hand, Table 4.1 shows the estimated reading time has a significant and negative effect on the probability of a user completing the SOD at the user level: when assigned readings that take longer, users are less likely to complete an average SOD. This is true and significant across all model specifications.   
#'  
#' We see that model results are similar with some variation in confidence intervals. The largest treatment effect in magnitude is reported by the IPW estimator with logistic propensities. The AIPW estimator with causal forests propensity model and a linear regression outcome model reports a more conservative treatment effect and also reports the minimum confidence interval length, while the difference in means estimator is a close second. However, we believe that the former is more reliable because it tries to control for confounding. We note that results across models are similar across treatments and outcomes. 
#' At the user-session level in Table 1.3, we do not see significant effects on time until next SOD started under any model. 
#' 
#' At the user-level in Table 4.3, our causal forest model has confidence intervals that include zero, while all other models indicate a negative effect of higher average reading time on time until the next SOD attempt. Otherwise, the results are analogous to the effect on SOD completion, with AIPW estimators with Causal Forest Propensities giving the smallest Confidence Intervals.
#' 
#' ### Word Count  
#' 
#' At the user-session level in Table 1.2, we see that longer SOD by word-count had a slightly negative (but statistically significant) effect on finishing SOD. This is the case for the difference in means estimate of the treatment effect as well as all AIPW methods, which are generally more reliable than the (conflicting evidence presented by) the IPW methods. This matches our initial hypothesis.  
#' In Table 1.4 we see some evidence for the effect of SOD word count on time until next login: while the causal forest confidence intervals include zero, all other AIPW models are significant and show that longer passages by word-length increases the average time until the next session.  
#' 
#' At the user-level in Table 4.2, users who were on average assigned longer readings by word count are generally less likely to finish their SOD. This is true of all AIPW methods other than the AIPW estimated with linear propensity scores. 
#' Table 4.4 shows that there is conflicting evidence for the effect of longer word-count SODs on time until the next SOD view: the causal forest shows a significant and positive relationship between the two, while almost all other methods have confidence intervals that include zero. 
#'  
#' # Conditional Average Treatment Effect  
#' 
#' We now examine the CATE of word count on time until the next SOD attempt.  These results are detailed in Table 5.1. 
#' 
#' We experimented with a number of models for estimating CATEs, but only present those CATE estimates which were done using causal forests. 
#' We first use causal forests to estimate the CATE, and use our CATE estimates to construct quartiles of user-sessions. We then report the ATE as estimated with AIPW from our causal forest estimate across quartiles. We find significant effects in the first and fourth quartiles in opposite directions: for those users in the bottom quartile of estimated CATE, we find that an increase in word count has a negative and statistically significant effect on time until next SOD attempt; for the highest quartile we observe a significant, positive relationship between these two quantities. If our goal is to increase total amount read, this indicates that we may be able to do so by targeting users with SODs of different length. 
#' 
#' 
#' # Optimal Policy  
#' 
#' In this section, we divide story word counts into 4 quartiles, leading to 4 different treatments. We use multi treatment causal forests to estimate treatment effects for each of these. We consider two outcomes; SOD completion and Number of Words Read by a user. We calculate treatment effects for these outcomes at both the user level as well as the user session level. As we learnt in class, our goal is to learn optimal policy assignments for students in this section. Here, a policy would denote a decision to show a user a story of a certain length.
#' 
#' ## Model  
#' 
#' We use the Policy Tree model. This model does a greedy search over the covariate space, in the style of a regression forest, but exhausting all possible split points (with some approximation if desired). We show depth 2 policy trees learnt for policy allocation with this model.
#' 
#' ## Results  
#' 
#' We note that the Policy Tree model tells us that for users with more experience, we should adopt the policy of showing them stories which are longer to maximize their reading utilities (section 5.2); the second level split shows this. This is an interesting finding which shows that more experienced users mature over time and prefer to read lengthier stories.
#' 
#' For number of words read as outcome, the model tells us that we should prefer to show lengthier stories to all users (quantile 3 and 4 in section 5.1) in general
#' 
#' The results at the user level are more mixed. As in section 5.3, the model suggests showing longer stories to users with higher qa accuracy to increase number of words read.
#' Probably the most interesting plot is in section 5.4. Here, we note that the policy prescribed for students with less than 55 questions is quartile 3 for students with lower qa accuracy, and quartile 4 for the rest. However, for students who have answered an unusually high number of questions, the policy prescribed is to show shorter utility stories to maximize SOD completion. We hypothesize that this is possible because this class represents unusually enthusiastic students who wish to maximize number of stories read/questions completed, and the easiest way for them to do this is to read shorter stories. This suggests corrective plans for such strategies which may not be pedagogically the most beneficial.
#' 
#' # Conclusion  
#' 
#' In this paper, we analyzed the impact of exposing children to stories of varied lengths, and testing effects of this treatment on various reading outcomes. We note that there is, as we hypothesized, an increased completion of stories when the length of the stories is shorter. This is not true, however, in bringing back users to the app more often. At the same time, we also find that this increase is not true ubiquitously for all users. As we see with the CATE estimates, there are subpopulations for whom the treatment effect may work in the opposite direction. 
#' Further, as the Policy Tree analysis shows, we may want to particularly target different subpopulations in different ways; to some beginners we want to show stories which are shorter, so that they complete these easier stories and ride their wave of accomplishment to read more of these. For more advanced users, it is more beneficial to expose them to longer, more complex reads, to maximize reading.
#' We utilized multiple methods for ATE, CATE and Optimal Policy Estimation, which we will return to in future work.  




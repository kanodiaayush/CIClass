#' ---
#' title: "Causal Inference Course Final Project"
#' author: "Ayush Kanodia and Mitchell Linegar"
#' toc: true
#' toc_depth: 2
#' always_allow_html: true
#' output:
#'   pdf_document:
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
  
  n1 <- sum(df[,"W"])     # Number of obs in treatment
  n0 <- sum(1 - df[,"W"]) # Number of obs in control
  
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
#### LOAD DATA ####

df <- fread(sprintf("%s/utility_dataset.csv", data_dir))

covariate_names <- colnames(df)[!colnames(df) %in% c("W", "Y")]

#### DATA WORK ####
# Drop rows containing missing values
df <- na.omit(df)

# Converting all columns to numerical
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))))

# coerce to data.table for future analysis
setDT(df)

#### RCT ANALYSIS ####
#' ## RCT Analysis
#' We now report the (presumably true) treatment effect $\hat{\tau}$ from the randomized experiment:  
#+ echo=TRUE
tauhat_rct <- difference_in_means(df)
print(tauhat_rct)

#### TESTING ASSUMPTIONS ####
#' ## Testing Assumptions
#' Here we test some of our traditional causal inference assumptions. 
#+ echo=FALSE
df_mod = copy(df)
setDT(df_mod)
Xmod = df_mod[,covariate_names, with=FALSE]
Ymod = df_mod$Y
Wmod = df_mod$W

#' \newpage
#' As a first step, we plot logistic predictions of the probabilities our treatment $pW$ and our outcome $pY$ (which is binary). 
#' We see that treatment assignment appears to follow a normal distribution, and that our outcome has an average unconditional probability of  `r mean(Ymod)`.  
#+ echo=TRUE
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = outcome_family)
pW_logistic <- predict(pW_logistic.fit, type = "response")
pW_logistic.fit.tidy <-  pW_logistic.fit %>% tidy()
hist(pW_logistic)

pY_logistic.fit <- glm(Ymod ~ as.matrix(Xmod), family = outcome_family)
pY_logistic <- predict(pY_logistic.fit, type = "response")
pY_logistic.fit.tidy <-  pY_logistic.fit %>% tidy()
hist(pY_logistic)

df_mod[, `:=`(p_Y = pY_logistic,
              p_W = pW_logistic)]

#' Some summary statistics of our data:
#+ results='asis'
stargazer(df_mod, header=FALSE)

#' We now produce a plot comparing predicted and actual treatment assignment. 
#' This plot is provided mostly for comparison (this is the plot the tutorial has);
#' future plots of this nature will be done with `ggplot` to make their options more explicit.  
#+ echo=TRUE
{plot(smooth.spline(pW_logistic, Wmod, df = 4))
  abline(0, 1)}

#### RCT ANALYSIS ####
#' ## RCT Analysis
#' We now report the (presumably true) treatment effect $\hat{\tau}$ from the randomized experiment:  
#+ echo=TRUE
tauhat_rct <- difference_in_means(df)
print(tauhat_rct)

#### LOGISTIC PROPENSITY SCORES ####
#+ echo=TRUE
# df_mod <- copy(df)
Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
Ymod = df_mod$Y
Wmod = df_mod$W
XWmod = cbind(Xmod, Wmod)

# Computing the propensity score by logistic regression of W on X.
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = outcome_family)
pW_logistic <- predict(pW_logistic.fit, type = "response")

df_mod[, logistic_propensity := pW_logistic]

#### OVERLAP ####
#' We now plot (logistic) propensity scores, showing that we still have overlap after removing observations. 
#' We may be somewhat concerned about the small number of observations with propensities close to one;
#' so remove all observations with propensity score outside of `r propensity_bound[1]` and `r propensity_bound[2]` to fix this. 
#' We then re-estimate the propensity model.  
#' We first show overlap before truncating.  
#+ echo=TRUE
overlap <- df_mod %>% ggplot(aes(x=logistic_propensity,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()
overlap

#' We now truncate, and plot truncated overlap. 
#+ echo=TRUE
df_mod <- df_mod[logistic_propensity %between% propensity_bound]

Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
Ymod = df_mod$Y
Wmod = df_mod$W
XWmod = cbind(Xmod, Wmod)

overlap <- df_mod %>% ggplot(aes(x=logistic_propensity,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()
overlap

#### PREDICTING PROPENSITIES AND OUTCOMES, ORIGINAL AND EXPANDED DATA ####
# some of this is used to calculate bias function, hence the ordering

# logistic model
# original data
pW_logistic.fit <- glm(Wmod ~ Xmod, family = outcome_family)
pW_logistic <- predict(pW_logistic.fit, type = "response")

# original data
pY_logistic.fit <- glm(Ymod ~ XWmod, family = outcome_family)
pY_logistic <- predict(pY_logistic.fit, type = "response")

# lasso expanded data, code provided by TA
# original data
pW_glmnet.fit.model = glmnet::cv.glmnet(Xmod, Wmod, lambda = lambda, family = outcome_family, type.measure = "class", keep=TRUE) 
pY_glmnet.fit.model = glmnet::cv.glmnet(Xmod, Ymod, lambda = lambda, family = outcome_family, type.measure = outcome_type, keep=TRUE) 
# expanded data

# demonstration of lasso fit across lambdas:
pW_lasso = pW_glmnet.fit.model$fit.preval[, pW_glmnet.fit.model$lambda == pW_glmnet.fit.model$lambda.min] %>% convert_to_prob()
pW_lasso.min = pW_glmnet.fit.model$fit.preval[, pW_glmnet.fit.model$lambda == min(pW_glmnet.fit.model$lambda)] %>% convert_to_prob()
pW_lasso.max = pW_glmnet.fit.model$fit.preval[, pW_glmnet.fit.model$lambda == max(pW_glmnet.fit.model$lambda)] %>% convert_to_prob()
pW_lasso.rand = pW_glmnet.fit.model$fit.preval[, pW_glmnet.fit.model$lambda == base::sample(pW_glmnet.fit.model$lambda, 1)] %>% convert_to_prob()

# random forest
pW_rf.fit = regression_forest(Xmod, Wmod, num.trees = 500)
pY_rf.fit = regression_forest(Xmod, Ymod, num.trees = 500)

# pW_rf = pW_rf.fit$predictions
# pY_rf = pY_rf.fit$predictions
pW_rf = predict(pW_rf.fit, newdata = Xmod) %>% as.matrix
pY_rf = predict(pY_rf.fit, newdata = Xmod) %>% as.matrix

# CF
cf = causal_forest(Xmod, Ymod, Wmod, num.trees = 500)

#### BIAS FUNCTION ####
#' Next we plot the bias function $b(X)$ following Athey, Imbens, Pham and Wager (AER P&P, 2017, Section IIID). 
#' We plot $b(x)$ for all units in the sample, and see that the bias seems evenly distributed around zero.  
#' We see that bias for most observations is close to zero.  
#+ echo=TRUE

mu_avg <- function(treated, df){df[W==treated, mean(Y)]}
mu <- function(treated, df){df[W==treated, mean(pY)]}

B <- function(df, treatment_model, outcome_model, outcome_type = "response"){
  # have to supply models so that can estimate counterfactual predictions given an alternative treatment assignment
  # note that this will NOT work for lasso model, attempt to warn of this misbehavior:
  if (grepl("lasso|rf|cf", deparse(quote(treatment_model)))){
    simpleMessage("The predict method appears to be broken for lasso models estimated using glmnet::cv.glmnet.
                  You may want to try another predictive model.")
  }
  df = copy(df)
  p = df[,mean(W)]
  mu0 <- df[W==0,mean(Y)]
  mu1 <- df[W==1,mean(Y)]
  
  pY_w0 <- predict(outcome_model, newdata = df[,.SD, .SDcols = !c('W', 'Y')][, W := 0], type = outcome_type)
  pY_w1 <- predict(outcome_model, newdata = df[,.SD, .SDcols = !c('W', 'Y')][, W := 1], type = outcome_type)
  pW <- predict(treatment_model, newdata = df[,.SD, .SDcols = !c('W', 'Y')], type = "response")
  df[, `:=`(W = NULL, pY_w0 = pY_w0, pY_w1 = pY_w1, pW = pW)]
  
  
  df[, b := (pW - p) * (p * (pY_w0 - mu0) + (1 - p) * (pY_w1 - mu1))]
  
  return(df[,.(b)])
}

df_mod_bias <- B(df_mod, pW_logistic.fit, pY_logistic.fit)
ggplot(df_mod_bias, aes(x = b)) + geom_histogram() + labs(title = "Histogram of per-observation b(x)")

#### ESTIMATING ATE INTRO ####
#' \newpage
#' ## Estimating the ATE  
#' In this section we explore various methods for estimating the ATE. We explore the following methods:  
#' 1. inverse propensity weighting via logistic regression  
#' 2. direct regression analysis via OLS  
#' 3. traditional double robust analysis via augmented inverse-propensity score weighting that combines the above two estimators.  
#' We also re-run the above methods after expanding the data to include all interactions of all of the covariates, 
#' and re-estimate outcome and proensity models using the original linear model, as well as running lasso and random forest models on the expanded data. 

#### ATE CALCULATIONS: ORIGINAL DATA ####

# linear models
tauhat_ols <- ate_condmean_ols(df_mod)

# linear models
tauhat_logistic_ipw <- ipw(df_mod, pW_logistic)
tauhat_pscore_ols <- prop_score_ols(df_mod, pW_logistic)
tauhat_lin_logistic_aipw <- aipw_ols(df_mod, pW_logistic)

# lasso
tauhat_lasso_ipw <- ipw(df_mod, pW_lasso)
tauhat_pscore_lasso <- prop_score_ols(df_mod, pW_lasso)
tauhat_lasso_logistic_aipw <- aipw_ols(df_mod, pW_lasso)
# FIXME: add this
# prior code to add:
# Xmod.for.lasso = cbind(Wmod, Xmod, (2 * Wmod - 1) * Xmod)
# glmnet.fit.outcome = amlinear:::crossfit.cv.glmnet(Xmod.for.lasso, Ymod,
#                                                    penalty.factor = c(0, rep(1, ncol(Xmod.for.lasso) - 1)))
# lasso.yhat.control = amlinear:::crossfit.predict(glmnet.fit.outcome,
#                                                  cbind(0, Xmod, -Xmod))
# lasso.yhat.treated = amlinear:::crossfit.predict(glmnet.fit.outcome,
#                                                  cbind(1, Xmod, Xmod))
# The lasso AIPW estimator. Here, the inference is justified via
# orthogonal moments.
# G = lasso.yhat.treated - lasso.yhat.control +
#   Wmod / pW_lasso * (Ymod - lasso.yhat.treated) -
#   (1 - Wmod) / (1 - pW_lasso) * (Ymod - lasso.yhat.control)
# tau.hat = mean(G)
# se.hat = sqrt(var(G) / length(G))
# tauhat_lasso_aipw = c(ATE=tau.hat,
#                       lower_ci=tau.hat-1.96*se.hat,
#                       upper_ci=tau.hat+1.96*se.hat)

# FIXME: add this
# balancing.weights = amlinear::balance_minimax(Xmod, Wmod, zeta = 0.5)
# G.balance = lasso.yhat.treated - lasso.yhat.control +
#   balancing.weights * (Ymod - Wmod * lasso.yhat.treated
#                        - (1 - Wmod) * lasso.yhat.control)
# tau.hat = mean(G.balance)
# se.hat = sqrt(var(G.balance) / length(G.balance))
# tauhat_lasso_balance = c(ATE=tau.hat,
#                          lower_ci=tau.hat-1.96*se.hat,
#                          upper_ci=tau.hat+1.96*se.hat)

# RF
tauhat_rf_ipw = ipw(df_mod, pW_rf)
ate_rf_aipw = average_treatment_effect(cf)
tauhat_rf_aipw = c(ATE=ate_rf_aipw["estimate"],
                   lower_ci=ate_rf_aipw["estimate"] - 1.96 * ate_rf_aipw["std.err"],
                   upper_ci=ate_rf_aipw["estimate"] + 1.96 * ate_rf_aipw["std.err"])
tauhat_ols_rf_aipw = aipw_ols(df_mod, pW_rf)

#### COMPARING ATE ACROSS MODELS ####
#' ## Comparing ATE Across Models with Original Data
#' Finally, we compare ATE across various models. We see that AIPW forest methods performs the best across the original and interacted data,
#' though propensity weighted regression performed on the original data.
#+ echo=FALSE
all_estimators = rbind(
  RCT_gold_standard = tauhat_rct,
  naive_observational = tauhat_confounded,
  linear_regression = tauhat_ols,
  propensity_weighted_regression = tauhat_pscore_ols,
  IPW_logistic = tauhat_logistic_ipw,
  AIPW_linear_plus_logistic = tauhat_lin_logistic_aipw,
  IPW_forest = tauhat_rf_ipw,
  AIPW_forest = tauhat_rf_aipw,
  AIPW_linear_plus_forest = tauhat_ols_rf_aipw,
  IPW_lasso = tauhat_lasso_ipw#,
  # FIXME: add this above
  # AIPW_lasso = tauhat_lasso_aipw,
  # FIXME: add this above
  # approx_residual_balance = tauhat_lasso_balance
)
all_estimators = data.frame(all_estimators)
all_estimators$ci_length = all_estimators$upper_ci - all_estimators$lower_ci
print(round(as.matrix(all_estimators), 3))

#### Honest Causal Tree ####
#' # Part 3: Heterogeneous Treatment Effects in Randomized Experiments
#' ## Data Definition
#' For this section, we return to `df`, our original, unaltered dataset. 
#' We use random sampling to divide the dataset into three datasets, call them df_tr, df_est, and df_test as in the tutorial.
#+ echo=TRUE
library(causalTree)
library(gt)

df <- df %>% mutate(id = row_number())
#Create training set
df_tr <- df %>% sample_frac(.50)
df_split  <- anti_join(df, df_tr, by = 'id')
df_est <- df_split %>% sample_frac(.50)
df_test <- anti_join(df_split, df_est, by = 'id')
#Create test set
nrow(df)
nrow(df_tr) + nrow(df_est) + nrow(df_test)

fmla_ct <- paste("factor(Y) ~", paste(covariate_names, collapse = " + "))

#' ## Honest Causal Tree
#' We now create a factor variable for the leaves in samples df_tr, df_est and df_test, and run linear regressions that estimate the treatment effect magnitudes and standard errors for each leaf in each sample.  
#' Our code follows that from the tutorial.  
#+ causal_tree, echo=TRUE
ct_unpruned <- honest.causalTree(
  formula=fmla_ct,            # Define the model
  data=df_tr,              # Subset used to create tree structure
  est_data=df_est,            # Which data set to use to estimate effects
  treatment=df_tr$W,       # Splitting sample treatment variable
  est_treatment=df_est$W,     # Estimation sample treatment variable
  split.Rule="CT",            # Define the splitting option
  cv.option="TOT",            # Cross validation options
  cp=0,                       # Complexity parameter
  split.Honest=TRUE,          # Use honesty when splitting
  cv.Honest=TRUE,             # Use honesty when performing cross-validation
  minsize=15,                 # Min. number of treatment and control cases in each leaf
  HonestSampleSize=nrow(df_est)) # Num obs used in estimation after building the tree
# Table of cross-validated values by tuning parameter.
ct_cptable <- as.data.frame(ct_unpruned$cptable)
# Obtain optimal complexity parameter to prune tree.
selected_cp <- which.min(ct_cptable$xerror)
optim_cp_ct <- ct_cptable[selected_cp, "CP"]
# Prune the tree at optimal complexity parameter.
ct_pruned <- prune(tree=ct_unpruned, cp=optim_cp_ct)
tauhat_ct_est <- predict(ct_pruned, newdata=df_est)
tauhat_ct_test <- predict(ct_pruned, newdata=df_test)
tauhat_ct_tr <- predict(ct_pruned, newdata=df_tr)
# Create a factor column 'leaf' indicating leaf assignment
num_leaves <- length(unique(tauhat_ct_est))  # There are as many leaves as there are predictions
df_est$leaf <- factor(tauhat_ct_est, labels = seq(num_leaves))
df_test$leaf <- factor(tauhat_ct_test, labels = seq(num_leaves))
df_tr$leaf <- factor(tauhat_ct_tr, labels = seq(num_leaves))
# Run the regression
ols_ct_est <- lm_robust(Y ~ 0 + leaf + W:leaf, data=df_est)
ols_ct_test <- lm_robust(Y ~ 0 + leaf + W:leaf, data=df_test)
ols_ct_tr <- lm_robust(Y ~ 0 + leaf + W:leaf, data=df_tr)
ols_ct_summary_est <- summary(ols_ct_est)
ols_ct_summary_test <- summary(ols_ct_test)
ols_ct_summary_tr <- summary(ols_ct_tr)
te_summary_est <- coef(ols_ct_summary_est)[(num_leaves+1):(2*num_leaves), c("Estimate")] %>% round(3)
te_summary_test <- coef(ols_ct_summary_test)[(num_leaves+1):(2*num_leaves), c("Estimate")] %>% round(3)
te_summary_tr <- coef(ols_ct_summary_tr)[(num_leaves+1):(2*num_leaves), c("Estimate")] %>% round(3)
te_summary_est_se <- coef(ols_ct_summary_est)[(num_leaves+1):(2*num_leaves), c( "Std. Error")] %>% round(3)
te_summary_test_se <- coef(ols_ct_summary_test)[(num_leaves+1):(2*num_leaves), c( "Std. Error")] %>% round(3)
te_summary_tr_se <- coef(ols_ct_summary_tr)[(num_leaves+1):(2*num_leaves),c("Std. Error")] %>% round(3)
row_names <-  sprintf("Leaf %d",1:num_leaves)
estimates <- cbind(row_names,paste0(te_summary_est,"(",te_summary_est_se,")"),paste0(te_summary_test,"(",te_summary_test_se,")"), paste0(te_summary_tr,"(",te_summary_tr_se,")")) %>% as.data.frame()
gt(estimates,rowname_col = "row_names")  %>%  
  tab_header(
    title = md("Per Leaf ATE"),
    subtitle = paste0("Updated to ",Sys.Date())) %>%
  cols_label(V2 = "df_est",
             V3 = "df_test",
             V4 = "df_tr") %>%
  as_latex()

#### Partial dependence plots ####
#' ### Partial Dependence Plots
#' All code here follows directly from tutorial.  
#' 
#' It may also be interesting to examine how our CATE estimates behave when we change a single covariate, while keeping all the other covariates at a some fixed value. In the plot below we evaluate a variable of interest across quantiles, while keeping all other covariates at their median (see the RMarkdown source for code).  
#' **Note:** It is important to recognize that in the following plots and tables, we may be evaluating the CATE at $x$ values in regions where there are few or no data points. Also, it may be the case that varying some particular variable while keeping others fixed may just not be very interesting. 
#' For example, in the _welfare_ dataset (the dataset we used last problem set!), we will not see a lot difference when we change `partyid` if we keep `polviews` fixed at their median value. It might be instructive to re-run this tutorial without using the variable `partyid`.


#' We see that there is some movement in the mean CATE across age groups, but they are all essentially the same. 
#' We observe a similar pattern across `hrs1`, our other covariate we use to form subgroups.  
#+ echo=TRUE
df_tr <- df %>% sample_frac(.50)
df_split  <- anti_join(df, df_tr, by = 'id')
cf <- causal_forest(
  X = as.matrix(df_tr[,covariate_names]),
  Y = df_tr$Y,
  W = df_tr$W,
  num.trees=200) 
if (dataset_name == "welfare") {
  var_of_interest = "age"
} else {
  # Selecting a continuous variable, if available, to make for a more interesting graph
  continuous_variables <- sapply(covariate_names, function(x) length(unique(df_tr[, x])) > 5)
  # Select variable for single variable plot
  var_of_interest <- ifelse(sum(continuous_variables) > 0,
                            covariate_names[continuous_variables][1],
                            covariate_names[1])
  # Select variables for two variable plot
  vars_of_interest <- c(var_of_interest,
                        ifelse(sum(continuous_variables) > 1,
                               covariate_names[continuous_variables][2],
                               covariate_names[covariate_names != var_of_interest][1]))
}
# Create a grid of values: if continuous, quantiles; else, plot the actual values
is_continuous <- (length(unique(df_tr[,var_of_interest])) > 5) # crude rule for determining continuity
if(is_continuous) {
  x_grid <- quantile(df_tr[,var_of_interest], probs = seq(0, 1, length.out = 5))
} else {
  x_grid <- sort(unique(df_tr[,var_of_interest]))
}
df_grid <- setNames(data.frame(x_grid), var_of_interest)
# For the other variables, keep them at their median
other_covariates <- covariate_names[which(covariate_names != var_of_interest)]
df_median <- data.frame(lapply(df_tr[,other_covariates], median))
df_eval <- crossing(df_median, df_grid)
# Predict the treatment effect
pred <- predict(cf, newdata=df_eval[,covariate_names], estimate.variance=TRUE)
df_eval$tauhat <- pred$predictions
df_eval$se <- sqrt(pred$variance.estimates)
# Change to factor so the plotted values are evenly spaced
df_eval[, var_of_interest] <- as.factor(round(df_grid[, var_of_interest], digits = 4))
# Descriptive labeling
label_description <- ifelse(is_continuous, '\n(Evaluated at quintiles)', '')

#+ echo=TRUE
# Plot
df_eval %>%
  mutate(ymin_val = tauhat-1.96*se) %>%
  mutate(ymax_val = tauhat+1.96*se) %>%
  ggplot() +
  geom_line(aes_string(x=var_of_interest, y="tauhat",group=1), color="red") +
  geom_errorbar(aes_string(x=var_of_interest,ymin="ymin_val", ymax="ymax_val", width=.2),color="blue") +
  xlab(paste0("Effect of ", var_of_interest, label_description)) +
  ylab("Predicted Treatment Effect") +
  theme_linedraw() +
  theme(axis.ticks = element_blank())


df_eval[,c(var_of_interest,"tauhat","se")] %>% gt() %>%
  as_latex()

#' Now for `hrs1`. This code exactly mirrors that above, so we omit it.  
#+ 
if (dataset_name == "welfare") {
  var_of_interest = "hrs1"
} else {
  # Selecting a continuous variable, if available, to make for a more interesting graph
  continuous_variables <- sapply(covariate_names, function(x) length(unique(df_tr[, x])) > 5)
  # Select variable for single variable plot
  var_of_interest <- ifelse(sum(continuous_variables) > 0,
                            covariate_names[continuous_variables][1],
                            covariate_names[1])
  # Select variables for two variable plot
  vars_of_interest <- c(var_of_interest,
                        ifelse(sum(continuous_variables) > 1,
                               covariate_names[continuous_variables][2],
                               covariate_names[covariate_names != var_of_interest][1]))
}
# Create a grid of values: if continuous, quantiles; else, plot the actual values
is_continuous <- (length(unique(df_tr[,var_of_interest])) > 5) # crude rule for determining continuity
if(is_continuous) {
  x_grid <- quantile(df_tr[,var_of_interest], probs = seq(0, 1, length.out = 5))
} else {
  x_grid <- sort(unique(df_tr[,var_of_interest]))
}
df_grid <- setNames(data.frame(x_grid), var_of_interest)
# For the other variables, keep them at their median
other_covariates <- covariate_names[which(covariate_names != var_of_interest)]
df_median <- data.frame(lapply(df_tr[,other_covariates], median))
df_eval <- crossing(df_median, df_grid)
# Predict the treatment effect
pred <- predict(cf, newdata=df_eval[,covariate_names], estimate.variance=TRUE)
df_eval$tauhat <- pred$predictions
df_eval$se <- sqrt(pred$variance.estimates)
# Change to factor so the plotted values are evenly spaced
df_eval[, var_of_interest] <- as.factor(round(df_grid[, var_of_interest], digits = 4))
# Descriptive labeling
label_description <- ifelse(is_continuous, '\n(Evaluated at quintiles)', '')

#' 
#+
df_eval %>%
  mutate(ymin_val = tauhat-1.96*se) %>%
  mutate(ymax_val = tauhat+1.96*se) %>%
  ggplot() +
  geom_line(aes_string(x=var_of_interest, y="tauhat", group = 1), color="red") +
  geom_errorbar(aes_string(x=var_of_interest,ymin="ymin_val", ymax="ymax_val", width=.2),color="blue") +
  xlab(paste0("Effect of ", var_of_interest, label_description)) +
  ylab("Predicted Treatment Effect") +
  theme_linedraw() +
  theme(axis.ticks = element_blank())

df_eval[,c(var_of_interest,"tauhat","se")] %>% gt() %>%
  as_latex()


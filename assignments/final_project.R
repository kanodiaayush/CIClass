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

outcome_family <- "binomial" # based on whether your outcome is binary or not; input to glm call
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

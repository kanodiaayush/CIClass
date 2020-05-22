#' ---
#' title: "Problem Set 2"
#' author: "Mitchell Linegar"
#' toc: true
#' toc_depth: 2
#' includes:
#'   pdf_document:
#'     number_sections: true
#'     df_print: paged
#'     toc: yes
#'     toc_depth: 2
#'     toc_float: true
#' ---
#### SETUP ####
#' echo=TRUE
# set global options
dataset_name <- "welfare"
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
# Pick any data set from the list above for parts I and II of the tutorial

if (dataset_name == "social") {
  # Load raw data
  df <- readr::read_csv(file = "https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Social/ProcessedData/socialpressnofact.csv", na = character())
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "outcome_voted"
  treatment_variable_name <- "treat_neighbors"
  covariate_names <-c("sex", "yob", "city", "hh_size", "totalpopulation_estimate",
                      "percent_male", "median_age", "percent_62yearsandover",
                      "percent_white", "percent_black", "percent_asian", "median_income",
                      "employ_20to64", "highschool", "bach_orhigher","percent_hispanicorlatino")
} else if (dataset_name == "charitable") {
  # Load raw data
  df <- readr::read_csv(file = "https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Charitable/ProcessedData/charitable_withdummyvariables.csv", na = character())
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "out_gavedum"
  treatment_variable_name <- "treatment"
  covariate_names <- c("pwhite", "pblack", "ave_hh_sz", "median_hhincome",  "pop_propurban", "female", "couple", "red0", "redcty")
  # This dataset encodes missing values as -999. Let's change these to NA.
  df[df == -999] <- NA
} else if (dataset_name == "mobilization") {
  # Create a temporary file and download the data to it.
  temp <- tempfile()
  download.file("https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Mobilization/ProcessedData/mobilization_with_unlisted.zip", temp)
  # # Load raw data
  df <- readr::read_csv((unz(temp, "mobilization_with_unlisted.csv")), na = character())
  unlink(temp)
  df <- subset(df, treatment ==1 | treatment==0)
  #  Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "vote02"
  treatment_variable_name <- "treatment"
  covariate_names <- c("comp_ia", "age", "female", "vote98", "vote00",
                       "state", "newreg", "comp_mi", # These covariates are commented out because
                       "county","st_sen","st_hse", # they are constant in the first 5000 rows.
                       "persons", "competiv")     
  # Thus, they can be included if using the entire dataset.sfor
} else if (dataset_name == "secrecy") {
  # Load raw data
  df <- readr::read_csv(file = "https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Secrecy/ProcessedData/ct_ballotsecrecy_processed.csv", na = character())
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "turnoutindex_12"
  treatment_variable_name <- "anysecrecytreatment"
  covariate_names <- c("v_cong_general_10", "v_pres_primary_12",	"v_cong_primary_12","v_pres_general_12",
                       "town1_block",	"town2_block", "town3_block","town4_block", "town5_block",
                       "town6_block",	"i_grp_addr_1",	"i_grp_addr_2",	"i_grp_addr_3",	"i_grp_addr_4",
                       "dem",	"rep",	"female",	"age_md",	"age_sq_md")
} else if (dataset_name == "welfare") {
  # Load raw data
  df <- readr::read_csv(file = "https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Welfare/ProcessedData/welfarenolabel3.csv", na = character())
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "y"
  treatment_variable_name <- "w"
  covariate_names <- c("hrs1", "partyid", "income", "rincome", "wrkstat", "wrkslf","age", "polviews",
                       "educ", "earnrs", "race","wrkslf",
                       "marital","sibs","childs", "occ80",  "prestg80", "indus80","res16","reg16","mobile16", "family16", "parborn","maeduc","degree","sex","race","born","hompop","babies","preteen","teens","adults")
  
  # This dataset encodes missing values as -999. Let's change these to NA.
  df[df == -999] <- NA
} else if (dataset_name == "advertising") {
  # Load raw data
  current_directory <- getwd()
  #  checking to see if the criteo data already exists in the directory. if not download the data
  if (!("adcontentworth_qje.csv" %in% list.files(current_directory))){
    url <- 'https://raw.githubusercontent.com/gsbDBI/ExperimentData/master/Advertising/adcontentworth_qje.csv'
    download.file(url, paste0(current_directory, "/", "adcontentworth_qje.csv"))
  }
  # Load raw data
  df <- read.csv(file = "adcontentworth_qje.csv")
  # Subsetting based on wave>1 yields the 53,194 observations used in paper
  df <- df[df$wave > 1, ]
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "applied" # alternatively, use "tookup"
  treatment_variable_name <- "speak_trt" # alternatively, use "oneln_trt" or another
  covariate_names <- c("offer4", "stripany",  "dphoto_black", "dphoto_female", "prize", "oneln_trt", "use_any", "intshown", "comploss_n", "comp_n", "risk", "waved3", "female", "race", "race_match" , "gender_match")
  # --- Content treatments from paper ---
  # speak_trt, stripany, dphoto_none, dphoto_black, dphoto_female,
  # gender_match, race_match, prize, oneln_trt, use_any, intshown, comploss_n ,comp_n
  # -------------------------------------
} else if (dataset_name == "criteo") {
  # Load raw data
  current_directory <- getwd()
  #  checking to see if the criteo data already exists in the directory. if not download the data
  if (!("criteo.csv.gz" %in% list.files(current_directory))){
    url <- 'https://s3.us-east-2.amazonaws.com/criteo-uplift-dataset/criteo-uplift.csv.gz'
    download.file(url, paste0(current_directory, "/", "criteo.csv.gz"))
  }
  df <- readr::read_csv("criteo.csv.gz", na = character())
  # Specify outcome, treatment, and covariate variable names to use
  outcome_variable_name <- "visit"
  treatment_variable_name <- "treatment"
  covariate_names <- c( "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "f10", "f11" )
} else {
  # Raise an error if the dataset is unknown
  print("Incorrect dataset is specified. Change 'dataset_name' to one of {charitable, mobilization, secrecy, social, welfare, criteo }")
  # Exiting knitr
  knitr::knit_exit()
  print("Exiting knitr without completion")
}

#### SETUP TEXT ####
#' # Part One  
#' This problem set examines the `r dataset_name` data set. Throughout, the treatment variable will be referred to as `r toupper(treatment_variable_name)` and the outcome variable will be referred to as `r toupper(outcome_variable_name)`.   
#' Note that models with the suffix `.int` refer to our interacted data; we work with this expanded dataset to demonstrate the usefulness of machine learning methods in higher dimensions.  

#' ## Collaborators
#' I worked closely on this problem set with Ayush Kanodia. I also worked with Kaleb Javier and Haviland Sheldahl-Thomason, and indicate where we collaborated on code.   

#### DATA WORK ####
# follows tutorial exactly
#' ## Pre-Processing the Data
#' I use code from a set of AtheyLab tutorials, which include the following note:  
#'        > The datasets in our [github webpage](https://github.com/gsbDBI/ExperimentData/) have been prepared for analysis so they will not require a lot of cleaning and manipulation, but let's do some minimal housekeeping. First, we will drop the columns that aren't outcomes, treatments or (pre-treatment) covariates, since we won't be using those.
#' Specifically, we keep only a subset of predictors and drop observations with missing information.  
all_variables_names <- c(outcome_variable_name, treatment_variable_name, covariate_names)
df <- df[, which(names(df) %in% all_variables_names)]
# Drop rows containing missing values
df <- na.omit(df)
# Rename variables
names(df)[names(df) == outcome_variable_name] <- "Y"
names(df)[names(df) == treatment_variable_name] <- "W"
# Converting all columns to numerical
df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))))

# coerce to data.table for future analysis
setDT(df)

# if option supplied to randomly subset, do so
df <- df[base::sample(1:nrow(df), prop_to_keep * nrow(df))]


#### RCT ANALYSIS ####
#' ## RCT Analysis
#' We now report the (presumably true) treatment effect $\hat{\tau}$ from the randomized experiment:  
#+ echo=TRUE
tauhat_rct <- difference_in_means(df)
print(tauhat_rct)

#### SAMPLING BIAS ####
#' ## Introducing Bias
#+ echo=FALSE
# these cutoffs courtesy of Kaleb (poc, .45 and .01)
prob = .4
prob_control = prob # how much of the sample are we dropping?
prob_treatment = prob
# prob_control = 0.45  # how much of the sample are we dropping?
# prob_treatment = 0.01

#' We now introduce sampling bias in order to simulate the situation we would be in if our data was from an observational study rather than a randomized experiment. 
#' This situation might arise due to sampling error or selection bias, and we will be able to see how various methods correct for this induced bias. 
#' To do, so, we under-sample treated units matching the following rule, and under-sample control units in its complement:   
#'           - Independents on closer to the Democratic party (`partyid < 4`)
#'           - Who have at least a college degree (`educ >= 16`)
#'           
# - Republican-leaning independent or closer to the Republican party (`partyid >= 4 & partyid < 7`)
# - More than slightly conservative (`polviews >= 5.5 & polviews < 8`)
# - with self-reported prestige of occuption above the bottom quartile (`prestg80 >= 34`)
#' We remove `r prob * 100` percent of observations in these sets. 
#' We picked our rule by using a tree to pick covariates highly correlated with the outcome, so that we could drop high outcomes and low outcome groups from treatment and control to confound results.  
#+ include=FALSE
setDT(df)
# college_educ_lib_dems = df[, partyid < 3 & reg16 <= 16 & degree <=1] #  & df$educ >= 16 & df$polviews < 3.5
# how to find a good sample? use a forest for splitting!

# library(rpart)
# library(rpart.plot)
# mypart <- rpart(Y ~ . - W, data = df, cp=0.001)
# rpart.plot::rpart.plot(mypart)

# alternative (better), courtesy of Kaleb:
# poc =   df$race != 1

# obs_to_remove = df[, partyid >= 4 & partyid < 7 & polviews >= 5.5 & polviews < 8 & prestg80 >= 34]
# Our method:
obs_to_remove = df[, partyid < 4 & educ >= 16]
#obs_to_remove = poc
# obs_to_remove = df[,((persons > 1) | (newreg > 0)) | ((competiv > 1)) | (vote00 == 1)]
mean(df$Y[obs_to_remove])
mean(df$Y[!obs_to_remove])


drop_from_treat <- base::sample(which(obs_to_remove & df$W == 1), round(prob_treatment * sum(obs_to_remove)))
drop_from_control <- base::sample(which(!obs_to_remove & df$W == 0), round(prob_control * sum(!obs_to_remove)))

df_mod <- df[-c(drop_from_treat, drop_from_control),]

#' The difference in means is now biased, and significantly outside the confidence interval indicated by the RCT.  
#' Check if difference in treatment effect estimates is substantial
tauhat_confounded <- difference_in_means(df_mod)
tauhat_confounded

Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
Ymod = df_mod$Y
Wmod = df_mod$W
XWmod = cbind(Xmod, Wmod)
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")
# linear models
tauhat_logistic_ipw <- ipw(df_mod, pW_logistic)
tauhat_pscore_ols <- prop_score_ols(df_mod, pW_logistic)
tauhat_lin_logistic_aipw <- aipw_ols(df_mod, pW_logistic)

# tauhat_logistic_ipw 
# tauhat_pscore_ols 
# tauhat_lin_logistic_aipw 

#### Q1 SIMULATION ####
#' # Problem 1
#' ## CATE Simulation
#' We first run the following simulation:
#+ echo=TRUE
n = 8000; p = 6
taufn = function(x) { 1 / 3 }
X = matrix(rnorm(n*p), n, p)
tau = apply(X, 1, taufn)
W = rbinom(n, 1, 1/(1 + exp(-X[,1] / 2)))
Y = pmax(0, X[,1] / 2 + X[,2]) + W * tau + rnorm(n)
df <- cbind(X, W, Y)
colnames(df) <- c(paste0("X", 1:ncol(X)), "W", "Y")
df <- as.data.frame(df)


#' We now train a causal forest, compate out-of-bag CATE estimates, and sort the simulation by estimated CATE into quartiles.  
#+ echo=TRUE
cf <- causal_forest(
  X = as.matrix(df[,grepl("X", colnames(df))]),
  Y = df$Y,
  W = df$W,
  num.trees=200)

oob_pred <- predict(cf, estimate.variance=TRUE)
df$cate <- oob_pred$predictions

# Manually creating subgroups
num_tiles <- 4  # ntiles = CATE is above / below the median
df$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))

# true ATE per quartile?
df$tau <- tau

# ATE estimates
df_qtile <- setDT(df)[,.(
  avg_tau = mean(tau),
  avg_cf_cate = mean(cate)
  ),.(ntile)][order(ntile)]

df_qtile

#' We see that the causal forest CATE estimates within each quartile are substantially different than the actual effect.  
#' Further, we see that there is substantial correlation between the quartile number and the average CATE estimate:
#+ echo=TRUE
cor(as.numeric(df_qtile$ntile), df_qtile$avg_cf_cate)

#' We now follow the HTE tutorial in using AIPW to estimate the CATE in each quartile. 
#+ echo=TRUE

estimated_aipw_ate <- lapply(
  seq(num_tiles), function(w) {
    ate <- average_treatment_effect(cf, subset = df$ntile == w)
  })
estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
estimated_aipw_ate$ntile <- factor(seq(num_tiles))
setDT(estimated_aipw_ate)
setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
# join CATE estimates together
df_qtile <- df_qtile[estimated_aipw_ate, on = .(ntile)]
df_qtile

#' We see that AIPW is generally much closer to the actual $\tau$, and that even when estimated within each quartile, $\tau$ lies within the confidence intervals implied by the standard error. 
#' Further, we see that the correlation between the AIPW estimate and the quartile number is lower than that of the causal forest estimate:  
#+ echo=TRUE
cor(as.numeric(df_qtile$ntile), df_qtile$aipw_estimate)

#' We now repeat the same analysis for the same experimental setup an additional `r n_sims` times. 
#+ echo=TRUE
taufn = function(x) { 1 / 3 }

# save things across sims
cf_cor_storage <- data.table()
df_qtile_storage <- data.table()
aipw_cor_storage <- data.table()

for (i in 1:n_sims){
  n = 8000; p = 6
  X = matrix(rnorm(n*p), n, p)
  tau = apply(X, 1, taufn)
  W = rbinom(n, 1, 1/(1 + exp(-X[,1] / 2)))
  Y = pmax(0, X[,1] / 2 + X[,2]) + W * tau + rnorm(n)
  df <- cbind(X, W, Y)
  colnames(df) <- c(paste0("X", 1:ncol(X)), "W", "Y")
  df <- as.data.frame(df)
  
  cf <- causal_forest(
    X = as.matrix(df[,grepl("X", colnames(df))]),
    Y = df$Y,
    W = df$W,
    num.trees=200)
  
  oob_pred <- predict(cf, estimate.variance=TRUE)
  df$cate <- oob_pred$predictions
  
  # Manually creating subgroups
  num_tiles <- 4  # ntiles = CATE is above / below the median
  df$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))
  
  # true ATE per quartile?
  df$tau <- tau
  
  # ATE estimates
  df_qtile <- setDT(df)[,.(
    avg_tau = mean(tau),
    avg_cf_cate = mean(cate)
  ),.(ntile)][order(ntile)]
  
  cf_cor <- cor(as.numeric(df_qtile$ntile), df_qtile$avg_cf_cate) %>% as.data.table()
  setnames(cf_cor, "cor")
  
  cf_cor_storage <- rbindlist(list(cf_cor_storage, cf_cor))
  
  estimated_aipw_ate <- lapply(
    seq(num_tiles), function(w) {
      ate <- average_treatment_effect(cf, subset = df$ntile == w)
    })
  estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
  estimated_aipw_ate$ntile <- factor(seq(num_tiles))
  setDT(estimated_aipw_ate)
  setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
  # join CATE estimates together
  df_qtile <- df_qtile[estimated_aipw_ate, on = .(ntile)]
  df_qtile$sim <- i
  
  df_qtile_storage <- rbindlist(list(df_qtile_storage, df_qtile))
  
  aipw_cor <- cor(as.numeric(df_qtile$ntile), df_qtile$aipw_estimate) %>% as.data.table()
  setnames(aipw_cor, "cor")
  aipw_cor_storage <- rbindlist(list(aipw_cor_storage, aipw_cor))
}

#' We now summarize the results over our simulations.  
#' We see that causal forests has a high degree of correlation between CATE and `ntile` over all simulations:  
#+ echo=TRUE
cf_cor_storage %>% summary()

#' On the other hand, the AIPW correlations are centered close to zero:  
#+ echo=TRUE 
aipw_cor_storage %>% summary()

#' Finally, we summarize the means over each variable over all simulations by `ntile`. 
#' We see again see the high degree of correlation between `ntile` and `avg_cf_cate`, and that `aipw_estimate` are close to the true $\tau$. 
#+ echo=TRUE
df_qtile_storage[,lapply(.SD, mean), .(ntile)]

#' ## Simulation Exercise with new $\tau$
#' In this section we repeat the analysis above, but we now use:
#+ echo=TRUE
taufn = function(x) { 1 / (1 + exp(-x[2]/2)) }

#+ echo=FALSE

# save things across sims
cf_cor_storage <- data.table()
df_qtile_storage <- data.table()
aipw_cor_storage <- data.table()

for (i in 1:n_sims){
  n = 8000; p = 6
  X = matrix(rnorm(n*p), n, p)
  tau = apply(X, 1, taufn)
  W = rbinom(n, 1, 1/(1 + exp(-X[,1] / 2)))
  Y = pmax(0, X[,1] / 2 + X[,2]) + W * tau + rnorm(n)
  df <- cbind(X, W, Y)
  colnames(df) <- c(paste0("X", 1:ncol(X)), "W", "Y")
  df <- as.data.frame(df)
  
  cf <- causal_forest(
    X = as.matrix(df[,grepl("X", colnames(df))]),
    Y = df$Y,
    W = df$W,
    num.trees=200)
  
  oob_pred <- predict(cf, estimate.variance=TRUE)
  df$cate <- oob_pred$predictions
  
  # Manually creating subgroups
  num_tiles <- 4  # ntiles = CATE is above / below the median
  df$ntile <- factor(ntile(oob_pred$predictions, n=num_tiles))
  
  # true ATE per quartile?
  df$tau <- tau
  
  # ATE estimates
  df_qtile <- setDT(df)[,.(
    avg_tau = mean(tau),
    avg_cf_cate = mean(cate)
  ),.(ntile)][order(ntile)]
  
  cf_cor <- cor(as.numeric(df_qtile$ntile), df_qtile$avg_cf_cate) %>% as.data.table()
  setnames(cf_cor, "cor")
  
  cf_cor_storage <- rbindlist(list(cf_cor_storage, cf_cor))
  
  estimated_aipw_ate <- lapply(
    seq(num_tiles), function(w) {
      ate <- average_treatment_effect(cf, subset = df$ntile == w)
    })
  estimated_aipw_ate <- data.frame(do.call(rbind, estimated_aipw_ate))
  estimated_aipw_ate$ntile <- factor(seq(num_tiles))
  setDT(estimated_aipw_ate)
  setnames(estimated_aipw_ate, c("aipw_estimate", "aipw_std.err", "ntile"))
  # join CATE estimates together
  df_qtile <- df_qtile[estimated_aipw_ate, on = .(ntile)]
  df_qtile$sim <- i
  
  df_qtile_storage <- rbindlist(list(df_qtile_storage, df_qtile))
  
  aipw_cor <- cor(as.numeric(df_qtile$ntile), df_qtile$aipw_estimate) %>% as.data.table()
  setnames(aipw_cor, "cor")
  aipw_cor_storage <- rbindlist(list(aipw_cor_storage, aipw_cor))
}

#' We now summarize the results over our simulations.  
#' We again see that causal forests has a high degree of correlation between CATE and `ntile` over all simulations:  
#+ echo=TRUE
cf_cor_storage %>% summary()

#' This time, AIPW estimates have a similar degree of correlation:
#+ echo=TRUE 
aipw_cor_storage %>% summary()

#' Again, we summarize the means over each variable over all simulations by `ntile`. 
#' This time we see a high degree of correlation between `ntile`, `avg_cf_cate`, and  `aipw_estimate`, 
#' and that our estimates of $\tau$ are badly biased, and outside our confidence intervals for our AIPW estimates.  
#+ echo=TRUE
df_qtile_storage[,lapply(.SD, mean), .(ntile)]

#' ## Testing Calibration
#' For this section, we run a single simulation with two different $\tau$ functions, but otherwise identical data so that our results are comparable. 
#' 
#+ echo=TRUE

taufn1 = function(x) { 1 / 3 }
taufn2 = function(x) { 1 / (1 + exp(-x[2]/2)) }

X = matrix(rnorm(n*p), n, p)
tau1 = apply(X, 1, taufn1)
tau2 = apply(X, 1, taufn2)
W = rbinom(n, 1, 1/(1 + exp(-X[,1] / 2)))
epsilon = rnorm(n)
Y1 = pmax(0, X[,1] / 2 + X[,2]) + W * tau1 + epsilon
Y2 = pmax(0, X[,1] / 2 + X[,2]) + W * tau2 + epsilon
df1 <- cbind(X, W, Y1)
df2 <- cbind(X, W, Y2)
colnames(df1) <- c(paste0("X", 1:ncol(X)), "W", "Y")
colnames(df2) <- c(paste0("X", 1:ncol(X)), "W", "Y")
df1 <- as.data.frame(df1)
df2 <- as.data.frame(df2)

cf1 <- causal_forest(
  X = as.matrix(df1[,grepl("X", colnames(df1))]),
  Y = df1$Y,
  W = df1$W,
  num.trees=200)

cf2 <- causal_forest(
  X = as.matrix(df2[,grepl("X", colnames(df2))]),
  Y = df2$Y,
  W = df2$W,
  num.trees=200)

tc1 <- test_calibration(cf1)
tc2 <- test_calibration(cf2)
# convert to one-sided p-values
#dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
#blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 4] / 2, blp.summary[, 4] / 2)
# convert to one-sided p-values
#dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
#blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[, 4] / 2, blp.summary[, 4] / 2)

tc1
tc2

#' We see that the $\alpha$ is close to 1 under both $\tau$, so in either case the average prediction from the forest is correct. 
#' From the tutorial we see that we should not interpret the coefficient on $\beta$ in the second case as it is negative, but in the first we capture some of the heterogeneity in the data. 

#' ### Recommended Best Practices
#' 

#### Q2 Heterogeneous Treatment Effects in Observational Studies ####
#' # Heterogeneous Treatment Effects in Observational Studies 
#' 
#' ## Different random forest based strategies for estimation heterogeneous treatment effects
#' # For this part, we include the artificial confounding from the first problem set in `df_mod`
#' ### S-learner
#' We first examine the S-learner, which is a single random forest, fitted to all of the data. 
#' We estimate $\hat\tau = $\hat\mu((x,1)) - \hat\mu((x,0))$.  
#+ echo=TRUE

# from lecture

s_learner <- function(df_mod){
  Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
  Ymod = df_mod$Y
  Wmod = df_mod$W
  XWmod = cbind(Xmod, Wmod)
  sf = regression_forest(cbind(Xmod, Wmod), Ymod)
  pred.sf.0 = predict(sf, cbind(Xmod, 0))$predictions
  pred.sf.1 = predict(sf, cbind(Xmod, 1))$predictions
  preds.sf.oob = predict(sf)$predictions
  pred.sf.0[Wmod==0] = preds.sf.oob[Wmod==0]
  pred.sf.1[Wmod==1] = preds.sf.oob[Wmod==1]
  
  preds.sf = pred.sf.1 - pred.sf.0
  # mean(preds.sf, na.rm = TRUE)
  preds.sf
}

# s_learner(df_mod)



#' ### T-learner
#' We now examine the T-learner, where we estimate two random forests
#' We follow the example code from lecture:
#+ echo=TRUE
t_learner <- function(df_mod){
  Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
  Ymod = df_mod$Y
  Wmod = df_mod$W
  XWmod = cbind(Xmod, Wmod)
  tf0 = regression_forest(Xmod[Wmod==0,], Ymod[Wmod==0])
  tf1 = regression_forest(Xmod[Wmod==1,], Ymod[Wmod==1])
  tf.preds.0 = predict(tf0, Xmod)$predictions
  tf.preds.1 = predict(tf1, Xmod)$predictions
  tf.preds.0[Wmod==0] = predict(tf0)$predictions  #OOB
  tf.preds.1[Wmod==1] = predict(tf1)$predictions  #OOB
  preds.tf = tf.preds.1 - tf.preds.0
  # mean(preds.tf, na.rm = TRUE)
  preds.tf
}

# t_learner(df_mod)

#' ### X-learner
#' We now estimate $\hat\tau$ using the X-learner, which
#+ echo=TRUE
x_learner <- function(df_mod){
  Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
  Ymod = df_mod$Y
  Wmod = df_mod$W
  XWmod = cbind(Xmod, Wmod)
  tf0 = regression_forest(Xmod[Wmod==0,], Ymod[Wmod==0])
  yhat0 = predict(tf0, Xmod[Wmod==1,])$predictions
  xf1 = regression_forest(Xmod[Wmod==1,], Ymod[Wmod==1]-yhat0)
  xf.preds.1 = predict(xf1, Xmod)$predictions
  xf.preds.1[Wmod==1] = predict(xf1)$predictions
  tf1 = regression_forest(Xmod[Wmod==1,], Ymod[Wmod==1])
  yhat1 = predict(tf1, Xmod[Wmod==0,])$predictions
  xf0 = regression_forest(Xmod[Wmod==0,], yhat1-Ymod[Wmod==0])
  xf.preds.0 = predict(xf0, Xmod)$predictions
  xf.preds.0[Wmod==0] = predict(xf0)$predictions
  propf = regression_forest(Xmod, Wmod, tune.parameters = "all")
  ehat = predict(propf)$predictions
  preds.xf = (1 - ehat) * xf.preds.1 + ehat * xf.preds.0
  preds.xf
}

# x_learner(df_mod)

#' ### Causal Forest
#' Finally, we estimate $\hat\tau$ with a causal forest. 
#+ echo=TRUE

cf_learner <- function(df_mod){
  Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
  Ymod = df_mod$Y
  Wmod = df_mod$W
  XWmod = cbind(Xmod, Wmod)
  cf = causal_forest(Xmod, Ymod, Wmod, tune.parameters = "all")
  preds.cf <- predict(cf)
  preds.cf$predictions
  # ate.hat = average_treatment_effect(cf,target.sample = "all")
  # # paste("95% CI:", as.character(round(ate.hat["estimate"], 3)),"+/-",  as.character(round(1.96 * ate.hat["std.err"], 3)))
  # as.numeric(ate.hat["estimate"])
}

# cf_learner(df_mod)

#' Here we aggregate and compare our results. 
#+ df_mod, echo=TRUE
p <- mean(df_mod$W)
Y_star <- ((df_mod$W - p)/(p*(1-p)))*df_mod$Y

# Compute test mse for all methods
mse_1 <- data.frame(
  Causal_Forest_Loss = (Y_star - cf_learner(df_mod))^2,
  X_Learner_Loss = (Y_star - x_learner(df_mod))^2,
  S_Learner_Loss = (Y_star - s_learner(df_mod))^2,
  T_Learner_Loss = (Y_star - s_learner(df_mod))^2)
mse_summary_1 <- describe(mse_1)[, c('mean', 'se')]
#+ results='asis'
kable_styling(kable(mse_summary_1,  "html", digits = 5,
                    caption="Estimate loss: comparison across methods"),
              bootstrap_options=c("striped", "hover", "condensed", "responsive"),
              full_width=FALSE)

#+ echo=TRUE
rloss_long <- mse %>% pivot_longer(cols = everything())

ggplot(rloss_long,aes(x=value)) + 
  geom_histogram() + 
  facet_grid(rows  = vars(name))   
#' ## Changing Size of Datasets
#' In this section, we change the size of the data along various dimensions, and compare results across our 4 models.  
#' ### Draw a random subset of 20% of your data (or 400 observations, whichever is greater)
#+ df_mod_20pct, echo=TRUE
df_mod_20pct <- df_mod %>% 
  dplyr::sample_frac(.2) %>% setDT()

p <- mean(df_mod_20pct$W)
Y_star <- ((df_mod_20pct$W - p)/(p*(1-p)))*df_mod_20pct$Y
mse_2 <- data.frame(
  Causal_Forest_Loss = (Y_star - cf_learner(df_mod_20pct))^2,
  X_Learner_Loss = (Y_star - x_learner(df_mod_20pct))^2,
  S_Learner_Loss = (Y_star - s_learner(df_mod_20pct))^2,
  T_Learner_Loss = (Y_star - s_learner(df_mod_20pct))^2)
mse_summary_2 <- describe(mse_2)[, c('mean', 'se')]
#+ results='asis'
kable_styling(kable(mse_summary_2,  "html", digits = 5,
                    caption="Estimate loss: comparison across methods"),
              bootstrap_options=c("striped", "hover", "condensed", "responsive"),
              full_width=FALSE)
#+ echo=TRUE
rloss_long <- mse %>% pivot_longer(cols = everything())

ggplot(rloss_long,aes(x=value)) + 
  geom_histogram() + 
  facet_grid(rows  = vars(name))   

#' ### Subset the data such that there are exactly the same number of treated and control units
#+ df_mod_ctrl_treat_balance, echo=TRUE
# Calculate number of treated and control units
n_ctrl <- df_mod[W==0,.N]
n_treat <- df_mod[W==1,.N]
# Sample 
df_mod_ctrl <- df_mod[W==0] %>% sample_n(min(n_ctrl, n_treat)) %>% setDT()
df_mod_treat <- df_mod[W==1] %>% sample_n(min(n_ctrl, n_treat)) %>% setDT()
df_mod_ctrl_treat_balance <- rbindlist(list(df_mod_ctrl, df_mod_treat))

p <- mean(df_mod_ctrl_treat_balance$W)
Y_star <- ((df_mod_ctrl_treat_balance$W - p)/(p*(1-p)))*df_mod_ctrl_treat_balance$Y
mse_3 <- data.frame(
  Causal_Forest_Loss = (Y_star - cf_learner(df_mod_ctrl_treat_balance))^2,
  X_Learner_Loss = (Y_star - x_learner(df_mod_ctrl_treat_balance))^2,
  S_Learner_Loss = (Y_star - s_learner(df_mod_ctrl_treat_balance))^2,
  T_Learner_Loss = (Y_star - s_learner(df_mod_ctrl_treat_balance))^2)
mse_summary_3 <- describe(mse_3)[, c('mean', 'se')]
#+ results='asis'
kable_styling(kable(mse_summary_3,  "html", digits = 5,
                    caption="Estimate loss: comparison across methods"),
              bootstrap_options=c("striped", "hover", "condensed", "responsive"),
              full_width=FALSE)
#+ echo=TRUE
rloss_long <- mse %>% pivot_longer(cols = everything())

ggplot(rloss_long,aes(x=value)) + 
  geom_histogram() + 
  facet_grid(rows  = vars(name))   
#' ### Subset the data such that there are 5x more control units than treated units 
#+ df_mod_more_ctrl, echo=TRUE

df_mod_ctrl <- df_mod[W==0]
df_mod_treat <- df_mod[W==1] %>% sample_n(n_ctrl / 5) %>% setDT()
df_mod_more_ctrl <- rbindlist(list(df_mod_ctrl, df_mod_treat))

p <- mean(df_mod_more_ctrl$W)
Y_star <- ((df_mod_more_ctrl$W - p)/(p*(1-p)))*df_mod_more_ctrl$Y
mse_4 <- data.frame(
  Causal_Forest_Loss = (Y_star - cf_learner(df_mod_more_ctrl))^2,
  X_Learner_Loss = (Y_star - x_learner(df_mod_more_ctrl))^2,
  S_Learner_Loss = (Y_star - s_learner(df_mod_more_ctrl))^2,
  T_Learner_Loss = (Y_star - s_learner(df_mod_more_ctrl))^2)
mse_summary_4 <- describe(mse)[, c('mean', 'se')]
#+ results='asis'
kable_styling(kable(mse_summary_4,  "html", digits = 5,
                    caption="Estimate loss: comparison across methods"),
              bootstrap_options=c("striped", "hover", "condensed", "responsive"),
              full_width=FALSE)

#+ echo=TRUE
rloss_long <- mse %>% pivot_longer(cols = everything())

ggplot(rloss_long,aes(x=value)) + 
  geom_histogram() + 
  facet_grid(rows  = vars(name))   

#### QUESTION 3 ####
#' ### Question 3: Heterogeneous Treatment Effects in Randomized Experiments
#' First split data into 
#+ echo=TRUE
library(causalTree)
library(gt)

df <- df %>% mutate(id = row_number())
#Create training set
df_tr <- df %>% sample_frac(.40)
df_split  <- anti_join(df, df_tr, by = 'id')
df_est <- df_split %>% sample_frac(.50)
df_test <- anti_join(df_split, df_est, by = 'id')
#Create test set
nrow(df)
nrow(df_tr) + nrow(df_est) + nrow(df_test)

fmla_ct <- paste("factor(Y) ~", paste(covariate_names, collapse = " + "))

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
ols_ct_test <- lm_robust(Y ~ 0 + leaf + W:leaf, data=df_est)
ols_ct_tr <- lm_robust(Y ~ 0 + leaf + W:leaf, data=df_est)
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
             V4 = "df_tr")

#' ### Partial Dependence Plots


#### Partial dependence plots

#' Directly from tutorial:
#' 
#' It may also be interesting to examine how our CATE estimates behave when we change a single covariate, while keeping all the other covariates at a some fixed value. In the plot below we evaluate a variable of interest across quantiles, while keeping all other covariates at their median (see the RMarkdown source for code).  


#' **Note:** It is important to recognize that in the following plots and tables, we may be evaluating the CATE at $x$ values in regions where there are few or no data points. Also, it may be the case that varying some particular variable while keeping others fixed may just not be very interesting. For example, in the _welfare_ dataset, we will not see a lot difference when we change `partyid` if we keep `polviews` fixed at their median value. It might be instructive to re-run this tutorial without using the variable `partyid`.



#### 90% ####
#' First 90% of the data
df_90 <- df %>% sample_frac(.90)
df_tr <- df_90 %>% sample_frac(.40)
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

#' Now we plot
#+
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


df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

#' Now for another variable:
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

#' Plot
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

df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

#### 70% ####
#' Now 70% of the data
df_70 <- df %>% sample_frac(.70)
df_tr <- df_70 %>% sample_frac(.40)
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

#' Now we plot
#+
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


df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

#' Now for another variable:
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

#' Plot
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

df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

#### 50% ####
#' Finally 50% of the data
df_50 <- df %>% sample_frac(.50)
df_tr <- df_50 %>% sample_frac(.40)
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

#' Now we plot
#+
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


df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

#' Now for another variable:
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

#' Plot
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

df_eval[,c(var_of_interest,"tauhat","se")] %>% gt()

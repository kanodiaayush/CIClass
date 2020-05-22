#' ---
#' title: "Problem Set 1"
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
# set global options
dataset_name <- "welfare"
outcome_family <- "binomial" # based on whether your outcome is binary or not; input to glm call
outcome_type <- "class"
n_sims <- 20
lambda <- c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5) #, 0.7, 1, 5, 10, 50, 100, 1000)
#lambda <- c(0.0001, 0.01)
#lambda <- c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1, 5, 10, 50, 100, 1000)
prop_to_keep <- 1.0 # if you want to only run on a random sample of the data, if want to run on full data set to 1.0

prop_drop_rf <- c(0.01, 0.02, 0.04, 0.1, 0.2, 0.3, 0.4, .5)

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

#### KNITR SETUP ####
#+ setup, include=FALSE
knitr::opts_chunk$set(
  eval= TRUE,
  echo = TRUE,
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

aipw_both <- function(dataset, p, outcome_model) {
  
#  ols.fit = lm(Y ~ W * ., data = dataset)
  
  dataset.treatall = dataset
  dataset.treatall$W = 1
  treated_pred = predict(outcome_model, dataset.treatall)
  
  dataset.treatnone = dataset
  dataset.treatnone$W = 0
  control_pred = predict(outcome_model, dataset.treatnone)
  
  actual_pred = predict(outcome_model, dataset)
  
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
#'
#' ## Collaborators
#' I worked closely on this problem set with Mitchell Linegar, who is my team mate for the course. We also worked with Kaleb Javier Pena and Haviland Sheldahl-Thomason, and indicate where we collaborated on code.   

#### DATA WORK ####
# follows tutorial exactly
#' ## Pre-Processing the Data
#' We use code from a set of AtheyLab tutorials, which include the following note:  
#' > The datasets in our [github webpage](https://github.com/gsbDBI/ExperimentData/) have been prepared for analysis so they will not require a lot of cleaning and manipulation, but let's do some minimal housekeeping. First, we will drop the columns that aren't outcomes, treatments or (pre-treatment) covariates, since we won't be using those.
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


#### Outcome And Treatment Summary ####
#' ## Outcome and Treatment Summary
#' Here is a summary of our outcome and treatment variables
#+ echo=FALSE
df_mod = copy(df)
setDT(df_mod)
Xmod = df_mod[,covariate_names, with=FALSE]
Ymod = df_mod$Y
Wmod = df_mod$W

#' \newpage
#' We plot a histogram of treatment propensities and outcome probabilities, since our outcome is binary 0-1. These are learnt using a logistic regression model.
#' Treatment assignment is roughly normally distributed as per this model, with a mean of `r mean(Wmod)`. Outcome is skewed right with a mean of  `r mean(Ymod)`.  
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")
pW_logistic.fit.tidy <-  pW_logistic.fit %>% tidy()
hist(pW_logistic)

pY_logistic.fit <- glm(Ymod ~ as.matrix(Xmod), family = "binomial")
pY_logistic <- predict(pY_logistic.fit, type = "response")
pY_logistic.fit.tidy <-  pY_logistic.fit %>% tidy()
hist(pY_logistic)

df_mod[, `:=`(p_Y = pY_logistic,
              p_W = pW_logistic)]

#' The following plot compares predicted and actual treatment assignment probabilities. We see that the 
#' predicted probability tracks the true probability well.
#+ echo=FALSE
{plot(smooth.spline(pW_logistic, Wmod, df = 4))
  abline(0, 1)}

#### RCT ANALYSIS ####
#' ## RCT Analysis
#' We now report the gold standard treatment effect and its confidence intervals
tauhat_rct <- difference_in_means(df)
print(tauhat_rct)

#### SAMPLING BIAS ####
#' ## Introducing Bias
#' Now we introduce bias, as per the assignment, to introduce confounding, as if our data is from an observational study.
#' Our biasing rule is simple. We analyzed the data using rpart (a tree) to check which covariates affect outcomes the most.
#' We took hints from it to drop higher outcomes from one group, and lower outcomes from the other group.
#' Then, we analyze how various methods correct for this confounding.
#+ echo=FALSE
prob = .4
prob_control = prob # how much of the sample are we dropping?
prob_treatment = prob
# prob_control = 0.45  # how much of the sample are we dropping?
# prob_treatment = 0.01

#' We under-sample treated units matching the following rule, and under-sample control units in its complement:  
#'         - Independents on closer to the Democratic party (`partyid < 4`)
#'         - Who have at least a college degree (`educ >= 16`)
      # - Republican-leaning independent or closer to the Republican party (`partyid >= 4 & partyid < 7`)
      # - More than slightly conservative (`polviews >= 5.5 & polviews < 8`)
      # - with self-reported prestige of occuption above the bottom quartile (`prestg80 >= 34`)
#' We remove `r prob * 100` percent of observations in these sets. 
#+ include=FALSE
setDT(df)
# college_educ_lib_dems = df[, partyid < 3 & reg16 <= 16 & degree <=1] #  & df$educ >= 16 & df$polviews < 3.5
# how to find a good sample? use a forest for splitting!

# library(rpart)
# library(rpart.plot)
# mypart <- rpart(Y ~ . - W, data = df, cp=0.001)
# rpart.plot::rpart.plot(mypart)

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

#' We report the new difference in means estimator and note that is signficantly biased. The gold standard difference in means
#' also lies outside the confidence intervals that are reported
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
# tauhat_lin_logistic_both <- aipw_both(df_mod, pW_logistic, pY_logistic.fit)

# tauhat_logistic_ipw 
# tauhat_pscore_ols 
# tauhat_lin_logistic_aipw 

#### LOGISTIC PROPENSITY SCORES ####
# df_mod <- copy(df)
Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
Ymod = df_mod$Y
Wmod = df_mod$W
XWmod = cbind(Xmod, Wmod)

# Computing the propensity score by logistic regression of W on X.
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")

# hist(pW_logistic)

df_mod[, logistic_propensity := pW_logistic]
#' We now show summary statistics of the modified dataset

library('stargazer')
stargazer(df_mod)

#### OVERLAP ####
#' We calculate propensity of treatment based on pre treatment covariates. We show that we have significant overlap.
#' There are some propensity scores close to 1. We remove all observations with propensity score > 0.95 to fix this
#' and relearn the propensity model. We show overlap both before and after.
#' Overlap before truncating
overlap <- df_mod %>% ggplot(aes(x=logistic_propensity,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()
overlap

df_mod <- df_mod[logistic_propensity %between% c(0.05, 0.95)]

Xmod = df_mod[,.SD, .SDcols = names(df_mod)[!names(df_mod) %in% c("Y", "W")]] %>% as.matrix()
Ymod = df_mod$Y
Wmod = df_mod$W
XWmod = cbind(Xmod, Wmod)

# Computing the propensity score by logistic regression of W on X.
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")

# hist(pW_logistic)
overlap <- df_mod %>% ggplot(aes(x=logistic_propensity,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()
#' Overlap after truncating
overlap

#### PREDICTING PROPENSITIES AND OUTCOMES, ORIGINAL AND EXPANDED DATA ####
# some of this is used to calculate bias function, hence the ordering
# expand data
Xmod.int = model.matrix(~ . * ., data = as.data.frame(Xmod))
XWmod.int = cbind(Xmod.int, Wmod)
# make expanded df
df_mod.int <- Xmod.int %>% as.data.frame %>% setDT()
df_mod.int[,`:=`(W = Wmod, Y = Ymod)]

# logistic model
# original data
pW_logistic.fit <- glm(Wmod ~ Xmod, family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")
# expanded data
pW_logistic.fit.int <- glm(Wmod ~ Xmod.int, family = "binomial")
pW_logistic.int <- predict(pW_logistic.fit.int, type = "response")

# original data
pY_logistic.fit <- glm(Ymod ~ XWmod, family = outcome_family)
pY_logistic <- predict(pY_logistic.fit, type = "response")
# pY_logistic2 <- predict(pY_logistic.fit, newdata = as.data.frame(XWmod), type = "response")
# expanded data
pY_logistic.fit.int <- glm(Ymod ~ Xmod.int, family = outcome_family)
pY_logistic.int <- predict(pY_logistic.fit.int, type = "response")

# lasso expanded data, code provided by TA
# original data
pW_glmnet.fit.propensity = glmnet::cv.glmnet(Xmod, Wmod, lambda = lambda, family = "binomial", type.measure = "class", keep=TRUE) 
pY_glmnet.fit.propensity = glmnet::cv.glmnet(Xmod, Ymod, lambda = lambda, family = outcome_family, type.measure = outcome_type, keep=TRUE) 
# expanded data
pW_glmnet.fit.propensity.int = glmnet::cv.glmnet(Xmod.int, Wmod, lambda = lambda, family = "binomial", type.measure = "class", keep=TRUE) 
pY_glmnet.fit.propensity.int = glmnet::cv.glmnet(Xmod.int, Ymod, lambda = lambda, family = outcome_family, type.measure = outcome_type, keep=TRUE) 

# demonstration of lasso fit across lambdas:
pW_lasso = pW_glmnet.fit.propensity$fit.preval[, pW_glmnet.fit.propensity$lambda == pW_glmnet.fit.propensity$lambda.min] %>% convert_to_prob()
pW_lasso.min = pW_glmnet.fit.propensity$fit.preval[, pW_glmnet.fit.propensity$lambda == min(pW_glmnet.fit.propensity$lambda)] %>% convert_to_prob()
pW_lasso.max = pW_glmnet.fit.propensity$fit.preval[, pW_glmnet.fit.propensity$lambda == max(pW_glmnet.fit.propensity$lambda)] %>% convert_to_prob()
pW_lasso.rand = pW_glmnet.fit.propensity$fit.preval[, pW_glmnet.fit.propensity$lambda == base::sample(pW_glmnet.fit.propensity$lambda, 1)] %>% convert_to_prob()

# glmnet.fit.propensity = glmnet::cv.glmnet(Xmod.int, Wmod, family = "binomial",keep=TRUE) 
pW_lasso.int = pW_glmnet.fit.propensity.int$fit.preval[, pW_glmnet.fit.propensity.int$lambda == pW_glmnet.fit.propensity.int$lambda.min] %>% convert_to_prob()
pW_lasso.int.min = pW_glmnet.fit.propensity.int$fit.preval[, pW_glmnet.fit.propensity.int$lambda == min(pW_glmnet.fit.propensity.int$lambda)] %>% convert_to_prob()
pW_lasso.int.max = pW_glmnet.fit.propensity.int$fit.preval[, pW_glmnet.fit.propensity.int$lambda == max(pW_glmnet.fit.propensity.int$lambda)] %>% convert_to_prob()
pW_lasso.int.rand = pW_glmnet.fit.propensity.int$fit.preval[, pW_glmnet.fit.propensity.int$lambda == base::sample(pW_glmnet.fit.propensity.int$lambda, 1)] %>% convert_to_prob()

# FIXME: 
# problems calculating probabilities: 
# this "should" work: 
# pW_lasso.int2 <- predict(pW_glmnet.fit.propensity.int, newx = Xmod.int, type = "response")
# however, it returns identical predictions for every entry. This is the same if we specify a lambda, though the predcition changes:
# pW_lasso.int3 <- predict(pW_glmnet.fit.propensity.int, newx = Xmod.int,  s=pW_glmnet.fit.propensity.int$lambda.min, type = "response")
# pW_lasso.int4 <- predict(pW_glmnet.fit.propensity.int, newx = Xmod.int,  s="lambda.min", type = "response")
# pW_lasso.int5 <- predict(pW_glmnet.fit.propensity.int, newx = Xmod.int,  s=pW_glmnet.fit.propensity.int$lambda[3], type = "response")
# pW_lasso.int6 <- predict(pW_glmnet.fit.propensity.int, newx = Xmod.int,  s=pW_glmnet.fit.propensity.int$lambda[1], type = "response")

# random forest
pW_rf.fit = regression_forest(Xmod, Wmod, num.trees = 500)
pY_rf.fit = regression_forest(Xmod, Ymod, num.trees = 500)

# pW_rf = pW_rf.fit$predictions
# pY_rf = pY_rf.fit$predictions
pW_rf = predict(pW_rf.fit, newdata = Xmod) %>% as.matrix
pY_rf = predict(pY_rf.fit, newdata = Xmod) %>% as.matrix

# random forest, expanded data
pW_rf.fit.int = regression_forest(Xmod.int, Wmod, num.trees = 500)
pY_rf.fit.int = regression_forest(Xmod.int, Ymod, num.trees = 500)

# pW_rf.int = pW_rf.fit.int$predictions
# pY_rf.int = pY_rf.fit.int$predictions
pW_rf.int = predict(pW_rf.fit.int, newdata = Xmod.int) %>% as.matrix
pY_rf.int = predict(pY_rf.fit.int, newdata = Xmod.int) %>% as.matrix

# hist(pW_rf)
#### BIAS FUNCTION ####
#' Next we plot the bias function $b(X)$ following Athey, Imbens, Pham and Wager (AER P&P, 2017, Section IIID). 
#+ echo=FALSE

mu_avg <- function(treated, df){df[W==treated, mean(Y)]}
mu <- function(treated, df){df[W==treated, mean(pY)]}

B <- function(df, treatment_model, outcome_model, outcome_type = "response"){
  # have to supply models so that can estimate counterfactual predictions given an alternative treatment assignment
  # note that this will NOT work for lasso model, attempt to warn of this misbehavior:
  if (grepl("lasso|rf|cf", deparse(quote(treatment_model)))){
    simpleMessage("The predict method appears to be broken for lasso models estimated using glmnet::cv.glmnet.
                  You may want to try another predictive model. 
                  I am unclear how predict works using predict.causal_forest, and so currently do not recommend it.")
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
  
  # 1/(p*(1-p)) * df[, mean(b)]
  return(df[,.(b)])
  }

df_mod_bias <- B(df_mod, pW_logistic.fit, pY_logistic.fit)
ggplot(df_mod_bias, aes(x = b)) + geom_histogram() + labs(title = "Histogram of per-observation b(x)")
#' We see skew in the bias, and note that the mean lies at `r mean(df_mod_bias)`

#### ESTIMATING ATE INTRO ####
#' \newpage
#' ## ATE on biased data  
#` We measure ATE and confidence intervals on our biased data using the following methods, to explore their performance.
#' 1. IPW with various propensity models  
#' 2. OLS regression with propensity weights using various propensity models  
#' 3. AIPW (Augmented inverse propensity score), doubly robust methods, using various propensity and outcome models.  

#### PLOTTING PROPENSITY PERFORMANCE ####
#' \newpage
#' Now we analyze each propensity model by plotting how they track the true propensity of treatment. Models close to the Y=X line are better.
#' Logistic and Lasso Propensity scores perform the best in our original dataset. 
#+ echo = FALSE

plot_prob(pW_logistic, Wmod, "Logistic")
plot_prob(pW_lasso, Wmod, "Lasso")
plot_prob(pW_rf, Wmod, "Random Forest")

#' We create an expanded dataset by adding two way interaction effects for all our pre treatment covariates
#' We then add plots for treatment propensity models using this dataset.
#' Once again, Logistic Regression and Lasso perform the best.
plot_prob(pW_logistic.int, Wmod, "Logistic", "Expanded")
plot_prob(pW_lasso.int, Wmod, "Lasso", "Expanded")
plot_prob(pW_rf.int, Wmod, "Random Forest", "Expanded")

#### EXPLORING LASSO ####
#' \newpage
#' ## Lasso Cross Validation
#+ echo = FALSE
#' We plot how the lasso propensity models track the true propensity with various cross validation parameters.
#' Here the lasso is plotted for the larger, interacted, covariate set.
#' The lasso with the best cross validation parameter performs the best.
plot_prob(pW_lasso.int.min, Wmod, "Lasso with Min Lambda", "Expanded")
plot_prob(pW_lasso.int.max, Wmod, "Lasso with Max Lambda", "Expanded")
plot_prob(pW_lasso.int.rand, Wmod, "Lasso with Random Lambda", "Expanded")

# {plot(smooth.spline(pW_logistic, Wmod, df = 4))
#   abline(0, 1)}

# likelihoods over lambdas
lambda_log_liks <- pW_glmnet.fit.propensity.int.lambda_preds[
  ,lapply(.SD, loglike, Wmod), .SDcols = names(pW_glmnet.fit.propensity.int.lambda_preds)]

colnames(lambda_log_liks) <- as.character(pW_glmnet.fit.propensity.int$lambda)
lambda_log_liks.long <- melt(lambda_log_liks, variable.name = "lambda", value.name = "llh")
lambda_log_liks.long[, lambda := as.numeric(as.character(lambda))]

#' We now plot log likelihood over different values of lambda on the expanded data. We see later that the IPW and OLS with propensity score
#' weighting models above perform the best when the log likelihoods of the lasso propensity models are high. This is consistent with our conclusions downstream.
#+ echo=FALSE
ggplot(lambda_log_liks.long, aes(x = lambda, y = llh)) + 
  geom_line() + 
  labs(title = "Log-Likelihood over Lambdas for Predicting W")

ggplot(lambda_log_liks.long, aes(x = lambda, y = llh)) + 
  geom_line() + 
  labs(title = "Log-Likelihood over Lambdas for Predicting W, Zoomed in") + 
  coord_cartesian(x = c(0,1))

#' We show the results from the different methods by using propensity models with different lambdas.
#' We note that we only vary the propensity model; even for the case with AIPW we use ordinary OLS along with our lasso propensity models for varying lambdas.
#' We tried to also very the outcome model; we could not get the predict function to work with our glmnet model; hence the limitation
# plot lasso over grid of lambdas
pW_glmnet.fit.propensity.int.lambda_preds <- as.data.table(pW_glmnet.fit.propensity.int$fit.preval)
pW_glmnet.fit.propensity.int.lambda_preds <- pW_glmnet.fit.propensity.int.lambda_preds[
  # see discussion in FIXME above about using convert_to_prob here
  ,lapply(.SD, convert_to_prob), .SDcols = names(pW_glmnet.fit.propensity.int.lambda_preds)]

# credit to Kaleb for this formulation
 tauhat_lasso_ipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
  data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=ipw(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
}))[, model := "ipw"]
 
 tauhat_lasso_prop_score.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
   data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=prop_score_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
 }))[, model := "prop_score"]
 
 tauhat_lasso_aipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
   data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=aipw_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
 }))[, model := "aipw"]

tauhat_lasso_estimates.lambdas <- rbindlist(list(tauhat_lasso_ipw.lambdas,
                                                tauhat_lasso_prop_score.lambdas,
                                                tauhat_lasso_aipw.lambdas))


ggplot(tauhat_lasso_estimates.lambdas, aes(x = lambda, y = ATE, color = model)) + 
  geom_line() + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"])) + 
  ggtitle("Tauhat estimates over Lambda")

ggplot(tauhat_lasso_estimates.lambdas, aes(x = lambda, y = ATE, color = model)) + 
  geom_line() + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"])) + 
  coord_cartesian(ylim = c(-1, 1)) + 
  ggtitle("Tauhat estimates over Lambdas, zoomed in")

# ggplot(tauhat_lasso_estimates.lambdas, aes(x = lambda, y = ATE, color = model)) + 
#   geom_line() + 
#   geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"])) + 
#   coord_cartesian(ylim = c(tauhat_rct["lower_ci"], tauhat_rct["upper_ci"])) + 
#   ggtitle("Tauhat estimates over Lambdas, zoomed in more")

#' We learn that as the lambdas are varied, the ipw and OLS with propensity score weighting models start getting unstable.
#' This is probably because the model is misspecified with extreme values of lambdas. If the propensity model is misspecified,
#' We expect to learn incorrect treatment effects.
#' However, with the AIPW, we conclude that since the model is "doubly robust", the underlying outcome model (ordinary OLS) is
#' always correctly specified. As a result, we get stable estimates across all lambdas for the propensity lasso model
#' This teaches us the benefits of using the AIPW estimator!

#### EXPLORING CF ####
#' \newpage
#' ## Machine Learning Model of our choice
#' As per the assignment, we now choose a machine learning model of our choice, causal forests, to understand how prediciton quality varies with sample size

# CF
#' First, we use out of bag predictions (cross fitting)
cf.int.fit = causal_forest(Xmod.int, Ymod, Wmod, num.trees = 500)
pW_cf.int = predict(cf.int.fit, newdata = Xmod.int) %>% as.matrix
ate_rf_aipw.int = average_treatment_effect(cf.int.fit)
tauhat_rf_aipw.int = c(ATE=ate_rf_aipw.int["estimate"],
                       lower_ci=ate_rf_aipw.int["estimate"] - 1.96 * ate_rf_aipw.int["std.err"],
                       upper_ci=ate_rf_aipw.int["estimate"] + 1.96 * ate_rf_aipw.int["std.err"])
tauhat_rf_aipw.int
tauhat_rf_ipw.int = ipw(df_mod.int, pW_cf.int)
tauhat_ols_rf_aipw.int = aipw_ols(df_mod.int, pW_cf.int)

tauhat_rf_ipw.int
tauhat_ols_rf_aipw.int


#' Now, we do not use cross fitting.
cf.int.fit = causal_forest(Xmod.int, Ymod, Wmod, num.trees = 500, compute.oob.predictions = FALSE)
pW_cf.int = predict(cf.int.fit, newdata = Xmod.int) %>% as.matrix
ate_rf_aipw.int = average_treatment_effect(cf.int.fit)
tauhat_rf_aipw.int = c(ATE=ate_rf_aipw.int["estimate"],
                       lower_ci=ate_rf_aipw.int["estimate"] - 1.96 * ate_rf_aipw.int["std.err"],
                       upper_ci=ate_rf_aipw.int["estimate"] + 1.96 * ate_rf_aipw.int["std.err"])
tauhat_rf_aipw.int
tauhat_rf_ipw.int = ipw(df_mod.int, pW_cf.int)
tauhat_ols_rf_aipw.int = aipw_ols(df_mod.int, pW_cf.int)

tauhat_rf_ipw.int
tauhat_ols_rf_aipw.int

#' In this setup, we do not see a difference with/without crossfitting

#' Now we vary data sizes with and without cross fitting
#' First, with cross fitting
tauhat_rf_list <- data.frame()
tauhat_ols_rf_aipw_list <- data.frame()

for (prob_temp in prop_drop_rf) {
  # drop same proportion of treated and control units
  drop_from_treat_temp <- base::sample(which(df_mod$W == 1), round(prob_temp * sum(df_mod$W == 1)))
  drop_from_control_temp <- base::sample(which(df_mod$W == 0), round(prob_temp * sum(df_mod$W == 0)))
  df_mod_temp <- df_mod[-c(drop_from_treat_temp, drop_from_control_temp),]
  # df_mod_temp <- copy(df_mod)
  Xmod_temp = df_mod_temp[,.SD, .SDcols = names(df_mod_temp)[!names(df_mod_temp) %in% c("Y", "W")]] %>% as.matrix()
  Ymod_temp = df_mod_temp$Y
  Wmod_temp = df_mod_temp$W
  XWmod_temp = cbind(Xmod_temp, Wmod_temp)
  #pW_rf_temp.fit = causal_forest(Xmod_temp, Wmod_temp, num.trees = 500)
  #pW_rf_temp.fit = regression_forest(Xmod_temp, Wmod_temp, num.trees = 500)
  # pY_rf_mod.fit = regression_forest(Xmod_temp, Ymod_temp, num.trees = 500)
  Xmod_temp.int = model.matrix(~ . * ., data = as.data.frame(Xmod_temp))
  #XWmod_temp.int = cbind(Xmod_temp.int, Wmod_temp)
  df_mod_temp.int <- Xmod_temp.int %>% as.data.frame %>% setDT()
  df_mod_temp.int[,`:=`(W = Wmod_temp, Y = Ymod_temp)]
  # pW_rf_temp.int.fit = regression_forest(Xmod_temp.int, Wmod_temp, num.trees = 500)
  pW_rf_temp.int.fit = causal_forest(Xmod_temp.int, Ymod_temp, Wmod_temp, num.trees = 500)
  # pY_rf_temp.int.fit = regression_forest(Xmod_temp.int, Ymod_temp, num.trees = 500)
  pW_rf_temp.int = predict(pW_rf_temp.int.fit, newdata = Xmod_temp.int) %>% as.matrix
  # pY_rf_temp.int = predict(pY_rf.fit.int, newdata = Xmod_temp.int)
  
  # pW_rf = pW_rf.fit$predictions
  
  tauhat_rf_temp.int = ipw(df_mod_temp.int, pW_rf_temp.int) %>% as.list() %>% data.frame()
  tauhat_rf_temp.int$"prop_dropped" <- prob_temp
  tauhat_rf_temp.int$model <- "rf_ipw"
  tauhat_ols_rf_aipw_temp.int = aipw_ols(df_mod_temp.int, pW_rf_temp.int) %>% as.list() %>% data.frame()
  tauhat_ols_rf_aipw_temp.int$"prop_dropped" <- prob_temp
  tauhat_ols_rf_aipw_temp.int$model <- "rf_aipw"
  
  print(tauhat_rf_temp.int)
  print(tauhat_ols_rf_aipw_temp.int)
  tauhat_rf_list <- rbind(tauhat_rf_list, tauhat_rf_temp.int)
  tauhat_ols_rf_aipw_list <- rbind(tauhat_ols_rf_aipw_list, tauhat_ols_rf_aipw_temp.int)
}

ggplot(tauhat_rf_list, aes(x = prop_dropped, y = ATE, color = model)) + geom_line() + 
  ggtitle("IPW ATE Estimate with Random Forests Given Varied Sample Sizes") + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"]))

ggplot(tauhat_ols_rf_aipw_list, aes(x = prop_dropped, y = ATE, color = model)) + geom_line() + 
  ggtitle("AIPW ATE Estimate with Random Forests Given Varied Sample Sizes") + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"]))

#' Now, without cross fitting
tauhat_rf_list_ncf <- data.frame()
tauhat_ols_rf_aipw_list_ncf <- data.frame()

for (prob_temp in prop_drop_rf) {
  # drop same proportion of treated and control units
  drop_from_treat_temp <- base::sample(which(df_mod$W == 1), round(prob_temp * sum(df_mod$W == 1)))
  drop_from_control_temp <- base::sample(which(df_mod$W == 0), round(prob_temp * sum(df_mod$W == 0)))
  df_mod_temp <- df_mod[-c(drop_from_treat_temp, drop_from_control_temp),]
  # df_mod_temp <- copy(df_mod)
  Xmod_temp = df_mod_temp[,.SD, .SDcols = names(df_mod_temp)[!names(df_mod_temp) %in% c("Y", "W")]] %>% as.matrix()
  Ymod_temp = df_mod_temp$Y
  Wmod_temp = df_mod_temp$W
  XWmod_temp = cbind(Xmod_temp, Wmod_temp)
  #pW_rf_temp.fit = causal_forest(Xmod_temp, Wmod_temp, num.trees = 500)
  #pW_rf_temp.fit = regression_forest(Xmod_temp, Wmod_temp, num.trees = 500)
  # pY_rf_mod.fit = regression_forest(Xmod_temp, Ymod_temp, num.trees = 500)
  Xmod_temp.int = model.matrix(~ . * ., data = as.data.frame(Xmod_temp))
  #XWmod_temp.int = cbind(Xmod_temp.int, Wmod_temp)
  df_mod_temp.int <- Xmod_temp.int %>% as.data.frame %>% setDT()
  df_mod_temp.int[,`:=`(W = Wmod_temp, Y = Ymod_temp)]
  # pW_rf_temp.int.fit = regression_forest(Xmod_temp.int, Wmod_temp, num.trees = 500)
  pW_rf_temp.int.fit = causal_forest(Xmod_temp.int, Ymod_temp, Wmod_temp, num.trees = 500, compute.oob.predictions = FALSE)
  # pY_rf_temp.int.fit = regression_forest(Xmod_temp.int, Ymod_temp, num.trees = 500)
  pW_rf_temp.int = predict(pW_rf_temp.int.fit, newdata = Xmod_temp.int) %>% as.matrix
  # pY_rf_temp.int = predict(pY_rf.fit.int, newdata = Xmod_temp.int)
  
  # pW_rf = pW_rf.fit$predictions
  
  tauhat_rf_temp.int = ipw(df_mod_temp.int, pW_rf_temp.int) %>% as.list() %>% data.frame()
  tauhat_rf_temp.int$"prop_dropped" <- prob_temp
  tauhat_rf_temp.int$model <- "rf_ipw"
  tauhat_ols_rf_aipw_temp.int = aipw_ols(df_mod_temp.int, pW_rf_temp.int) %>% as.list() %>% data.frame()
  tauhat_ols_rf_aipw_temp.int$"prop_dropped" <- prob_temp
  tauhat_ols_rf_aipw_temp.int$model <- "rf_aipw"
  
  print(tauhat_rf_temp.int)
  print(tauhat_ols_rf_aipw_temp.int)
  tauhat_rf_list_ncf <- rbind(tauhat_rf_list_ncf, tauhat_rf_temp.int)
  tauhat_ols_rf_aipw_list_ncf <- rbind(tauhat_ols_rf_aipw_list_ncf, tauhat_ols_rf_aipw_temp.int)
}

#tauhat_rf_list <- rbind(tauhat_rf_list, tauhat_ols_rf_aipw_list)

ggplot(tauhat_rf_list_ncf, aes(x = prop_dropped, y = ATE, color = model)) + geom_line() + 
  ggtitle("IPW ATE Estimate with Random Forests Given Varied Sample Sizes") + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"]))

ggplot(tauhat_ols_rf_aipw_list_ncf, aes(x = prop_dropped, y = ATE, color = model)) + geom_line() + 
  ggtitle("AIPW ATE Estimate with Random Forests Given Varied Sample Sizes") + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"]))

#' We notice that cross fitting improves ATE estimates with increase in size, as the estimates are significantly off and ustable without cross fitting for especially AIPW estimator 
#' The AIPW estimator is almost always better than the IPW estimator in this experiment

cf = causal_forest(Xmod, Ymod, Wmod, num.trees = 500)
cf.int = causal_forest(Xmod.int, Ymod, Wmod, num.trees = 500)

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

#### ATE CALCULATIONS: INTERACTED DATA ####

# linear models
tauhat_ols.int <- ate_condmean_ols(df_mod.int)

# linear models
tauhat_logistic_ipw.int <- ipw(df_mod.int, pW_logistic.int)
tauhat_pscore_ols.int <- prop_score_ols(df_mod.int, pW_logistic.int)
tauhat_lin_logistic_aipw.int <- aipw_ols(df_mod.int, pW_logistic.int)

# lasso
tauhat_lasso_ipw.int <- ipw(df_mod.int, pW_lasso.int)
tauhat_pscore_lasso.int <- prop_score_ols(df_mod.int, pW_lasso.int)
tauhat_lasso_logistic_aipw.int <- aipw_ols(df_mod.int, pW_lasso.int)
# FIXME: add this
# prior code to add:
# Xmod.for.lasso = cbind(Wmod, Xmod.int, (2 * Wmod - 1) * Xmod.int)
# glmnet.fit.outcome = amlinear:::crossfit.cv.glmnet(Xmod.for.lasso, Ymod,
#                                                    penalty.factor = c(0, rep(1, ncol(Xmod.for.lasso) - 1)))
# lasso.yhat.control = amlinear:::crossfit.predict(glmnet.fit.outcome,
#                                                  cbind(0, Xmod.int, -Xmod.int))
# lasso.yhat.treated = amlinear:::crossfit.predict(glmnet.fit.outcome,
#                                                  cbind(1, Xmod.int, Xmod.int))
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
# balancing.weights = amlinear::balance_minimax(Xmod.int, Wmod, zeta = 0.5)
# G.balance = lasso.yhat.treated - lasso.yhat.control +
#   balancing.weights * (Ymod - Wmod * lasso.yhat.treated
#                        - (1 - Wmod) * lasso.yhat.control)
# tau.hat = mean(G.balance)
# se.hat = sqrt(var(G.balance) / length(G.balance))
# tauhat_lasso_balance = c(ATE=tau.hat,
#                          lower_ci=tau.hat-1.96*se.hat,
#                          upper_ci=tau.hat+1.96*se.hat)

#### COMPARING ATE ACROSS MODELS ####
#' ## Comparing ATE Across Models with Original Data
#' We now show ATE estimates and Confidence Intervals across models.
#' AIPW with forest methods provides the most accurate measurements with small confidence intervals.
#' This is true for both the original as well as interacted data
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

#' ## Comparing ATE Across Models with Interacted Data
#+ echo=FALSE
# FIXME: add identical setup for non-interacted models (remove .int, copy code from above and add to own section)
all_estimators.int = rbind(
  RCT_gold_standard = tauhat_rct,
  naive_observational = tauhat_confounded,
  linear_regression = tauhat_ols.int,
  propensity_weighted_regression = tauhat_pscore_ols.int,
  IPW_logistic = tauhat_logistic_ipw.int,
  AIPW_linear_plus_logistic = tauhat_lin_logistic_aipw.int,
  IPW_forest = tauhat_rf_ipw.int,
  AIPW_forest = tauhat_rf_aipw.int,
  AIPW_linear_plus_forest = tauhat_ols_rf_aipw.int,
  IPW_lasso = tauhat_lasso_ipw.int#,
  # FIXME: add this above
  # AIPW_lasso = tauhat_lasso_aipw.int,
  # FIXME: add this above
  # approx_residual_balance = tauhat_lasso_balance
)
all_estimators.int = data.frame(all_estimators.int)
all_estimators.int$ci_length = all_estimators.int$upper_ci - all_estimators.int$lower_ci
print(round(as.matrix(all_estimators.int), 3))


#### PROPENSITY STRATIFICATION INTRODUCTION ####
#' # Part Two
#' 
#' Why does propensity stratification work?
#' We need the unconfoundedness assumption, that is the potential outcomes are uncorrelated with the treatment assignments given the pretreatment covariates
#' In Propensity Score weighting, we deduce that for units with the same propensity score, once we know the propensity score, the unconfoundedness assumption holds.
#' In other words, instead of equalising on the covariates, we can equalise on the propensity score given covariates.
#' The idea is that even if the covariates are different, but the propensity score of treatment is the same, the latter is a complete statistic summarising all the
#' confounding due to covariates. Put differently and simplistically, if there are two different Xs, if their propensity of treatment is the same, then we are again
#' in the completely randomized design setup where the probability of treatment is independent of covariates (assume there are only these two covariates)
#' As such, stratification on propensity scores resembles a high dimensional version of a stratified random experiment, where the strata are formed over higher dimensions
#' using the Propensity Score, learned from some Machine Learning predictive model.
#' Now, as we grow the number of samples, we should improve upon our estimates even with the same number of strata.
#' If we increase the number of strata, we do finer stratification, reducing inter strata differences as in a stratified random design.
#' Hence, we reduce bias, tending to no bias, as we increase the number of strata with more number of points. More points are required for more strata because otherwise,
#' the effects we calculate in each strata might be too noisy with too few points.

#### PROPENSITY STRATIFICATION FUNCTION ####
#' ## Propensity Stratification Function
#' Here we present a function to calculate propensity scores at a strata-level. 
#' The user supplies either a dataframe with treatment column $W$ and outcome column $Y$, 
#' as well as a model with which to estimate propensity scores 
#' (we don't supply a pre-calculated set of propensity scores in case we want to examine the usefulness of propensity stratification for new data). 
#' The function takes options for a number of strata and the function to estimate on the strata - 
#' the user could supply something more complex than the simple difference-in-means, for example.  
#' 
#' What do we do if some strata has only treatment or only control or no units? 
#' We exclude it from our calculations entirely
#' 
<<<<<<< HEAD:assignments/ps1_mlinegar.R
#' Note that the true effect is:
#+ eval=FALSE, include=TRUE 
=======
#' Note that the true effect is   
>>>>>>> 9a80fbe252345f0649cc6ab2f8118c272b18b57a:assignments/ps1_akanodia.R
n = 1000; p = 20
X = matrix(rnorm(n * p), n, p)
Y1 = pmax(X[,1] + X[,2], 0)
Y0 = pmax(X[,1], 0)
mean(Y1 - Y0)
#+ echo=TRUE

propensity_stratification <- function(df, treatment_model, n_strata = 10, 
                                      tau_estimator = difference_in_means){
  df <- copy(df)
  df[, pW := predict(treatment_model, newdata = df[,.SD, .SDcols = !c('W', 'Y')], type = "response")]
  df[, pW_strata := ntile(pW, n_strata)]
  # if any strata is empty, automatically ignored
  # if all
  strata_count <- df[, .N, by = .(pW_strata, W)]
  strata_to_keep <- strata_count[, .(n_strata_W_nonempty = sum(N > 0)), by = .(pW_strata)][
    # keep strata if has observations in both treatment
    n_strata_W_nonempty == 2, pW_strata] %>% unique 
  
  # sloppy but can't figure out right solution atm
  strata_tau_est <- df[pW_strata %in% strata_to_keep, .(
    tauhat = tau_estimator(.SD)[1],
    lower_ci = tau_estimator(.SD)[2],
    upper_ci = tau_estimator(.SD)[3]), by = .(pW_strata)][
      ,.(tauhat = mean(tauhat),
        lower_ci = mean(lower_ci),
        upper_ci = mean(upper_ci))]
  return(c(ATE = strata_tau_est$tauhat, 
           lower_ci = strata_tau_est$lower_ci, 
           upper_ci = strata_tau_est$upper_ci))
}

#### SIMULATION ####
#' ## Simulation Exercise
#' For all discussion that follows, we use the following code to generate a simulated dataset: 
#+ echo=TRUE
make_simulation <- function(){
  n = 1000; p = 20
  X = matrix(rnorm(n * p), n, p)
  propensity = pmax(0.2, pmin(0.8, 0.5 + X[,1]/3)) 
  W = rbinom(n, 1, propensity)
  Y = pmax(X[,1] + W * X[,2], 0) + rnorm(n)
  df = cbind(X, W, Y) %>% as.data.frame() %>% setDT
  df
}

#'
# df <- make_simulation()
# 
# difference_in_means(df)
# 
# tauhat_propensity_stratification <- propensity_stratification(df_mod, pW_logistic.fit)
# tauhat_propensity_stratification
# # to get identical behavior (to be able to supply predictions to newdata for the propensity score model), run this line first:
# # pW_logistic <- predict(pW_logistic.fit, newdata = df_mod, type = "response")
# tauhat_logistic_ipw <- ipw(df_mod, pW_logistic)
# 
# 
# sim_pW_logistic.fit <- glm(W ~ ., data = df %>% select(-Y), family = "binomial")
# sim_pW_logistic <- predict(sim_pW_logistic.fit, type = "response")
# 
# difference_in_means(df)
# ipw(df, sim_pW_logistic)
# propensity_stratification(df, sim_pW_logistic.fit)

<<<<<<< HEAD:assignments/ps1_mlinegar.R
#' We now run `r n_sims` simulations of the type described above for each number of strata, 
#' and report results over each simulation.  
#' We see that propensity stratification and IPW perform similarly, but propensity stratification has lower average MSE.  
#' We settle on 10 strata, as it minimizes bias.  
#+ echo=TRUE
=======
#' We now run `r n_sims` simulations of the type described above, and report results over each simulation.  
#' We also run these simulations for different number of strata
#' We see that propensity stratification and IPW perform similarly, but propensity stratification has lower average MSE.  
#' For 1 stratum, we exactly track the difference in means estimator
#' As we increase the number of strata, we start approaching the ipw estimator. Bias falls.
#' We settle on 10 as the ideal number of strata. The reason is that it provides ATE with the least bias.
#' At a high level, having too many strata increases the amount of noise with a fixed number of points, so we start at
#' 1 and reach an optimum
#+ echo=FALSE

nstratas <- c(1, 2, 5, 10, 20, 50)
for (nstrata in nstratas) {
  sim_storage <- data.table()
  for (i in 1:n_sims){
    df <- make_simulation()
    Y1 = pmax(df$V1 + df$V2, 0)
    Y0 = pmax(df$V1, 0)
    true_tau = mean(Y1 - Y0)
    
    sim_pW_logistic.fit <- glm(W ~ ., data = df %>% select(-Y), family = "binomial")
    sim_pW_logistic <- predict(sim_pW_logistic.fit, type = "response")
    
    tau_ests <- rbind(
      c(true_tau, NA, NA),
      difference_in_means(df),
      ipw(df, sim_pW_logistic),
      propensity_stratification(df, sim_pW_logistic.fit, n_strata=nstrata)  
    ) %>% as.data.frame()
    
    tau_ests$sim_num <- i
    tau_ests$est <- c("true", "difference_in_means", "ipw", "propensity_stratification")
    sim_storage <- rbindlist(list(sim_storage, tau_ests), fill = TRUE)
  }

  plt <- ggplot(sim_storage, aes(x = sim_num, y = ATE, color = est)) + geom_line() +
    ggtitle(sprintf("Number of Strata = %d", nstrata))
  print(plt)
  
}
# ggplot(sim_storage, aes(x = sim_num, y = ATE, color = est)) + geom_line()

##### COPY OF OTHER FUNCTIONS ####
#' # APPENDIX: Additional Functions Used  
#' In this section, we include a number of other functions necessary to run this code. We do not include any functions defined in any of the tutorials. 
#+ echo=TRUE
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

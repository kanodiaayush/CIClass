#' ---
#' title: "Problem Set 1"
#' author: "Mitchell Linegar"
#' toc: true
#' toc_depth: 2
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
lambda <- c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1, 5, 10, 50, 100, 1000)

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
#'
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
#' This problem set examines the `r dataset_name` data set. Throughout, the treatment variable will be referred to as `r treatment_variable` and the outcome variable will be referred to as `r outcome_variable`. 

#### DATA WORK ####
# follows tutorial exactly
#' ## Pre-Processing the Data
#' I use code from a set of AtheyLab tutorials, which include the following note:  
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
#### TESTING ASSUMPTIONS ####
df_mod = copy(df)
setDT(df_mod)
Xmod = df_mod[,covariate_names, with=FALSE]
Ymod = df_mod$Y
Wmod = df_mod$W

# Computing the propensity score by logistic regression of W on X.
pW_logistic.fit <- glm(Wmod ~ as.matrix(Xmod), family = "binomial")
pW_logistic <- predict(pW_logistic.fit, type = "response")
pW_logistic.fit.tidy <-  pW_logistic.fit %>% tidy()
hist(pW_logistic)

logistic.fit <- glm(Ymod ~ as.matrix(Xmod), family = "binomial")
logistic <- predict(logistic.fit, type = "response")
logistic.fit.tidy <-  logistic.fit %>% tidy()
hist(logistic)

df_mod[, `:=`(p_Y = logistic,
              p_W = pW_logistic)]

#' 
{plot(smooth.spline(pW_logistic, Wmod, df = 4))
  abline(0, 1)}
#### RCT ANALYSIS ####
#' ## RCT Analysis
#' We first report the (presumably true) treatment effect $\hat{\tau}$ from the randomized experiment:  
tauhat_rct <- difference_in_means(df)
print(tauhat_rct)

#### SAMPLING BIAS ####
#' ## Introducing Bias
#+ echo=false
# FIXME: this doesn't seem to add enough bias
# these cutoffs courtesy of Kaleb (poc, .45 and .01)
prob = .5
# prob_control = prob # how much of the sample are we dropping?
# prop_treatment = prob
prob_control = 0.45  # how much of the sample are we dropping?
prob_treatment = 0.01

#' We now introduce sampling bias in order to simulate the situation we would be in if our data was from an observational study rather than a randomized experiment. 
#' This situation might arise due to sampling error or selection bias, and we will be able to see how various methods correct for this induced bias. 
#' To do, so, we under-sample treated units matching the following rule, and under-sample control units in its complement:  
#' FIXME: replace with POC discussion
#'       - Republican-leaning independent or closer to the Republican party (`partyid >= 4 & partyid < 7`)  
#'       - More than slightly conservative (`polviews >= 5.5 & polviews < 8`)  
#'       - with self-reported prestige of occuption above the bottom quartile (`prestg80 >= 34`)  
#' We remove `r prob * 100` percent of observations in these sets. 
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

tauhat_logistic_ipw 
tauhat_pscore_ols 
tauhat_lin_logistic_aipw 

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

#### OVERLAP ####
#' We now plot (logistic) propensity scores, showing that we still have overlap after removing observations: 
overlap <- df_mod %>% ggplot(aes(x=logistic_propensity,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()
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
pW_rf = predict(pW_rf.fit, newdata = Xmod)
pY_rf = predict(pY_rf.fit, newdata = Xmod)

# random forest, expanded data
pW_rf.fit.int = regression_forest(Xmod.int, Wmod, num.trees = 500)
pY_rf.fit.int = regression_forest(Xmod.int, Ymod, num.trees = 500)

# pW_rf.int = pW_rf.fit.int$predictions
# pY_rf.int = pY_rf.fit.int$predictions
pW_rf.int = predict(pW_rf.fit.int, newdata = Xmod)
pY_rf.int = predict(pY_rf.fit.int, newdata = Xmod)

# hist(pW_rf)
#### BIAS FUNCTION ####
#' Next we plot the bias function $b(X)$ following Athey, Imbens, Pham and Wager (AER P&P, 2017, Section IIID). 
#' We plot $b(x)$ for all units in the sample, and see that the bias seems evenly distributed around zero.  
#+ echo=false

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

#### ESTIMATING ATE INTRO ####
#' \newpage
#' # Estimating the ATE  
#' In this section we explore various methods for estimating the ATE. We explore the following methods:  
#' 1. inverse propensity weighting via logistic regression  
#' 2. direct regression analysis via OLS  
#' 3. traditional double robust analysis via augmented inverse-propensity score weighting that combines the above two estimators.  
#' We also re-run the above methods after expanding the data to include all interactions of all of the covariates, 
#' and re-estimate outcome and proensity models using the original linear model, as well as running lasso and random forest models on the expanded data. 

#### PLOTTING PROPENSITY PERFORMANCE ####
#' \newpage
#' We first plot propensity scores against treatment status on the original and expanded set of coefficients. 
#' Models closer to the 45-degree line are better.  
#+ echo = false

plot_prob(pW_logistic, Wmod, "Logistic")
plot_prob(pW_lasso.min, Wmod, "Lasso")
plot_prob(pW_rf, Wmod, "Random Forest")
plot_prob(pW_logistic.int, Wmod, "Logistic", "Expanded")
plot_prob(pW_lasso.int.min, Wmod, "Lasso", "Expanded")
plot_prob(pW_rf.int, Wmod, "Random Forest", "Expanded")

#### EXPLORING LASSO ####
### Exploring the Lasso Model Along Lambda
#' \newpage
#' To show how cross-validating lambda is important for the lasso, we compare predicted and actual treatment status for the minimum, maximum, a randomly selected, and the lambda that minimizes MSE (already shown above). 
#' The lasso with the best lamda is the one closest to the 45-degree line. 
#+ echo = false
plot_prob(pW_lasso.int.min, Wmod, "Lasso with Min Lambda", "Expanded")
plot_prob(pW_lasso.int.max, Wmod, "Lasso with Max Lambda", "Expanded")
plot_prob(pW_lasso.int.rand, Wmod, "Lasso with Random Lambda", "Expanded")
plot_prob(pW_lasso.int, Wmod, "Lasso with Best Lambda", "Expanded")

# {plot(smooth.spline(pW_logistic, Wmod, df = 4))
#   abline(0, 1)}

# plot lasso over grid of lambdas
pW_glmnet.fit.propensity.int.lambda_preds <- as.data.table(pW_glmnet.fit.propensity.int$fit.preval)
pW_glmnet.fit.propensity.int.lambda_preds <- pW_glmnet.fit.propensity.int.lambda_preds[
  # see discussion in FIXME above about using convert_to_prob here
  ,lapply(.SD, convert_to_prob), .SDcols = names(pW_glmnet.fit.propensity.int.lambda_preds)]

# credit to Kaleb for this formulation
<<<<<<< HEAD
 tauhat_lasso_ipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
  data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=ipw(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
}))[, model := "ipw"]
 
 tauhat_lasso_prop_score.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
   data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=prop_score_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
 }))[, model := "prop_score"]
 
 tauhat_lasso_aipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
   data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=aipw_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
 }))[, model := "aipw"]
=======
tauhat_lasso_ipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
  data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=ipw(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
}))[, model := "ipw"]

tauhat_lasso_prop_score.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
  data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=prop_score_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
}))[, model := "prop_score"]

tauhat_lasso_aipw.lambdas <- rbindlist(lapply(1:ncol(pW_glmnet.fit.propensity.int.lambda_preds), function(p){
  data.frame(lambda=pW_glmnet.fit.propensity.int$lambda[p],"ATE"=aipw_ols(df_mod.int, as.matrix(pW_glmnet.fit.propensity.int.lambda_preds[,..p]))["ATE"])
}))[, model := "aipw"]
>>>>>>> 876f2a3088166c76d799775073465eb50fa45441

tauhat_lass_estimates.lambdas <- rbindlist(list(tauhat_lasso_ipw.lambdas,
                                                tauhat_lasso_prop_score.lambdas,
                                                tauhat_lasso_aipw.lambdas))

#' We now plot our lasso estimates of $\hat{\tau}$ over our lambdas. 
#+ echo=FALSE
ggplot(tauhat_lass_estimates.lambdas, aes(x = lambda, y = ATE, color = model)) + 
  geom_line() + 
  geom_abline(aes(slope = 0, intercept = tauhat_rct["ATE"]))

# likelihoods over lambdas
lambda_log_liks <- pW_glmnet.fit.propensity.int.lambda_preds[
<<<<<<< HEAD
    ,lapply(.SD, loglike, Wmod), .SDcols = names(pW_glmnet.fit.propensity.int.lambda_preds)]
=======
  ,lapply(.SD, loglike, Wmod), .SDcols = names(pW_glmnet.fit.propensity.int.lambda_preds)]
>>>>>>> 876f2a3088166c76d799775073465eb50fa45441

colnames(lambda_log_liks) <- as.character(pW_glmnet.fit.propensity.int$lambda)
lambda_log_liks.long <- melt(lambda_log_liks, variable.name = "lambda", value.name = "llh")
lambda_log_liks.long[, lambda := as.numeric(as.character(lambda))]

#' We now plot log likelihood over different values of lambda on the expanded data.  
#+ echo=false
ggplot(lambda_log_liks.long, aes(x = lambda, y = llh)) + 
  geom_line() + 
  labs(title = "Log-Likelihood over Lambdas for Predicting W")

#### EXPLORING RF ####
# FIXME: copy analysis from above for RF, ideally varying N
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

# RF
tauhat_rf_ipw.int = ipw(df_mod.int, pW_rf.int)
ate_rf_aipw.int = average_treatment_effect(cf.int)
tauhat_rf_aipw.int = c(ATE=ate_rf_aipw.int["estimate"],
                       lower_ci=ate_rf_aipw.int["estimate"] - 1.96 * ate_rf_aipw.int["std.err"],
                       upper_ci=ate_rf_aipw.int["estimate"] + 1.96 * ate_rf_aipw.int["std.err"])
tauhat_ols_rf_aipw.int = aipw_ols(df_mod.int, pW_rf.int)

#### COMPARING ATE ACROSS MODELS ####
#' ## Comparing ATE Across Models with Original Data
#+ echo=false
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
#+ echo=false
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
#' ## Justification of Propensity Stratification
#' 


#### PROPENSITY STRATIFICATION FUNCTION ####
#' ## Propensity Stratification Function
#' Here we present a function to calculate propensity scores at a strata-level. 
#' The user supplies either a dataframe with treatment column $W$ and outcome column $Y$, 
#' as well as a model with which to estimate propensity scores 
#' (we don't supply a pre-calculated set of propensity scores in case we want to examine the usefulness of propensity stratification for new data). 
#' The function takes options for a number of strata and the function to estimate on the strata - 
#' the user could supply something more complex than the simple difference-in-means, for example.  
#' 
#' The function checks that each strata has both treatment groups; if it does not then it is not included in the ATE calculation.  
#+ echo=TRUE

propensity_stratification <- function(df, treatment_model, n_strata = 10, tau_estimator = difference_in_means){
  df <- copy(df)
  df[, pW := predict(treatment_model, newdata = df[,.SD, .SDcols = !c('W', 'Y')], type = "response")]
  df[, pW_strata := ntile(pW, n_strata)]
  # if any strata is empty, automatically ignored
  # if all
  strata_count <- df[, .N, by = .(pW_strata, W)]
  strata_to_keep <- strata_count[, .(n_strata_W_nonempty = sum(N > 0)), by = .(pW_strata)][
    n_strata_W_nonempty == 2, pW_strata] %>% unique # keep strata if has observations in both treatment
  
  # sloppy but can't figure out right solution atm
  strata_tau_est <- df[pW_strata %in% strata_to_keep, .(tauhat = tau_estimator(.SD)[1],
                                                        lower_ci = tau_estimator(.SD)[2],
                                                        upper_ci = tau_estimator(.SD)[3]), by = .(pW_strata)][
                                                          ,.(tauhat = mean(tauhat),
                                                             lower_ci = mean(lower_ci),
                                                             upper_ci = mean(upper_ci))]
  return(c(ATE = strata_tau_est$tauhat, lower_ci = strata_tau_est$lower_ci, upper_ci = strata_tau_est$upper_ci))
}

#### SIMULATION ####
#' ## Simulation Exercise
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

df <- make_simulation()
difference_in_means(df)

tauhat_propensity_stratification <- propensity_stratification(df_mod, pW_logistic.fit)
tauhat_propensity_stratification
#' compare this to 
# to get identical behavior (to be able to supply predictions to newdata for the propensity score model), run this line first:
# pW_logistic <- predict(pW_logistic.fit, newdata = df_mod, type = "response")
tauhat_logistic_ipw <- ipw(df_mod, pW_logistic)


sim_pW_logistic.fit <- glm(W ~ ., data = df %>% select(-Y), family = "binomial")
sim_pW_logistic <- predict(sim_pW_logistic.fit, type = "response")

difference_in_means(df)
ipw(df, sim_pW_logistic)
propensity_stratification(df, sim_pW_logistic.fit)

#' We now run `r n_sims` simulations of the type described above, and report results over each simulation.  
#+ echo=FALSE

sim_storage <- data.table()
for (i in 1:n_sims){
  df <- make_simulation()
  
  sim_pW_logistic.fit <- glm(W ~ ., data = df %>% select(-Y), family = "binomial")
  sim_pW_logistic <- predict(sim_pW_logistic.fit, type = "response")
  
  tau_ests <- rbind(
    c(df[,mean(W * V2)], NA, NA),
    difference_in_means(df),
    ipw(df, sim_pW_logistic),
    propensity_stratification(df, sim_pW_logistic.fit)  
  ) %>% as.data.frame()
  
  tau_ests$sim_num <- i
  tau_ests$est <- c("true", "difference_in_means", "ipw", "propensity_stratification")
  sim_storage <- rbindlist(list(sim_storage, tau_ests), fill = TRUE)
}

ggplot(sim_storage, aes(x = sim_num, y = ATE, color = est)) + geom_line()

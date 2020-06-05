library(CBPS)

#### MY FUNCS ####

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
    p_W_CBPS <- CBPS(Wmod ~ as.matrix(Xmod), method = "exact")$weights
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

run_all_models_on_df <- function(df, df_propensity = NULL, continuous_treatment = FALSE){
  # calculate pW_rf
  if(is.null(df_propensity)){
    df_propensity <- calculate_propensities(df, continuous_treatment = continuous_treatment)
  }
  
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
  # print(round(as.matrix(all_estimators), 3))
}

#### LOAD DATA ####


freadom <- fread(sprintf("%s/treatment_effects.csv", data_dir))

# cleaning
# some W not really continuous treatment - just different levels
freadom[,W_reading_time] %>% hist()
freadom$W_wordcount %>% hist()

# for now only look at 'medium' and 'high' levels of W_reading_time
freadom <- freadom[W_reading_time==12.5 | W_reading_time==7.5]
# create binary variable
freadom[, W_reading_time_high := ifelse(W_reading_time==max(W_reading_time), 1, 0)]

covariate_names <- colnames(freadom)[!grepl("W|Y", colnames(freadom), ignore.case = FALSE)]

# construct various data sets we'll need
df_time_finish <- freadom[,c(covariate_names, "W_reading_time_high", "Y_utility"), with=FALSE]
df_time_finish %>% setnames_W_Y()

df_wcount_finish <- freadom[,c(covariate_names, "W_wordcount", "Y_utility"), with=FALSE]
df_wcount_finish %>% setnames_W_Y()

df_time_lag <- freadom[,c(covariate_names, "W_reading_time_high", "Y_next_view_lag"), with=FALSE]
df_time_lag %>% setnames_W_Y()

df_wcount_lag <- freadom[,c(covariate_names, "W_wordcount", "Y_next_view_lag"), with=FALSE]
df_wcount_lag %>% setnames_W_Y()

df_finish_lag <- freadom[,c(covariate_names, "W_utility", "Y_next_view_lag"), with=FALSE]
df_finish_lag %>% setnames_W_Y()


#### ANALYSIS PER DATASET ####

# Calculate and Plot Overlap
df_time_finish_overlap <- calculate_propensities(df_time_finish)
df_time_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()

df_wcount_finish_overlap <- calculate_propensities(df_wcount_finish, continuous_treatment = TRUE)
df_wcount_finish_overlap %>% ggplot(aes(x=p_W_rf,color=as.factor(W),fill=as.factor(W)))+ geom_histogram()

# Subset if necessary

# Run models




#### READING TIME FINISH SOD ANALYSIS ####
#' # Rate at Which Users Finish Stories of the Day by Length {#finish_sod_by_length}
#+ finish_sod_ate

df_time_finish_models <- run_all_models_on_df(df_time_finish, df_time_finish_overlap)
print(round(as.matrix(df_time_finish_models), 3))


df_wcount_finish_overlap <- calculate_propensities(df_wcount_finish, continuous_treatment = TRUE)
df_wcount_finish_models <- run_all_models_on_df(df_wcount_finish, df_wcount_finish_overlap, continuous_treatment = TRUE)
print(round(as.matrix(df_wcount_finish_models), 3))


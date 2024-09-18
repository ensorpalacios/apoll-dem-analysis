#' Cox regression
#'
#' Run survival analysis on dementia data;
#' address_buff specifies minimum years participants must have lived at baseline address;
#' dem_buff specifies minimum years before incidence of dementia is accepted;
#' outcome: either date_dem (all cause dementia), date_ad (Alzheimer) or date_vd (vascular 
#' dementia)
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#f
#' ----------------------------------------------------------------------------

# install and import libraries
install.packages('here')
install.packages('finalfit')
install.packages('lme4')
install.packages('kableExtra')
install.packages('coxme')
install.packages('rms')
install.packages('glmnet')

library(data.table)
library(tidyverse)
library(here)
library(survival)
library(coxme)
library(finalfit)
library(kableExtra)
library(glmnet)
library(MASS)
library(devtools)
library(remotes)
install_version('rms', version = '6.8-1', repo = 'https://lib.stat.cmu.edu/R/CRAN', upgrade = F)
library(rms)

source(paste0(here(), '/src/select.R'))

# Load and select data --------------------------------------------------------
# Selection parameters
address_buff <- 1 # discard participants with < x years at baseline address
dem_buff <- 0 # discard participants with dem incidence < x years after baseline
list_covariates <- c('age', # for main model (2)
                     'sex', 
                     'ethnic', 
                     'ea',
                     'income', # causes ~ 70k to drop!
                     'pop_density',
                     'mdi_eng')
# Get data
data_path <- paste0(here(), '/data/processed/df_ukbb.rds')
df_all <- readRDS(file = data_path)
dummy <- select_data(df_all, 
                     address_buff = address_buff,
                     dem_buff = dem_buff,
                     list_covariates)
df_sel <- dummy[[1]]

# Save numbers of complete case analysis
cca <- dummy[c('n_all',
               'apoll_keep', 'apoll_drop',
               'dem_keep', 'dem_drop',
               'cov_keep_all', 'cov_drop_all', 
               'cov_drop', 'x100mdi_eng',
               'address_keep', 'address_drop' # address_keep: total n for analysis
               )]

# Compute time and status variables for survival analysis ---------------------
# Because end_study (follow-up) comes earlier that last records, don't consider
# dem incidence after end_study when computing status
for (outcome in c('dem', 'ad', 'vd')) {
    ap_var = paste0('date_', outcome)
    df_sel <- df_sel %>% mutate('time_{outcome}' := pmin(get(ap_var),
                                                         date_lost_fu,
                                                         date_death,
                                                         end_study,
                                                         na.rm = T) - date_baseline,
                                'time_{outcome}' := time_length(get(str_glue('time_{outcome}')), 'years'),
                                # status 1 if dem occurs, except if it date comes after end of study
                                'status_{outcome}' := case_when(get(ap_var) <= end_study ~ 1L,
                                                                get(ap_var) > end_study ~ 0L,
                                                                is.na(get(ap_var)) ~ 0L)
    )
}

# Prepare lists and model idx -------------------------------------------------
# Setup lists and parameters
list_fit <- as.list(rep(NA, 500)) # initialise lists with enough elements
list_aic <- as.list(rep(NA, 500))
list_test <- as.list(rep(NA, 500))
list_res <- as.list(rep(NA, 500))
list_spline <- as.list(rep(NA, 500))
list_vif <- as.list(rep(NA, 500))

idx <- 1 # global index for list_fit and list_test; run models in order!
ap_list <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'no')
ap_all <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'no', 'all', 'appc')

# Run basic cox regression ----------------------------------------------------
# (1) Main model (include appc_all - only for model 1)
for (outcome in c('dem', 'ad', 'vd')) {
    for (ap_type in c(ap_all, 'appc_all')) {
        for (ap_unit in c('s', 'q')) {
            # Define exposure
            if (ap_type== 'all') { # For model with all ap included
                apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                               collapse = ' + ')
            } else if (ap_type %in% ap_list) { # for single pollutant models
                apoll = str_glue('{ap_type}_{ap_unit}')
            } else { # for appc and appc_all
                apoll <- ap_type
            }

            if (!(grepl('appc', ap_type) & ap_unit == 'q')) { # do not run for appc-q (not existing)
                # Define Model
                model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~ 
                                  age + 
                                  sex + 
                                  ethnic + 
                                  ea + 
                                  income + 
                                  pop_density +
                                  {apoll}')
                model <- as.formula(model)

                # Fit and save model
                fit <- coxph(model, data = df_sel, method = 'efron')
                name <- str_glue('model1_{outcome}_{ap_unit}_{ap_type}')
                list_fit[[idx]] <- fit
                names(list_fit)[idx] <- name

                # AIC
                list_aic[[idx]] <- AIC(fit)
                names(list_aic)[idx] <- name

                # Proportionality tests
                list_test[[idx]] <- fit %>% cox.zph
                names(list_test)[idx] <- name 

                # Deviance residuals
                list_res[[idx]] <- fit %>% residuals(type = 'deviance')
                names(list_res)[idx] <- name

                # Compute variation inflation factor for 'all' model
                if (ap_type == 'all') {
                    list_vif[[idx]] <- vif(fit)
                    names(list_vif)[idx] <- name
                }

                # Update global index
                idx = idx + 1
            }
        }
    }
}

# (2) Same as (1) but add multiple deprivation index quartiles for England (mdi_eng_q)
for (outcome in c('dem', 'ad', 'vd')) {
    for (ap_type in ap_all) {
        for (ap_unit in c('s', 'q')) {
            # Define exposure
            if (ap_type== 'all') { # For model with all ap included
                apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                               collapse = ' + ')
            } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
                apoll = str_glue('{ap_type}_{ap_unit}')
            } else { # for appc
                apoll <- ap_type
            }

            if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
                # Define Model
                model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~ 
                                  age + 
                                  sex + 
                                  ethnic + 
                                  ea + 
                                  income + 
                                  pop_density +
                                  mdi_eng_q +
                                  {apoll}')
                model <- as.formula(model)

                # Fit and save model
                fit <- coxph(model, data = df_sel, method = 'efron')
                name <- str_glue('model2_{outcome}_{ap_unit}_{ap_type}')
                list_fit[[idx]] <- fit
                names(list_fit)[idx] <- name

                # AIC
                list_aic[[idx]] <- AIC(fit)
                names(list_aic)[idx] <- name

                # Deviance residuals
                list_res[[idx]] <- fit %>% residuals(type = 'deviance')
                names(list_res)[idx] <- name

                # Compute variation inflation factor for 'all' model
                if (ap_type == 'all') {
                    list_vif[[idx]] <- vif(fit)
                    names(list_vif)[idx] <- name
                }

                # Update global index
                idx = idx + 1
            }
        }
    }
}

# (3) Same as (2) but drop people who lived less than 5 years at baseline address
# - get data with address_buff = 5
address_buff <- 5 # discard participants with < x years at baseline address
dummy5 <- select_data(df_all, 
                     address_buff = address_buff,
                     dem_buff = dem_buff,
                     list_covariates)
df_sel5 <- dummy5[[1]]

# - add address keep/drop with 5 years to cca
cca[['address_keep5']] <- dummy5[['address_keep']] # total n for 5 years buffer analysis
cca[['address_drop5']] <- dummy5[['address_drop']]

# - Compute time and status variables for survival analysis
for (outcome in c('dem', 'ad', 'vd')) {
    ap_var = paste0('date_', outcome)
    df_sel5 <- df_sel5 %>% mutate('time_{outcome}' := pmin(get(ap_var),
                                                         date_lost_fu,
                                                         date_death,
                                                         end_study,
                                                         na.rm = T) - date_baseline,
                                'time_{outcome}' := time_length(get(str_glue('time_{outcome}')), 'years'),
                                # status 1 if dem occurs, except if it date comes after end of study
                                'status_{outcome}' := case_when(get(ap_var) <= end_study ~ 1L,
                                                                get(ap_var) > end_study ~ 0L,
                                                                is.na(get(ap_var)) ~ 0L)
    )
}


# - run analysis
for (outcome in c('dem', 'ad', 'vd')) {
    for (ap_type in ap_all) {
        for (ap_unit in c('s', 'q')) {
            # Define exposure
            if (ap_type== 'all') { # For model with all ap included
                apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                               collapse = ' + ')
            } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
                apoll = str_glue('{ap_type}_{ap_unit}')
            } else { # for appc
                apoll <- ap_type
            }

            if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
                # Define Model
                model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~ 
                                  age + 
                                  sex + 
                                  ethnic + 
                                  ea + 
                                  income + 
                                  pop_density +
                                  mdi_eng_q +
                                  {apoll}')
                model <- as.formula(model)

                # Fit and save model
                fit <- coxph(model, data = df_sel5, method = 'efron')
                name <- str_glue('model3_{outcome}_{ap_unit}_{ap_type}')
                list_fit[[idx]] <- fit
                names(list_fit)[idx] <- name

                # AIC
                list_aic[[idx]] <- AIC(fit)
                names(list_aic)[idx] <- name

                # Compute variation inflation factor for 'all' model
                if (ap_type == 'all') {
                    list_vif[[idx]] <- vif(fit)
                    names(list_vif)[idx] <- name
                }

                # Update global index
                idx = idx + 1
            }
        }
    }
}


# Random effect model for recruitment centre ----------------------------------
# (4) Same as (2), include mdi_eng; exclude model with all 
# ap (-7); run only for apoll coded as continuous
for (ap_type in ap_all[-7]) {
    if (ap_type %in% ap_list) { # for single pollutant models
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Model
    model <- str_glue('Surv(time_dem, status_dem) ~
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      ({apoll} + 1 | rec_centre)')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxme(model, data = df_sel)
    name <- str_glue('model4_dem_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Update global index
    idx = idx + 1
}

# Save n participants for each rec_centre
cca[['n_rec_centre']] <- df_sel[['rec_centre']] %>% table

# (5) Same as (4) but exclude participants who lived less than 5 years at baseline address
for (ap_type in ap_all[-7]) {
    if (ap_type %in% ap_list) { # for single pollutant models
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Model
    model <- str_glue('Surv(time_dem, status_dem) ~
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      ({apoll} + 1 | rec_centre)')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxme(model, data = df_sel5)
    name <- str_glue('model5_dem_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Update global index
    idx = idx + 1
}

# Save n participants for each rec_centre (address_buff = 5)
cca[['n_rec_centre5']] <- df_sel5[['rec_centre']] %>% table

## Run Mundlack regression for spatial confounding -----------------------------
## (6) Mundlack regression with recruitment centres; include mdi_eng
## - ap mean within recruitment centre
#for (apoll in c('pm25_s', 'pmcoarse_s', 'pmabs_s', 'pm10_s', 'no2_s', 'no_s', 'appc')) {
#    df_sel <- df_sel %>%
#        group_by(rec_centre) %>% 
#        mutate('{apoll}_mean' := mean(get(apoll), na.rm = T)) %>%
#        ungroup() # necessary for use of glmnet as groupby adds attribute 'groups'
#}
#
## - run analysis; exclude model with all ap (-7); run only for apoll coded as continuous
#for (outcome in c('dem')) {
#    for (ap_type in ap_all[-7]) {
#        if (ap_type %in% ap_list) { # for single pollutant models
#            apoll = str_glue('{ap_type}_s')
#        } else { # for appc
#            apoll <- ap_type
#        }
#
#        # Model
#        model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~
#                          age + 
#                          sex + 
#                          ethnic + 
#                          ea + 
#                          income + 
#                          pop_density +
#                          mdi_eng_q +
#                          {apoll} +
#                          {apoll}_mean +
#                          ({apoll} + 1 | rec_centre)')
#        model <- as.formula(model)
#
#        # Fit and save model
#        fit <- coxme(model, data = df_sel)
#        name <- str_glue('model6_{outcome}_s_{ap_type}')
#        list_fit[[idx]] <- fit
#        names(list_fit)[idx] <- name
#
#        # AIC
#        list_aic[[idx]] <- AIC(fit)
#        names(list_aic)[idx] <- name
#
#        # Update global index
#        idx = idx + 1
#    }
#}
#
## (7) Same as (6) but exclude participants who lived less than 5 years at baseline address
## - ap mean within recruitment centre
#for (apoll in c('pm25_s', 'pmcoarse_s', 'pmabs_s', 'pm10_s', 'no2_s', 'no_s', 'appc')) {
#    df_sel5 <- df_sel5 %>%
#        group_by(rec_centre) %>% 
#        mutate('{apoll}_mean' := mean(get(apoll), na.rm = T)) %>%
#        ungroup() # necessary for use of glmnet as groupby adds attribute 'groups'
#}
#
## - run analysis
#for (outcome in c('dem')) {
#    for (ap_type in ap_all[-7]) {
#        if (ap_type %in% ap_list) { # for single pollutant models
#            apoll = str_glue('{ap_type}_s')
#        } else { # for appc
#            apoll <- ap_type
#        }
#
#        # Model
#        model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~
#                          age + 
#                          sex + 
#                          ethnic + 
#                          ea + 
#                          income + 
#                          pop_density +
#                          mdi_eng_q +
#                          {apoll} +
#                          {apoll}_mean +
#                          ({apoll} + 1 | rec_centre)')
#        model <- as.formula(model)
#
#        # Fit and save model
#        fit <- coxme(model, data = df_sel5)
#        name <- str_glue('model7_{outcome}_s_{ap_type}')
#        list_fit[[idx]] <- fit
#        names(list_fit)[idx] <- name
#
#        # AIC
#        list_aic[[idx]] <- AIC(fit)
#        names(list_aic)[idx] <- name
#
#        # Update global index
#        idx = idx + 1
#    }
#}

# Interaction analysis --------------------------------------------------------
# (8) Model (1) + age/sex/income time interactions
# restrict to date_dem (all cause dementia) and appc
# Split data for different time groups; choose cut based on Schoenfeld residual plots
df_split <- survSplit(Surv(time_dem, status_dem) ~ ., data = df_sel, cut = c(12), episode = 'tgroup')
# list_tint <- c('age', 'sex', 'appc') # variables interacting with time
list_tint <- c('age', 'sex', 'income') # variables interacting with time
for (var_tint in list_tint) {
    # Model
    model <- str_glue('Surv(time_dem, status_dem) ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      appc + 
                      {var_tint}:strata(tgroup)')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxph(model, data = df_split, method='efron')
    name <- str_glue('model8_dem_s_appc_{var_tint}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Update global index
    idx = idx + 1
}


# Positive control analysis (copd) --------------------------------------------
# (9) Same as (1) but with copd (no mdi_eng)
# - get data for copd (address buffer 1)
address_buff <- 1
dummy_copd <- select_data(df_all, 
                         address_buff = address_buff,
                         dem_buff = dem_buff,
                         copd = T,
                         list_covariates)
df_sel_copd <- dummy_copd[[1]]

# - Compute time and status variables for survival analysis
df_sel_copd <- df_sel_copd %>% mutate('time_copd' = pmin(date_copd,
                                                         date_lost_fu,
                                                         date_death,
                                                         end_study,
                                                         na.rm = T) - date_baseline,
                                      'time_copd' = time_length(time_copd, 'years'),
                                      # status 1 if copd occurs, except if it date comes after end of study
                                      'status_copd' = case_when(date_copd <= end_study ~ 1L,
                                                                date_copd > end_study ~ 0L,
                                                                is.na(date_copd) ~ 0L)
)

# - add copd numbers to cca
cca[['copd_keep']] <- dummy_copd[['copd_keep']]
cca[['copd_drop']] <- dummy_copd[['copd_drop']]
cca[['copd_cases']] <- dummy_copd[['n_copd']]
cca[['copd_address_keep']] <- dummy_copd[['address_keep']] # total n for copd analysis

# - run analysis
for (ap_type in ap_all) {
    for (ap_unit in c('s', 'q')) {
        # Define exposure
        if (ap_type== 'all') { # For model with all ap included
            apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                           collapse = ' + ')
        } else if (ap_type %in% ap_list) { # for single pollutant models
            apoll = str_glue('{ap_type}_{ap_unit}')
        } else { # for appc
            apoll <- ap_type
        }

        if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
            # Define Model
            model <- str_glue('Surv(time_copd, status_copd) ~ 
                              age + 
                              sex + 
                              ethnic + 
                              ea + 
                              income + 
                              pop_density +
                              {apoll}')
            model <- as.formula(model)

            # Fit and save model
            fit <- coxph(model, data = df_sel_copd, method = 'efron')
            name <- str_glue('model9_copd_{ap_unit}_{ap_type}')
            list_fit[[idx]] <- fit
            names(list_fit)[idx] <- name

            # Update global index
            idx = idx + 1
        }
    }
}


# (10) Same as (9) but with mdi_eng
for (ap_type in ap_all) {
    for (ap_unit in c('s', 'q')) {
        # Define exposure
        if (ap_type== 'all') { # For model with all ap included
            apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                           collapse = ' + ')
        } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
            apoll = str_glue('{ap_type}_{ap_unit}')
        } else { # for appc
            apoll <- ap_type
        }

        if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
            # Define Model
            model <- str_glue('Surv(time_copd, status_copd) ~ 
                              age + 
                              sex + 
                              ethnic + 
                              ea + 
                              income + 
                              pop_density +
                              mdi_eng_q +
                              {apoll}')
            model <- as.formula(model)

            # Fit and save model
            fit <- coxph(model, data = df_sel_copd, method = 'efron')
            name <- str_glue('model10_copd_{ap_unit}_{ap_type}')
            list_fit[[idx]] <- fit
            names(list_fit)[idx] <- name

            # Update global index
            idx = idx + 1
        }
    }
}

# (11) Same as (10) but exclude participants who lived less than 5 years at baseline address
# - get data for copd (address buffer 5)
address_buff <- 5
dummy_copd5 <- select_data(df_all, 
                         address_buff = address_buff,
                         dem_buff = dem_buff,
                         copd = T,
                         list_covariates)
df_sel5_copd <- dummy_copd5[[1]]

# - Compute time and status variables for survival analysis
df_sel5_copd <- df_sel5_copd %>% mutate('time_copd' = pmin(date_copd,
                                                           date_lost_fu,
                                                           date_death,
                                                           end_study,
                                                           na.rm = T) - date_baseline,
                                        'time_copd' = time_length(time_copd, 'years'),
                                        # status 1 if copd occurs, except if it date comes after end of study
                                        'status_copd' = case_when(date_copd <= end_study ~ 1L,
                                                                  date_copd > end_study ~ 0L,
                                                                  is.na(date_copd) ~ 0L)
)

# - add address keep/drop with 5 years to cca
cca[['copd_address_keep5']] <- dummy5[['address_keep']] # total n for 5 years buffer copd analysis
cca[['copd_address_drop5']] <- dummy5[['address_drop']]

# - run analysis
for (ap_type in ap_all) {
    for (ap_unit in c('s', 'q')) {
        # Define exposure
        if (ap_type== 'all') { # For model with all ap included
            apoll <- paste(paste0(ap_list, str_glue('_{ap_unit}')),
                           collapse = ' + ')
        } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
            apoll = str_glue('{ap_type}_{ap_unit}')
        } else { # for appc
            apoll <- ap_type
        }

        if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
            # Define Model
            model <- str_glue('Surv(time_copd, status_copd) ~ 
                              age + 
                              sex + 
                              ethnic + 
                              ea + 
                              income + 
                              pop_density +
                              mdi_eng_q +
                              {apoll}')
            model <- as.formula(model)

            # Fit and save model
            fit <- coxph(model, data = df_sel5_copd, method = 'efron')
            name <- str_glue('model11_copd_{ap_unit}_{ap_type}')
            list_fit[[idx]] <- fit
            names(list_fit)[idx] <- name

            # Update global index
            idx = idx + 1
        }
    }
}

# Negative control analysis (oily fish intake) -----------------------------------------
# (12) Same as (1) but using ordinal logistic regression on oily fish intake frequency (no mdi_eng)
# - get data for oily fish consumption
address_buff <- c(1, 200) # 200 is upper bound time lived at baseline address (to select all)
dummy_fish <- select_data(df_all, 
                            address_buff = address_buff,
                            dem_buff = dem_buff,
                            fish = T,
                            list_covariates)
df_sel_fish <- dummy_fish[[1]]

# - add oily fish intake numbers to cca
cca[['fish_keep']] <- dummy_fish[['fish_keep']]
cca[['fish_drop']] <- dummy_fish[['fish_drop']]
cca[['fish_address_keep']] <- dummy_fish[['address_keep']] # total n for copd analysis

# - run analysis
for (ap_type in ap_all) {
    # Define exposure
    if (ap_type== 'all') { # For model with all ap included
        apoll <- paste(paste0(ap_list, str_glue('_s')),
                       collapse = ' + ')
    } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Define Model
    model <- str_glue('fish ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      {apoll}')
    model <- as.formula(model)

    # Fit and save model
    fit <- polr(model, data = df_sel_fish, Hess = T)
    name <- str_glue('model12_fish_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # Update global index
    idx = idx + 1
}

# (13) Same as (12) but with mdi_eng_q
# - run analysis
for (ap_type in ap_all) {
    # Define exposure
    if (ap_type== 'all') { # For model with all ap included
        apoll <- paste(paste0(ap_list, str_glue('_s')),
                       collapse = ' + ')
    } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Define Model
    model <- str_glue('fish ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      {apoll}')
    model <- as.formula(model)

    # Fit and save model
    fit <- polr(model, data = df_sel_fish, Hess = T)
    name <- str_glue('model13_fish_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # Update global index
    idx = idx + 1
}


# (14) Same as (12) but restrict to participants with 1 <= time_cuaddress <= 2
# - get data for oily fish consumption
address_buff <- c(1, 2)
dummy_fish12 <- select_data(df_all, 
                            address_buff = address_buff,
                            dem_buff = dem_buff,
                            fish = T,
                            list_covariates)
df_sel_fish12 <- dummy_fish12[[1]]

# - add oily fish intake numbers to cca
cca[['fish_address_keep12']] <- dummy_fish12[['address_keep']] # total n for copd analysis
# - run analysis
for (ap_type in ap_all) {
    # Define exposure
    if (ap_type== 'all') { # For model with all ap included
        apoll <- paste(paste0(ap_list, str_glue('_s')),
                       collapse = ' + ')
    } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Define Model
    model <- str_glue('fish ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      {apoll}')
    model <- as.formula(model)

    # Fit and save model
    fit <- polr(model, data = df_sel_fish12, Hess = T)
    name <- str_glue('model14_fish_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # Update global index
    idx = idx + 1
}


# (15) Same as (13) but restrict to participants with 1 <= time_cuaddress <= 2
# - run analysis
for (ap_type in ap_all) {
    # Define exposure
    if (ap_type== 'all') { # For model with all ap included
        apoll <- paste(paste0(ap_list, str_glue('_s')),
                       collapse = ' + ')
    } else if (ap_type %in% ap_list) { # for single pollutant models (not existing)
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Define Model
    model <- str_glue('fish ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      {apoll}')
    model <- as.formula(model)

    # Fit and save model
    fit <- polr(model, data = df_sel_fish12, Hess = T)
    name <- str_glue('model15_fish_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # Update global index
    idx = idx + 1
}

# Lasso regression for variable selection -------------------------------------
# (16) As (2) + lasso penalisation (only for model with all pollutants
# included); address_buff = 1
x <- df_sel %>% 
    dplyr::select(age, sex, ethnic, ea, income, pop_density, mdi_eng_q, pm25, pmcoarse, pmabs, pm10, no2, no) %>%
    data.matrix()
lasso <- map(c('dem', 'ad', 'vd'), \(z) {
                  y <- with(df_sel, 
                            Surv(get(str_glue('time_{z}')), get(str_glue('status_{z}')), type = 'right'))
                  cv.glmnet(x, 
                            y, 
                            family = 'cox', 
                            penalty.factor = c(rep(0, 7), rep(1, 6)), 
                            type.measure = 'C')}) %>% setNames(c('dem', 'ad', 'vd'))
# - Add lasso to list_fit
walk(names(lasso), \(x) {
         list_fit[[idx]] <- lasso[[x]]
         names(list_fit)[idx] <- str_glue('model16_lasso_{x}')
         idx <- idx + 1
         #  Assign to global environment
         assign('list_fit', list_fit, envir = .GlobalEnv)
         assign('idx', idx, envir = .GlobalEnv)
                  })

# (17) As (3) + lasso penalisation (only for model with all pollutants 
#' included); address_buff = 5
x <- df_sel5 %>% 
    dplyr::select(age, sex, ethnic, ea, income, pop_density, mdi_eng_q, pm25, pmcoarse, pmabs, pm10, no2, no) %>%
    data.matrix()
lasso5 <- map(c('dem', 'ad', 'vd'), \(z) {
                  y <- with(df_sel5, 
                            Surv(get(str_glue('time_{z}')), get(str_glue('status_{z}')), type = 'right'))
                  cv.glmnet(x, 
                            y, 
                            family = 'cox', 
                            penalty.factor = c(rep(0, 7), rep(1, 6)), 
                            type.measure = 'C')}) %>% setNames(c('dem', 'ad', 'vd'))
# - Add lasso to list_fit
walk(names(lasso5), \(x) {
         list_fit[[idx]] <- lasso5[[x]]
         names(list_fit)[idx] <- str_glue('model17_lasso_{x}')
         idx <- idx + 1
         #  Assign to global environment
         assign('list_fit', list_fit, envir = .GlobalEnv)
         assign('idx', idx, envir = .GlobalEnv)
                  })

#aa <- cv.glmnet(x, y, family = 'cox', penalty.factor = c(rep(0, 6), rep(1, 6)), type.measure = 'C')
#aa <- cv.glmnet(x, y, family = 'cox', type.measure = 'C')
#aa %>% str
#exp(coef(aa, s=0.0002209))
#exp(coef(aa, s=0.0015582))
#aa[['glmnet.fit']][['beta']]
#aa$beta


# (18) Same as (3) plus noise pollution; exclude model with all ap (-7)
# - run analysis
for (ap_type in ap_all[-7]) {
    for (ap_unit in c('s', 'q')) {
        # Define exposure
        if (ap_type %in% ap_list) { # for single pollutant models
            apoll = str_glue('{ap_type}_{ap_unit}')
        } else { # for appc
            apoll <- ap_type
        }

        if (!(ap_type == 'appc' & ap_unit == 'q')) { # do not run for appc-q (not existing)
            # Define Model
            model <- str_glue('Surv(time_dem, status_dem) ~ 
                              age + 
                              sex + 
                              ethnic + 
                              ea + 
                              income + 
                              pop_density +
                              mdi_eng_q +
                              noise_poll_s +
                              {apoll}')
            model <- as.formula(model)

            # Fit and save model
            fit <- coxph(model, data = df_sel5, method = 'efron')
            name <- str_glue('model18_dem_{ap_unit}_{ap_type}')
            list_fit[[idx]] <- fit
            names(list_fit)[idx] <- name

            # AIC
            list_aic[[idx]] <- AIC(fit)
            names(list_aic)[idx] <- name

            # Compute variation inflation factor for 'all' model
            if (ap_type == 'all') {
                list_vif[[idx]] <- vif(fit)
                names(list_vif)[idx] <- name
            }

            # Update global index
            idx = idx + 1
        }
    }
}

# Interaction analysis --------------------------------------------------------
# (19) Same as (2) plus apoll-mdi interaction; exclude model with all ap (-7)
# - run analysis
for (ap_type in ap_all[-7]) {
    # Define exposure
    if (ap_type %in% ap_list) { # for single pollutant models
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

    # Define Model
    model <- str_glue('Surv(time_dem, status_dem) ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      {apoll} +
                      {apoll}:mdi_eng_q')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxph(model, data = df_sel, method = 'efron')
    name <- str_glue('model19_dem_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Deviance residuals
    list_res[[idx]] <- fit %>% residuals(type = 'deviance')
    names(list_res)[idx] <- name

    # Update global index
    idx = idx + 1
}

# (20) Same as (3) but plus apoll-mdi interaction; exclude model with all ap (-7)
# - run analysis
for (ap_type in ap_all[-7]) {
    # Define exposure
    if (ap_type %in% ap_list) { # for single pollutant models
        apoll = str_glue('{ap_type}_s')
    } else { # for appc
        apoll <- ap_type
    }

        # Define Model
    model <- str_glue('Surv(time_dem, status_dem) ~ 
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      mdi_eng_q +
                      {apoll} +
                      {apoll}:mdi_eng_q')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxph(model, data = df_sel5, method = 'efron')
    name <- str_glue('model20_dem_s_{ap_type}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Update global index
    idx = idx + 1
}

## Spline model ---------------------------------------------------------------
## (16) Model (2) + apoll non-linear effect
## do not run for ap quartiles
#for (outcome in c('dem', 'ad', 'vd')) {
#    # Model
#    model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~
#                      age + 
#                      sex + 
#                      ethnic + 
#                      ea + 
#                      income + 
#                      pop_density +
#                      pspline(appc, df = 3)')
#    model <- as.formula(model)
#
#    # Fit and save model
#    fit <- coxph(model, data = df_sel, method='efron')
#    name <- str_glue('model6_{outcome}_s_appc')
#    list_fit[[idx]] <- fit
#    names(list_fit)[idx] <- name
#
#    # AIC
#    list_aic[[idx]] <- AIC(fit)
#    names(list_aic)[idx] <- name
#
#    # Plot spline fit
#    list_spline[[idx]] <- termplot(fit, se = T, plot = F)
#    names(list_spline)[idx] <- name
#
#    # Update global index
#    idx = idx + 1
#}

# Save results ----------------------------------------------------------------
save_path <- here(str_glue('output/results/cox/savedata/'))
if (!file.exists(save_path)) {
    dir.create(save_path, recursive = T)
}

saveRDS(list_fit, file = paste0(save_path, 'fits.rds'))
saveRDS(list_test, file = paste0(save_path, 'tests.rds'))
saveRDS(list_res, file = paste0(save_path, 'res.rds'))
saveRDS(list_aic, file = paste0(save_path, 'aic.rds'))
saveRDS(list_spline, file = paste0(save_path, 'spline.rds'))
saveRDS(list_vif, file = paste0(save_path, 'vif.rds'))
saveRDS(cca, file = paste0(save_path, 'cca.rds'))



#' Cox regression
#'
#' Run survival analysis on dementia data;
#' address_buff specifies minimum years participants must have lived at baseline address;
#' dem_buff specifies minimum years before incidence of dementia is accepted;
#' outcome: either date_dem (all cause dementia), date_ad (Alzheimer) or date_vd (vascular 
#' dementia)
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#'
#' ----------------------------------------------------------------------------

# install and import libraries
install.packages('here')
install.packages('finalfit')
install.packages('lme4')
install.packages('kableExtra')
install.packages('coxme')
#install.packages('rstanarm')
#install.packages('rstan', repos = c('https://mc-stan.org/r-packages/', getoption('repos')))
#install.packages('brms')
#install.packages('splines2')

library(data.table)
library(tidyverse)
library(here)
library(survival)
library(coxme)
library(finalfit)
library(kableExtra)
#install.packages('lme4') # warhin about Matrix ABI version when loading finalfit - ignored
#library(brms)
#library(rstanarm)

source(paste0(here(), '/src/select.R'))

# Load and select data --------------------------------------------------------
# Selection parameters
address_buff <- 1 # discard participants with < x years at baseline address
dem_buff <- 0 # discard participants with dem incidence < x years after baseline
drop_nowhite <- F # restrict analysis to 'white'
list_covariates <- c('age', # for main model (2)
                     'sex', 
                     'ethnic', 
                     'ea',
                     'income', # causes ~ 70k to drop!
                     'pop_density')
# Get data
data_path <- paste0(here(), '/data/processed/df_ukbb.rds')
df_all <- readRDS(file = data_path)
dummy <- select_data(df_all, 
                     address_buff = address_buff,
                     dem_buff = dem_buff,
                     drop_nowhite = drop_nowhite,
                     list_covariates)
df_sel <- dummy[[1]]

# Compute time and status variables for survival analysis ---------------------
# Because end_study (follow-up) comes earlier that last records, don't consider
# dem incidence after end_study when computing status
tscale <- 'years'
for (outcome in c('dem', 'ad', 'vd')) {
    ap_var = paste0('date_', outcome)
    df_sel <- df_sel %>% mutate('time_{outcome}' := pmin(get(ap_var),
                                                         date_lost_fu,
                                                         date_death,
                                                         end_study,
                                                         na.rm = T) - date_baseline,
                                'time_{outcome}' := time_length(get(str_glue('time_{outcome}')), tscale),
                                'status_{outcome}' := case_when(get(ap_var) <= end_study ~ 1L,
                                                                get(ap_var) > end_study ~ 0L,
                                                                is.na(get(ap_var)) ~ 0L)
    )
}

# Run analysis and save in lists ----------------------------------------------
# Analysis involve:
# - Fit Cox model
# - Compute AIC
# - Run proportionality test
# - Plot spline fit

# Setup lists and parameters
list_fit <- as.list(rep(NA, 500)) # initialise lists with enough elements
list_aic <- as.list(rep(NA, 500))
list_test <- as.list(rep(NA, 500))
list_spline <- as.list(rep(NA, 500))
idx <- 1 # global index for list_fit and list_test; run models in order!
ap_list <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'nox')
ap_all <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'nox', 'all', 'appc')

# (1) Main model
for (outcome in c('dem', 'ad', 'vd')) {
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

                # Update global index
                idx = idx + 1
            }
        }
    }
}

# (2) Main model - restricted to England to compare with model (3) including mdi_eng
df_sel_eng <- df_sel[df_sel$mdi_eng %>% complete.cases, ]
for (outcome in c('dem', 'ad', 'vd')) {
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
                fit <- coxph(model, data = df_sel_eng, method = 'efron')
                name <- str_glue('model2_{outcome}_{ap_unit}_{ap_type}')
                list_fit[[idx]] <- fit
                names(list_fit)[idx] <- name

                # AIC
                list_aic[[idx]] <- AIC(fit)
                names(list_aic)[idx] <- name

                # Proportionality tests
                list_test[[idx]] <- fit %>% cox.zph
                names(list_test)[idx] <- name 

                # Update global index
                idx = idx + 1
            }
        }
    }
}

# (3) Area level covariates: multiple deprivation index quartiles (England)
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
                fit <- coxph(model, data = df_sel_eng, method = 'efron')
                name <- str_glue('model3_{outcome}_{ap_unit}_{ap_type}')
                list_fit[[idx]] <- fit
                names(list_fit)[idx] <- name

                # AIC
                list_aic[[idx]] <- AIC(fit)
                names(list_aic)[idx] <- name

                # Update global index
                idx = idx + 1
            }
        }
    }
}

# (4) Spatial confounding: recruitment centres (Mundlack regression)
# - ap mean within recruitment centre
for (apoll in c('pm25_s', 'pmcoarse_s', 'pmabs_s', 'pm10_s', 'no2_s', 'nox_s', 'appc')) {
    df_sel <- df_sel %>%
        group_by(rec_centre) %>% 
        mutate('{apoll}_mean' := mean(get(apoll), na.rm = T))
}

# - Mundlack regression
# exclude model with all ap; run only for apoll coded as continuous
for (outcome in c('dem', 'ad', 'vd')) {
    for (ap_type in ap_all[-7]) {
        if (ap_type %in% ap_list) { # for single pollutant models
            apoll = str_glue('{ap_type}_s')
        } else { # for appc
            apoll <- ap_type
        }

        # Model
        model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~
                          age + 
                          sex + 
                          ethnic + 
                          ea + 
                          income + 
                          pop_density +
                          {apoll} +
                          {apoll}_mean +
                          ({apoll} + 1 | rec_centre)')
        model <- as.formula(model)

        # Fit and save model
        fit <- coxme(model, data = df_sel)
        name <- str_glue('model4_{outcome}_s_{ap_type}')
        list_fit[[idx]] <- fit
        names(list_fit)[idx] <- name

        # AIC
        list_aic[[idx]] <- AIC(fit)
        names(list_aic)[idx] <- name

        # Update global index
        idx = idx + 1
    }
}

# (5) Model (1) + age/sex/apoll time interactions
# restrict to date_dem (all cause dementia) and appc
# Split data for different time groups; choose cut based on Schoenfeld residual plots (df = 3) 
df_split <- survSplit(Surv(time_dem, status_dem) ~ ., data = df_sel, cut = c(8, 12), episode = 'tgroup')
list_tint <- c('age', 'sex', 'appc') # variables interacting with time
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
    name <- str_glue('model5_dem_s_appc_{var_tint}')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Update global index
    idx = idx + 1
}

# (6) Model (2) + apoll non-linear effect
# do not run for ap quartiles
# for knots use .05 .275 .5 .725 .95 quantiles (Harrell, Springer 2001,
for (outcome in c('dem', 'ad', 'vd')) {
    # Model
    model <- str_glue('Surv(time_{outcome}, status_{outcome}) ~
                      age + 
                      sex + 
                      ethnic + 
                      ea + 
                      income + 
                      pop_density +
                      pspline(appc, df = 3)')
    model <- as.formula(model)

    # Fit and save model
    fit <- coxph(model, data = df_sel, method='efron')
    name <- str_glue('model6_{outcome}_s_appc')
    list_fit[[idx]] <- fit
    names(list_fit)[idx] <- name

    # AIC
    list_aic[[idx]] <- AIC(fit)
    names(list_aic)[idx] <- name

    # Plot spline fit
    list_spline[[idx]] <- termplot(fit, se = T, plot = F)
    names(list_spline)[idx] <- name

    # Update global index
    idx = idx + 1
}

# Save results ----------------------------------------------------------------
save_path <- here(str_glue('output/results/cox/savedata/'))
if (!file.exists(save_path)) {
    dir.create(save_path, recursive = T)
}

saveRDS(list_fit, file = paste0(save_path, 'fits.rds'))
saveRDS(list_test, file = paste0(save_path, 'tests.rds'))
saveRDS(list_aic, file = paste0(save_path, 'aic.rds'))
saveRDS(list_spline, file = paste0(save_path, 'spline.rds'))



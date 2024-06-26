#' Display reults cox analysis
#'
#' Generate tables for results of cox regression; produce plots of Schoenfeld residuals;
#' compare models with and without interactions
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#'
#' ----------------------------------------------------------------------------

# install and import libraries
install.packages('here')
install.packages('finalfit')
install.packages('kableExtra')
install.packages('coxme')
install.packages('stargazer')

library(data.table)
library(tidyverse)
library(here)
library(finalfit)
library(survival)
library(coxme)
library(stargazer)

source(paste0(here(), '/src/forest_plot.R'))

# Load Data -------------------------------------------------------------------
data_path <- str_glue('{here()}/output/results/cox/savedata/')
list_fit_all <- readRDS(file = paste0(data_path, 'fits.rds'))
list_aic_all <- readRDS(file = paste0(data_path, 'aic.rds'))
list_test_all <- readRDS(file = paste0(data_path, 'tests.rds'))
list_spline_all <- readRDS(file = paste0(data_path, 'spline.rds'))

# Discard empty elements in lists
list_fit_all <- list_fit_all[list_fit_all %>% names %>% complete.cases]
list_aic_all <- list_aic_all[list_aic_all %>% names %>% complete.cases]
list_test_all <- list_test_all[list_test_all %>% names %>% complete.cases]
list_spline_all <- list_spline_all[list_spline_all %>% names %>% complete.cases]

# Drop models with all pollutants for forest plots
list_forest <- list_fit_all[!grepl('all', names(list_fit_all))]

# Exclude drop_var from aggregate tables
drop_var <- c('age', 'sex', 'ethnic', '\\bea\\b', 'income', 'pop_density', 'tdi', 'md') %>% 
    paste(collapse = '|')

# Plot and save forest plots --------------------------------------------------
# Aggregate air pollutant HR 
list_aggregate <- as.list(rep(NA, 23))
idx = 1

# Aggregate across pollutants (model 1-17)
for (model in as.character(1:5)) {
    for (outcome in c('dem', 'ad', 'vd')) {
        for (ap_unit in c('s', 'q')) {
            idx_cat <- grep(str_glue('model{model}_{outcome}_{ap_unit}'), names(list_forest))
            if (!is_empty(idx_cat)) {
                # Get aggregated data
                names_cat <- names(list_forest)[idx_cat]
                list_tidy <- map(list_forest[names_cat], finalfit::fit2df, condense = F)
                table_fits <- rbindlist(list_tidy)
                var_mask <- map(table_fits$explanatory, ~ grepl(drop_var, .x)) %>% unlist
                table_fits <- table_fits[!var_mask, ]

                # Save in list
                name <- str_glue('model{model}_{outcome}_{ap_unit}')
                list_aggregate[[idx]] <- table_fits
                names(list_aggregate)[idx] <- name

                # Update global index
                idx = idx + 1

            }
        }
    }
}

# Aggregate across air pollutant units (scaled by IQR and quartiles; model 18-23) 
# Only for model 1 and 2
combine_aptype <- function(x) {
    idx_combine <- grep(x, names(list_aggregate))
    new_order <- order(c(2,3,4, # pm25_q
                         6,7,8, # pmcoarse_q
                         10, 11, 12, # pmabs_q
                         14, 15, 16, # pm10_q
                         18, 19, 20, # no2_q
                         22, 23, 24, # no2_q
                         1, # pm25_s
                         5, # pmcoars_s
                         9, # pmabs_s
                         13, # pm10_s
                         17, # no2_s
                         21, # nox_s
                         25)) # appc
    table_ordered <- rbind(list_aggregate[[idx_combine[2]]], list_aggregate[[idx_combine[1]]])[new_order, ]
    table_ordered
}

name_au <- map(c('1', '2'), 
               function(x) map(c('dem', 'ad', 'vd'), 
                               function(y) str_glue('model{x}_{y}'))) %>% unlist

names(list_aggregate)[17 + 1:6] <- map(1:6, # model >17 save _s and _q combined
                                       function(x) names(list_aggregate)[[17 + x]] <- name_au[x])
list_aggregate[17 + 1:6] <- map(1:6, 
                                function(x) list_aggregate[[17 + x]] <- combine_aptype(name_au[x]))

# Plot forests
save_forest <-  here(str_glue('output/results/cox/display/forests/'))
if (!file.exists(save_forest)) {
    dir.create(save_forest, recursive = TRUE)
}

walk(13:23, # model1-2 with _s and _q combined
     function (x) plot_forest(list_aggregate[[x]],
                              str_glue('{save_forest}{names(list_aggregate)[[x]]}')))

# Generate and save Tables ----------------------------------------------------
save_table <-  here(str_glue('output/results/cox/display/tables/'))
if (!file.exists(save_table)) {
    dir.create(save_table, recursive = TRUE)
}
name_tables <- list_fit_all %>% names
walk(.x = name_tables, ~ finalfit::fit2df(list_fit_all[[.x]], condense = F) %>%
         mutate(# Recode p-value
                'p-value' = case_when(
                                      p < 0.001 ~ '<0.001',
                                      round(p, 2) == .05 ~ as.character(round(p, 3)),
                                      p < .01 ~ str_pad(
                                                        as.character(round(p, 3)),
                                                        width = 4,
                                                        pad = '0',
                                                        side = 'right'
                                                        ),
                              T ~ str_pad(
                                          as.character(round(p, 2)),
                                          width = 4,
                                          pad = '0',
                                          side = 'right')
                              ),
                p = NULL
                ) %>%
         kable(format = 'latex',
               digits = 2,
               booktabs = T) %>%
         kable_styling(latex_options = "striped") %>% 
         kableExtra::save_kable(str_glue('{save_table}{.x}.tex')))

#map(grep('model5', name_tables), ~ ranef(list_fit_all[[.x]])) %>% rbindlist

# Plot splines ----------------------------------------------------------------
attributes(list_spline_all[[1]])
df_spline <- list_spline_all[[1]]$appc
center <- with(df_spline, y[x == min(x)])
ytemp <- exp(df_spline[['y']] + outer(df_spline[['se']], c(0, -1, 1), '*'))
#ytemp <- df_spline[['y']] + outer(df_spline[['se']], c(0, -1, 1), '*')

matplot(df_spline[['x']], (ytemp), type = 'l', lty = c(1, 2, 2), col = 1)
#matplot(df_spline[['x']], (ytemp), log = 'y', type = 'l', lty = c(1, 2, 2), col = 1)


# Plot Shoefield residuals ----------------------------------------------------


#######################################################################

lss %>% length
lss[1]
name_tables[1] 
# Create 
# Loop through:
for (model in as.character(1:(n_model + n_model_ml))) {
    if (idx < n_model * 2 * 6 + 1) { # unilevel models
        for (ap_unit in c('s', 'q')) {
            for (apoll in c('pm_25', 'pm_coarse', 'pm_abs', 'pm_10', 'no_2', 'no_x')) {
                coxfit <- str_glue('fit{model}_{ap_unit}_{apoll}')
                # Results cox fit
                file_name <- str_glue('{save_path}tbl_{coxfit}.tex')
                list_forest[[coxfit]] %>% 
                    finalfit::fit2df() %>%
                    knitr::kable(format = 'latex') %>% 
                    kableExtra::save_kable(file_name)
                # Proportionality check
                file_name <- str_glue('{save_path}ptest_{coxfit}.tex')
                list_test[[coxfit]]$table %>% 
                    knitr::kable(format = 'latex') %>% 
                    kableExtra::save_kable(file_name)
                # Update index
                idx <- idx + 1
            }
        }
    } else { # multilevel models
        for (apoll in c('pm_25', 'pm_coarse', 'pm_abs', 'pm_10', 'no_2', 'no_x')) {
            coxfit <- str_glue('fit{model}_{apoll}')
            # Results cox fit
            file_name <- str_glue('{save_path}tbl_{coxfit}.tex')
            list_forest[[coxfit]] %>% 
                finalfit::fit2df() %>%
                knitr::kable(format = 'latex') %>% 
                kableExtra::save_kable(file_name)
            # Proportionality check
            file_name <- str_glue('{save_path}ptest_{coxfit}.tex')
            list_test[[coxfit]]$table %>% 
                knitr::kable(format = 'latex') %>% 
                kableExtra::save_kable(file_name)
        }
    }
}

# Function for plotting and saving Shoenfield residual
plot_res <- function(x, res = TRUE) {
    data_res <- list_test[[x]]
    # Initialise subplots
    n_subplots <- data_res$table %>% dim %>% .[1] - 1
    n_row <- n_subplots %>% sqrt %>% ceiling
    par(mfrow = c(n_subplots = c(n_row, n_row))) # set grid subplots
    # Plot
    for (ax in 1:n_subplots) {
        plot(data_res[ax], col = 'red', df = 3, resid = res); plot1 <- recordPlot()
    }
    # Save plot
    setEPS()
    postscript(paste0(save_path, x, '.eps'))
    print(plot1) # save prints environment; without print(plot) env is empty
    dev.off()
}

# Schoenfield residual plots for selected models
# Time axis is survival_function(t) using Keplan-Meier method, and transforms
# equally spaced time points in equally spaced survival probabilities, thus 
# evenly distributing events across x axis.
plot_res('fit1_q_pm_25', res = F)
plot_res('fit1_s_pm_25')
plot_res('fit2_s_pm_25')
plot_res('fit3_s_pm_25')
plot_res('fit4_s_pm_25')
plot_res('fit5_s_pm_25')
plot_res('fit6_s_pm_25')

# Merge apoll results in one df - keep only apoll coef
filter_fit <- funciont(x, idx) {
    list_x_f = as.list(rep(NA, 6))
    idx_ap = 1
    for (apoll in c('pm_25', 'pm_coarse', 'pm_abs', 'pm_10', 'no_2', 'no_x')) {
             list_x_f[[idx_ap]] <- list_fit_all[[str_glue('x_{apoll}')]] %>%
                 finalfit::fit2df() %>%
                 .[mask[idx], ]
             names(list_x_f)[idx_ap] = apoll
             idx_ap = idx_ap + 1
    }
    list_x_f
}
2*(list_fit_all[['fit1_s_pm_25']]$loglik - list_fit_all[['fit1_q_pm_25']]$loglik)[2]
                         
    x_m <- list_fit_all[[x]]
for (model in as.character(1:(n_model + n_model_ml))) {
    if (idx < n_model * 2 + 1) { # unilevel models
        for (ap_unit in c('s', 'q')) {
            coxfit <- str_glue('fit{model}_{ap_unit}')
            # Merge
            df_pm_25 <- > list_fit_all[[str_glue('{coxfit}_pm_25')]] %>% 
                finalfit::fit2df()
            assign(df_name,
                   rbind(list_fit_all[[str_glue('{coxfit}_pm_25')]][mask, ],
                         get(str_glue('{coxfit}_pm_coarse'))[mask, ],
                         get(str_glue('{coxfit}_pm_abs'))[mask, ],
                         get(str_glue('{coxfit}_pm_10'))[mask, ],
                         get(str_glue('{coxfit}_no_2'))[mask, ],
                         get(str_glue('{coxfit}_no_x'))[mask, ]
                   )
            )
            # Results cox fit
            file_name <- str_glue('{save_path}tbl_{coxfit}.tex')
            list_fit_all[[coxfit]] %>% 
                finalfit::fit2df() %>%
                knitr::kable(format = 'latex') %>% 
                kableExtra::save_kable(file_name)
            # Update index
            idx <- idx + 1
        }
    } else { # multilevel models
        coxfit <- str_glue('fit{model}_{apoll}')
        # Results cox fit
        file_name <- str_glue('{save_path}tbl_{coxfit}.tex')
        list_fit_all[[coxfit]] %>% 
            finalfit::fit2df() %>%
            knitr::kable(format = 'latex') %>% 
            kableExtra::save_kable(file_name)
    }
}

for (model in as.character(1:n_model)) {
    for (apoll_unit in c('s', 'q')) {
        # Mask to keep only apoll coefficients
        if (model == '1') {
            mask = c(-1:-2)
        } else if (model == '2') {
            mask = c(-1:-9)
        } else if (model == '3') {
            mask = c(-1:-12)
        }
        # Merge
        df_name = str_glue('df{model}_{apoll_unit}')
        assign(df_name,
               rbind(get(str_glue('{df_name}_pm_25'))[mask, ],
                     get(str_glue('{df_name}_pm_coarse'))[mask, ],
                     get(str_glue('{df_name}_pm_abs'))[mask, ],
                     get(str_glue('{df_name}_pm_10'))[mask, ],
                     get(str_glue('{df_name}_no_2'))[mask, ],
                     get(str_glue('{df_name}_no_x'))[mask, ]
                     )
               )
    }
}

# Save merged df as tables
for (model in c('1', '2', '3')) {
    get(str_glue('df{model}_s')) %>% 
        knitr::kable(format = 'latex') %>% 
        kableExtra::save_kable(str_glue('{save_path}tbl_fit{model}_s.tex'))
    get(str_glue('df{model}_q')) %>% 
        knitr::kable(format = 'latex') %>% 
        kableExtra::save_kable(str_glue('{save_path}tbl_fit{model}_q.tex'))
}

# Model 1 - covariate age, sex
#dependent <- 'Surv(time, status)'
#explanatory1 <- c('age', 'sex')
#formula1 <- 'Surv(time, status) ~ age + sex'
#for (apoll in c('pm_25', 'pm_coarse', 'pm_abs', 'pm_10', 'no_2', 'no_x')) {
#    # fit with apoll as quartiles
#    fit <- df %>%
#        finalfit(dependent, append(explanatory1, paste0(apoll, '_q')))
#    assign(paste0('fit1_q_', apoll), fit)
#    # fit with apoll scaled by IQR
#    fit <- df %>% 
#        finalfit(dependent, append(explanatory1, paste0(apoll, '_s')))
#    assign(paste0('fit1_s_', apoll), fit)
#}

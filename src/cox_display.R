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

install.packages('flextable') # needs apt install libfontconfig1-dev and libcairo2-dev
install.packages('gdtools') # needs apt install libfontconfig1-dev and libcairo2-dev
install.packages('interactionR') # needs apt install libfontconfig1-dev and libcairo2-dev
library(interactionR) 

library(data.table)
library(plyr)
library(tidyverse)
library(here)
library(finalfit)
library(kableExtra)
library(survival)
library(survminer)
library(coxme)
library(stargazer)
library(glmnet)
library(xtable)
library(rlang)
library(vctrs)

source(paste0(here(), '/src/forest_plot.R'))

# Load Data -------------------------------------------------------------------
data_path <- str_glue('{here()}/output/results/cox/savedata/')

list_fit <- readRDS(file = paste0(data_path, 'fits.rds'))
list_aic <- readRDS(file = paste0(data_path, 'aic.rds'))
list_tests <- readRDS(file = paste0(data_path, 'tests.rds'))
list_vif <- readRDS(file = paste0(data_path, 'vif.rds'))
list_cca <- readRDS(file = paste0(data_path, 'cca.rds'))

# Discard empty elements in lists
list_fit <- list_fit[list_fit %>% names %>% complete.cases]
list_aic <- list_aic[list_aic %>% names %>% complete.cases]
list_tests <- list_tests[list_tests %>% names %>% complete.cases]
list_vif <- list_vif[list_vif %>% names %>% complete.cases]
list_cca <- list_cca[list_cca %>% names %>% complete.cases]

# Drop models with all pollutants for forest plots and negative control later
list_aggregate <- list_fit[!grepl('all', names(list_fit))]


# Plot and save forest plots --------------------------------------------------
# Exclude drop_var from aggregate tables
drop_var <- c('age', 'sex', 'ethnic', '\\bea\\b', 'income', 'pop_density', 'tdi', 'mdi', 'noise_poll_s') %>% 
    paste(collapse = '|')

# Aggregate air pollutant HR  
list_forest <- as.list(rep(NA, 100))
idx = 1

# Aggregate across pollutants
# model 1:  basic
# model 2:  (1) + mdi_eng
# model 3:  (2) + address_buff = 5
# model 9:  (1) but for copd positive control
# model 10: (9) + mdi_eng
# model 11: (10) + address_buff = 5
# model 18: (3) + noise pollution
for (model in as.character(c(1:3, 9:11, 18))) {
    for (outcome in c('dem', 'ad', 'vd', 'copd')) {
        for (ap_unit in c('s', 'q')) {
            idx_cat <- grep(str_glue('model{model}_{outcome}_{ap_unit}'), names(list_aggregate))
            if (!is_empty(idx_cat)) {
                # Get aggregated data
                names_cat <- names(list_aggregate)[idx_cat]
                list_tidy <- map(list_aggregate[names_cat], finalfit::fit2df, condense = F)
                table_fits <- rbindlist(list_tidy)
                var_mask <- map(table_fits$explanatory, ~ grepl(drop_var, .x)) %>% unlist
                table_fits <- table_fits[!var_mask, ]

                # Save in list
                name <- str_glue('model{model}_{outcome}_{ap_unit}')
                list_forest[[idx]] <- table_fits
                names(list_forest)[idx] <- name

                # Update global index
                idx = idx + 1

            }
        }
    }
}

# Aggregate across air pollutant units (scaled by IQR and quartiles)
# - helper function
combine_aptype <- function(x) {
    idx_combine <- grep(x, names(list_forest))
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
                         21, # no_s
                         25)) # appc
    table_ordered <- rbind(list_forest[[idx_combine[2]]], list_forest[[idx_combine[1]]])[new_order, ]
    table_ordered
}

# - new names (all units)
name_au_dem <- map(as.character(c(1:3)), 
    function(x) map(c('dem', 'ad', 'vd'), 
        function(y) str_glue('model{x}_{y}'))) %>% unlist
name_au_copd <- map(as.character(c(9:11)), 
    function(x) map(c('copd'), 
        function(y) str_glue('model{x}_{y}'))) %>% unlist
name_au_noise <- 'model18_dem'

name_au <- c(name_au_dem, name_au_copd, name_au_noise)

# - add combined models
names(list_forest)[27 + 1:13] <- map(1:13, # model >17 save _s and _q combined
                                       function(x) names(list_forest)[[27 + x]] <- name_au[x])
list_forest[27 + 1:13] <- map(1:13, 
                                function(x) list_forest[[27 + x]] <- combine_aptype(name_au[x]))

list_forest <- list_forest[list_forest %>% names %>% complete.cases]

# - add 'all' models (with all single pollutants in one model) to list_forest
list_fit_all <- list_fit[grepl('all', names(list_fit)) & !grepl('fish', names(list_fit))]
list_fit_all <- map(list_fit_all, finalfit::fit2df, condense = F)
list_fit_all <- map(list_fit_all, ~ .x[!grepl(drop_var, .x$explanatory), ])
list_forest <- c(list_forest, list_fit_all)

# Saving directory
save_forest <-  here(str_glue('output/results/cox/display/forests/'))
if (!file.exists(save_forest)) {
    dir.create(save_forest, recursive = TRUE)
}

# Generate and save forests
walk(27:length(list_forest), # model1-2 with _s and _q combined
     function (x) plot_forest(list_forest[[x]],
        str_glue('{save_forest}{names(list_forest)[[x]]}')))


# # Compare results for global and restricted pollution scores
# list_fit %>% names %>% head(60) 
# list_fit[['model1_dem_s_pmcoarse']]
# list_fit[['model1_dem_s_pm10']]
# list_fit[['model1_dem_s_appc_all']] %>% finalfit::fit2df()
# list_fit[['model1_dem_s_appc']] %>% finalfit::fit2df()


# Tables for cox models -------------------------------------------------------
# Saving directory
save_table_cox <-  here(str_glue('output/results/cox/display/tables_cox/'))
if (!file.exists(save_table_cox)) {
    dir.create(save_table_cox, recursive = TRUE)
}

# Generate tables for cox models (excluding time interaction models)
list_tables <- c('model1_',  # main model
                 'model2_',  # adjusted for mdi
                 'model3_',  # +5
                 'model9_',  # positive control
                 'model10_', # adjusted for mdi
                 'model11_', # +5
                 'model18_', # noise pollution
                 'model19_', # interaction analysis
                 'model20_') %>% paste(collapse = '|') # +5
list_tables <- grep(list_tables, list_fit %>% names, value = T)

walk(.x = list_tables, ~ finalfit::fit2df(list_fit[[.x]], condense = F) %>%
    mutate(# Recode p-value
        'p-value' = case_when(
            p < 0.001 ~ '$<0.001$',
            p >= 0.001 ~ str_pad(
                as.character(round(p, 3)),
                width = 4,
                pad = '0',
                side = 'right') %>% paste0('$', ., '$'),
            TRUE ~ str_pad(
                as.character(round(p, 2)),
                width = 4,
                pad = '0',
                side = 'right') %>% paste0('$', ., '$')
        ),
        p = NULL,
        'explanatory' = case_when(
           grepl('age', explanatory) ~ c('age'),
           grepl('sexmale', explanatory) ~ c('male'),
           grepl('ethnicother', explanatory) ~ c('ethnic. other'), 
           grepl('ea', explanatory) ~ c('education'), 
           grepl('income18,000-30,999', explanatory) ~ c('inc. 18-31k'),
           grepl('income31,000-51,999', explanatory) ~  c('inc. 31-52k'),
           grepl('income52,000-100,000', explanatory) ~  c('inc. 52-100k'), 
           grepl('income>100,000', explanatory) ~ c('inc. $>$ 100k'), 
           grepl('pop_densityrural', explanatory) ~ c('rural area'), 
           grepl('mdi_eng_q2nd:', explanatory) ~ c('mdi 2q int.'), # set before mdi/air pollutants
           grepl('mdi_eng_q3rd:', explanatory) ~ c('mdi 3q int.'), # set before mdi/air pollutants
           grepl('mdi_eng_q4th:', explanatory) ~ c('mdi 4q int.'), # set before mdi/air pollutants
           grepl('mdi_eng_q2nd', explanatory) ~ c('mdi 2q'), 
           grepl('mdi_eng_q3rd', explanatory) ~ c('mdi 3q'), 
           grepl('mdi_eng_q4th', explanatory) ~ c('mdi 4q'), 
           grepl('pm25_s', explanatory) ~ c('PM$_{2.5} iqr$'), 
           grepl('pm25_q2nd', explanatory) ~ c('PM$_{2.5}$ 2q'), 
           grepl('pm25_q3rd', explanatory) ~ c('PM$_{2.5}$ 3q'),
           grepl('pm25_q4th', explanatory) ~ c('PM$_{2.5}$ 4q'), 
           grepl('pmcoarse_s', explanatory) ~ c('PM$_{coarse} iqr$') , 
           grepl('pmcoarse_q2nd', explanatory) ~ c('PM$_{coarse}$ 2q'), 
           grepl('pmcoarse_q3rd', explanatory) ~ c('PM$_{coarse}$ 3q'), 
           grepl('pmcoarse_q4th', explanatory) ~ c('PM$_{coarse}$ 4q'), 
           grepl('pmabs_s', explanatory) ~ c('PM$_{abs} iqr$'), 
           grepl('pmabs_q2nd', explanatory) ~ c('PM$_{abs}$ 2q'), 
           grepl('pmabs_q3rd', explanatory) ~ c('PM$_{abs}$ 3q'), 
           grepl('pmabs_q4th', explanatory) ~ c('PM$_{abs}$ 4q'), 
           grepl('pm10_s', explanatory) ~ c('PM$_{10} iqr$'), 
           grepl('pm10_q2nd', explanatory) ~ c('PM$_{10}$ 2q'), 
           grepl('pm10_q3rd', explanatory) ~ c('PM$_{10}$ 3q'), 
           grepl('pm10_q4th', explanatory) ~ c('PM$_{10}$ 4q'), 
           grepl('no2_s', explanatory) ~ c('NO$_{2} iqr$'), 
           grepl('no2_q2nd', explanatory) ~ c('NO$_{2}$ 2q'), 
           grepl('no2_q3rd', explanatory) ~ c('NO$_{2}$ 3q'), 
           grepl('no2_q4th', explanatory) ~ c('NO$_{2}$ 4q'), 
           grepl('no_s', explanatory) ~ c('NO iqr'), 
           grepl('no_q2nd', explanatory) ~ c('NO 2q'), 
           grepl('no_q3rd', explanatory) ~ c('NO 3q'), 
           grepl('no_q4th', explanatory) ~ c('NO 4q'), 
           grepl('appc', explanatory) ~ c('pollution score'), 
           grepl('appc_all', explanatory) ~ c('pollution score (all)'), 
           grepl('noise_poll_s', explanatory) ~ c('noise pollution')
           # TRUE ~ paste0('$', explanatory, '$')
        ) 
    ) %>%
    kable(format = 'latex',
        digits = 2,
        escape = F,
        booktabs = T) %>%
    kable_styling(latex_options = "striped") %>% 
    kableExtra::save_kable(str_glue('{save_table_cox}{.x}.tex')))


# Tables for Lasso models -----------------------------------------------------
# Saving directory
save_table_lasso <-  here(str_glue('output/results/cox/display/tables_lasso/'))
if (!file.exists(save_table_lasso)) {
    dir.create(save_table_lasso, recursive = TRUE)
}
# 2: Generate tables for lasso cox models
mask_lasso <- grep('lasso', list_fit %>% names, value = T) %>%
    grep('dem', ., value = T)
list_lasso <- list_fit[mask_lasso]

df_lasso <- lmap(list_lasso, \(x) {
         model_ = x[[1]]
         lambda_min = model_[['lambda.min']]
         if (grepl('16', x %>% names)) {
             col_name = '$\\geq1$ year' 
         } else {
             col_name = '$\\geq5$ year'
         }
         df_ = model_ %>% coef(s = lambda_min) %>% exp %>% .[, 1] %>% as.data.frame
         df_ = df_ %>% rename(., !!col_name:= .)
         # df_[['index']] = df_ %>% rownames
         return(df_)
    }
)
df_lasso <- df_lasso  %>% bind_cols() 
df_lasso <- df_lasso %>% slice(8:n())
rownames(df_lasso) <- c('PM$_{2.5}$', 'PM$_{coarse}$', 'PM$_{abs}$', 'PM$_{10}$', 'NO$_{2}$', 'NO') %>%
    as.expression

# Save as tables.tex
df_lasso %>%
    kable(format = 'latex', digits = 4, escape = F, booktabs = T) %>%
    kable_styling(latex_options = "striped") %>% 
    kableExtra::save_kable(str_glue('{save_table_lasso}lasso.tex'))


# Plots/tables for Random effect models -----------------------------------------------------
# Saving directory
save_re <-  here(str_glue('output/results/cox/display/rec_centres/'))
if (!file.exists(save_re)) {
    dir.create(save_re, recursive = TRUE)
}

# Generate tables for random effects
mask_re <- grepl(c('model4|model5'), list_fit %>% names)
list_re <- list_fit[mask_re]

df_rc_slope <- map(c('model4', 'model5'), function(x) {
                    mask_re <- grepl(x, list_fit %>% names)
                    list_re <- list_fit[mask_re] %>% map(ranef)
                    list_slope = list_re %>%
                        map(function(y) {y$rec_centre[, 2] %>% data.frame}) %>%
                        bind_cols
                    names(list_slope) = c('PM$_2.5$', 'PM$_coarse$', 'PM$_abs', 'PM$_10$', 'no$_2$', 'no', 'pollution score')
                    if (grepl('4', list_re %>% names) %>% any) {
                        list_slope['model'] = '>=1 year'
                    } else {
                        list_slope['model'] = '>=5 year'
                    }
                    list_slope
})



# Prepare df for plotting 
new_names <- expression(atop(pollution,score), no, no[2], pm[10], pm[abs], pm[coarse], pm[25]) # attention, axis reversed
fig <- df_rc_slope %>% 
    bind_rows %>% 
    data.table %>% 
    melt(measure.vars = names(.) %>% head(-1)) %>%
    ggplot(aes(x = fct_rev(variable),
               y = exp(value),
               color = model)) +
    geom_boxplot(outlier.shape = NA) +
    scale_x_discrete(name = 'Pollutants',
                     labels = new_names) +
    scale_y_continuous(name = 'Hazard Ratio') +
                     # expand = expansion(mul=0.01)) +
    geom_point(position = position_jitterdodge(), size = 1.5) +
    geom_hline(yintercept = 1, linetype="dotted", linewidth = 0.5) +
    coord_flip() +
    theme_gray(base_size = 15)

# Save plot
fig %>% ggsave(file = paste0(save_re, 'boxplot_rc.eps'))

# Generate tables
df_rc <-  map(c('model4', 'model5'), function(x) {
                  mask_re = grepl(x, list_fit %>% names)
                  list_re = list_fit[mask_re] %>% 
                      map(\(x) {x$vcoef %>% .[[1]] %>% .[2, 2] %>% sqrt %>% round(3)}) %>% 
                      unlist %>% 
                      data.frame %>%
                      rename(., 'sd' = .) %>%
                      mutate( 'sd (HR scale)' = exp(sd))
                  list_re['HR mean'] = list_fit[mask_re] %>% 
                      map(\(x) {x$frail %>% .[[1]] %>% .[, 2] %>% exp %>% mean %>% round(3)}) %>% 
                      unlist %>%
                      data.frame
                  list_re
                  }) %>% bind_rows



# Rename indexes
rc_names <- map(c('$\\geq1$ year', '$\\geq5$ years'), \(x) {
                    map(c('PM$_{2.5}$', 'PM$_{coarse}$', 'PM$_{abs}$', 'PM$_{10}$', 'NO$_{2}$', 'NO', 'pollution score'), \(y) {
                            str_glue('{y} {x}')
                       }) %>% unlist
                     }
) %>% unlist
rownames(df_rc) <- rc_names

# Reorder indexes
n_poll <- nrow(df_rc) / 2
reorder <- rep(c(1:n_poll), each=2)
reorder[seq(2, n_poll * 2, 2)] <- reorder[seq(2, n_poll * 2, 2)] + n_poll
df_rc <- df_rc[reorder,]

# Save as tables.tex
df_rc %>%
    kable(format = 'latex', digits = 4, escape = F, booktabs = T) %>%
    kable_styling(latex_options = "striped") %>% 
    kableExtra::save_kable(str_glue('{save_re}tc_table.tex'))


# Proportional Hazards check results ------------------------------------------------
# - saving directory
save_pha <-  here(str_glue('output/results/cox/display/ph_assumption/'))
if (!file.exists(save_pha)) {
    dir.create(save_pha, recursive = TRUE)
}

# - helper function for Schoenfeld redisuals
plot_schoef <- function(x, res = T, df = 2) {
    data_schoef <- list_tests[[x]]
    # Initialise subplots
    schoef_p <- data_schoef$table[, 'p'] %>% round(3)
    n_subplots <- schoef_p %>% length(.) - 1
    n_row <- n_subplots %>% sqrt %>% ceiling
    par(mfrow = c(n_row, n_row)) # set grid subplots
    # Plot
    for (ax in 1:n_subplots) {
        plot(data_schoef[ax], 
            col = 'red', 
            df = df, 
            resid = res, 
            main = str_glue('p = {schoef_p[ax]}')); plot1 <- recordPlot()
    }
    # Save plot
    if (res) {
        name_schoef = paste0(x, '_res') 
    } else {
        name_schoef = paste0(x, '_nores', '_', df)
    }
    setEPS()
    postscript(paste0(save_pha, name_schoef, '.eps'))
    print(plot1) # save prints environment; without print(plot) env is empty
    dev.off()
}

# Schoenfeld residual plots for selected models
# Time axis is survival_function(t) using Keplan-Meier method, and transforms
# equally spaced time points in equally spaced survival probabilities, thus 
# evenly distributing events across x axis.
for (apoll in c('pm25', 'pmabs', 'pmcoarse', 'pm10', 'no2', 'no', 'all')) {
    plot_schoef(str_glue('model1_dem_s_{apoll}'), res = T)
    plot_schoef(str_glue('model1_dem_s_{apoll}'), res = F, df = 2)
    plot_schoef(str_glue('model1_dem_s_{apoll}'), res = F, df = 3)
}

# Generate and save table for models with time interaction
names_tint <- paste(c('\\bmodel1_dem_s_appc\\b', 'model8_dem_s_appc_age', 'model8_dem_s_appc_sex','model8_dem_s_appc_income'), collapse = '|') 
list_tint <- list_fit[grep(names_tint, list_fit %>% names, value = T)]
df_tint <- lmap(list_tint, \(x) {
                    model_name = x %>% names
                    df_ <- x[[1]] %>% 
                        finalfit::fit2df(condense = F) %>% .[.$explanatory == 'appc', ]
                    df_ = df_ %>%
                        mutate(
                               explanatory = NULL,
                               model = model_name)
                    df_ = df_[, c('model', 'HR', 'L95', 'U95', 'p')]
            }
) %>% rbindlist

df_tint <- df_tint %>% 
    mutate(model = case_when(
                             grepl('model1', model) ~ 'pollution score',
                             grepl('_age', model) ~ 'pollscore-age',
                             grepl('_sex', model) ~ 'pollution score-sex', 
                             grepl('_income', model) ~ 'pollution score-income'),
           'p-value' = case_when(
                         p < 0.01 ~ '$<$0.001',
                         TRUE ~ str_pad(
                                        as.character(round(p, 2)),
                                        width = 4,
                                        pad = '0',
                                        side = 'right')
                         ),
           p = NULL
    )
df_tint %>% 
    kable(format = 'latex', digits = 4, escape = F, booktabs = T) %>%
    kable_styling(latex_options = "striped") %>% 
    kableExtra::save_kable(paste0(save_pha, 'table_appc_tint.tex'))


# Save tables for negative control analyses -----------------------------------
# - saving directory
save_negative <-  here(str_glue('output/results/cox/display/negative_control/'))
if (!file.exists(save_negative)) {
    dir.create(save_negative, recursive = TRUE)
}

# - aggregate across pollutants for models:
# model 11: Negative control with alcohol consumption frequency (without mdi_eng)
# model 12: (11) + mdi_eng
# model 13: (12) + address_buff = 5
list_apoll <- paste(c('pm25', 'pmabs', 'pmcoarse', 'pm10', 'no2', 'no', 'appc'), collapse = '|')
list_negative <- as.list(rep(NA, 20))
idx <- 1
#for (model in as.character(c(12:15))) {
for (model in as.character(c(12:13))) {
    idx_cat <- grep(str_glue('model{model}_fish_s'), names(list_aggregate))
    if (!is_empty(idx_cat)) {
        # Get aggregated data
        names_cat <- names(list_aggregate)[idx_cat]
        list_tidy <- map(list_aggregate[names_cat], \(x) {
                             df_confint = confint.default(x) %>% exp %>% as.data.frame
                             df_confint['OR'] = x$coef %>% exp
                             df_confint})
        list_tidy <- list_tidy %>% reduce(rbind) %>% round(digits = 2)
        colnames(list_tidy) <- c('L95', 'U95', 'OR')
        var_mask <-  list_tidy %>% rownames %>% grep(list_apoll, .)
        table_fits <- list_tidy[var_mask,]

        # Save in list
        name <- str_glue('model{model}_fish_s')
        list_negative[[idx]] <- table_fits
        names(list_negative)[idx] <- name

        # Update global index
        idx = idx + 1
    }
}

list_negative <- list_negative[list_negative %>% names %>% complete.cases]

# Combine in dataframe
df_negative <- lmap(list_negative, function(x) {
                        model_ = x[[1]]
                        if (grepl('12', x %>% names)) {
                            model_name = 'mdi unudjusted' 
                        } else {
                            model_name = 'mdi udjusted'
                        }
                        poll_names = model_ %>% 
                            rownames %>% 
                            map(., function(x) {
                                    paste0(x, '-', substr(model_name, 
                                                          start = 5, 
                                                          stop = 100L)) }) %>% 
                            unlist
                        model_ = model_ %>%
                            mutate(model = model_name,
                                   pollutant = poll_names)
                        model_ = model_ %>% .[, c('model', 'pollutant', 'OR', 'L95', 'U95')]
    }
) %>% bind_rows

# Prepare df indexes and columns
rownames(df_negative) <- df_negative[, 'pollutant']
df_negative <- df_negative  %>% mutate(pollutant = NULL)
new_idx <- vec_rep_each(1:7, 2)
new_idx[seq(2, 14, 2)] <- new_idx[seq(2, 14, 2)] + 7
df_negative <- df_negative[new_idx, ]

df_negative %>% 
    kable(format = 'latex', digits = 4, escape = F, booktabs = T) %>%
    kable_styling(latex_options = "striped") %>% 
    kableExtra::save_kable(paste0(save_negative, 'table_results.tex'))


# Analyse and save interaction model results ----------------------------------
# - saving directory
save_int <-  here(str_glue('output/results/cox/display/interaction/'))
if (!file.exists(save_int)) {
    dir.create(save_int, recursive = TRUE)
}

# Create list of df with interactions
list_interaction <- list_fit[grep(paste(c('model19', 'model20'), collapse = '|'), list_fit %>% names)]
interactions_ <- map(list_interaction, function(x) { 
                         coef_fit = x %>% coef
                         ap_exposure = coef_fit %>% names %>% .[length(coef_fit) - 3]
                         int_ = interactionR(x, exposure_names = c('mdi_eng_q4th', ap_exposure))
                         int_ = int_$dframe[c(1, 2, 3, 4, 7, 8),]
                         return(int_)
            }
)

# Create names for saving
model_names <- map(list_interaction %>% names, function(x) {
        apoll_ = substr(x, start = 14, stop = 100L)
        if (grepl('19', x)) {
            model_names_ = paste0('1year', apoll_)
        } else {
            model_names_ = paste0('5year', apoll_)
        }
        return(model_names_)}
)

# Save as individual tables
walk(seq(1, length(interactions_)), function(x) {
         interactions_[[x]] %>% 
             kable(format = 'latex', digits = 3, escape = F, booktabs = T) %>%
             kable_styling(latex_options = "striped") %>% 
             kableExtra::save_kable(str_glue('{save_int}{model_names[[x]]}.tex'))
})

# Compare Dementia and COPD results (5 vs 1 years with IMD) -------------------
# Saving directory
save_positive <-  here(str_glue('output/results/cox/display/positive_control/'))
if (!file.exists(save_positive)) {
    dir.create(save_positive, recursive = TRUE)
}

# Make dataframe
diff_dem <- list_forest$model3_dem$HR - list_forest$model2_dem$HR %>% as.data.frame()
colnames(diff_dem)  <-  'dem'

diff_copd <- list_forest$model11_copd$HR - list_forest$model10_copd$HR %>% as.data.frame()
colnames(diff_copd)  <-  'copd'

diff_df <- cbind(diff_copd, diff_dem)

# Plot histograms
hist_plot <- ggplot(diff_df) + 
    geom_histogram(aes(x=dem, fill='dem'), alpha=1) +
    geom_histogram(aes(x=copd, fill='copd'), alpha=0.5) +
    geom_vline(aes(xintercept=mean(dem)), color='deepskyblue2', linetype="dashed", linewidth=1) +
    geom_vline(aes(xintercept=mean(copd)), color='orangered2', linetype="dashed",linewidth=1) +
    scale_fill_discrete(name="Outcome", labels = c('COPD', 'all-cause dementia'))

# Save figure
hist_plot %>% ggsave(file = paste0(save_positive, 'hist_5vs1.eps'), device = cairo_ps)

# Test for equivalence in 5 vs 1 difference across outcomes
wilctest <- wilcox.test(diff_df$copd, diff_df$dem, paired=F)
saveRDS(wilctest, file = paste0(save_positive, 'wilctest.rdata.rds'))


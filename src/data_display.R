#' Display data
#' 
#' Plot distribution of dementia cases and exposure range; create tables; use
#' clea_data function to load without missing values (for exposure, outcome
#' and covariates.
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#'
#' ----------------------------------------------------------------------------

# Install and import libraries ------------------------------------------------
install.packages('here')
install.packages('gtsummary')
install.packages('kableExtra')

library(data.table)
library(tidyverse)
library(here)
library(gtsummary)
library('kableExtra')

# Import helper functions -----------------------------------------------------
source(paste0(here(), '/src/select.R'))

# Data selection parameters ---------------------------------------------------
address_buff <- 1 # buffer time for time at baseline address in years
dem_buff <- 1 # buffer time between baseline and dem incidence in years

# Load and select data --------------------------------------------------------
data_path <- paste0(here(), '/data/processed/df_ukbb.rds')

df_all <- readRDS(file = data_path)
dummy <- select_data(df_all, address_buff = address_buff, dem_buff = dem_buff)
df <- dummy[[1]]
eid_drop <- dummy[[2]]

# Convert dates to occurrance (factor yes or no)
list_sel = c('date_dem', 'date_ad', 'date_vd', 'date_lost_fu', 'date_sleep_disorder',
             'date_hypet', 'date_diab', 'date_atherosc', 'date_stroke', 'date_1cancer',
             'date_asthma', 'date_copd', 'date_mood_disorder', 'date_hearloss', 
             'date_angina', 'date_infarct', 'date_death')

df[, list_sel] = df %>% select(all_of(list_sel)) %>% sapply(function(x) {
                                                                mask_na = x %>% is.na
                                                                x = x %>% as.character
                                                                x[mask_na] = 'no'
                                                                x[!mask_na] = 'yes'
                                                                x
                   })

# Generate tables grouped by air pollution ------------------------------------
# List variables to include
list_sel = c('date_dem', 'date_ad', 'date_vd', 'age', 'sex', 'ethnic', 'ea', 'income',
             'tdi_q', 'pop_density', 'rec_centre', 'alcohol', 'smoking', 'activity_score',
             't_out', 'noise_poll', 'sleep_time', 'prs_ad', 'obese', 'chol',
             'date_lost_fu', 'date_death')
#list_sel = c('date_dem', 'date_ad', 'date_vd', 'age', 'sex', 'ethnic', 'ea', 'income', 'tdi',
#             'tdi_q', 'pop_density', 'rec_centre', 'alcohol', 'smoking', 'activity_score',
#             't_out', 'noise_poll', 'cog_score', 'sleep_time', 'date_sleep_disorder', 
#             'prs_ad', 'obese', 'date_hypet', 'date_diab', 'chol', 'date_atherosc', 
#             'date_stroke', 'date_angina', 'date_infarct', 'date_1cancer', 'date_asthma', 
#             'date_copd','date_mood_disorder', 'date_hearloss', 'date_lost_fu', 'date_death')
list_apoll = c('pm25_q', 'pmcoarse_q', 'pmabs_q', 'pm10_q', 'no2_q', 'nox_q')
list_tbl = vector(mode = 'list', length = length(list_apoll))

# Make tables
for (x in seq_along(list_apoll)) {
    dummy_sel = list_sel %>% append(list_apoll[[x]])
    header = list_apoll[[x]] %>% str_sub(1, -3)
    list_tbl[[x]] = df %>% 
        select(all_of(dummy_sel)) %>% 
        tbl_summary(
                    by = list_apoll[x],
                    statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                          '{min}-{max}'),
                                     all_categorical() ~ '{p}% ({n})'),
                    digits = list(all_continuous() ~ 1,
                                  all_categorical() ~ c(1, 0)),
                    type = list(all_continuous() ~ 'continuous2',
                                all_categorical() ~ 'categorical'),
                    label = list(date_dem ~ 'dementia',
                                 date_ad ~ 'Alzheimer',
                                 date_vd ~ 'vascular dementia',
                                 prs_ad ~ 'prs Alzheimer',
                                 tdi_q ~ 'tdi quartiles',
                                 ethnic ~ 'ethnicity',
                                 ea ~ 'education attainment',
                                 pop_density ~ 'population density',
                                 rec_centre ~ 'assessment centre',
                                 alcohol ~ 'alcohol freq.',
                                 smoking ~ 'smoking status',
                                 activity_score ~ 'activity score',
                                 t_out ~ 'time spent outdoor',
                                 noise_poll ~ 'noise pollution',
                                 chol ~ 'cholestol',
#                                 cog_score ~ 'cognitive score',
#                                 date_sleep_disorder ~ 'sleep disorder',
#                                 date_hypet ~ 'hypertension',
#                                 date_diab ~ 'diabetes',
#                                 chol ~ 'cholestol',
#                                 date_atherosc ~ 'atherosclerosis',
#                                 date_stroke ~ 'stroke',
#                                 date_angina ~ 'angina',
#                                 date_infarct ~ 'infarct',
#                                 date_1cancer ~ 'cancer',
#                                 date_asthma ~ 'asthma',
#                                 date_copd ~ 'copd',
#                                 date_mood_disorder ~ 'mood disorder',
#                                 date_hearloss ~ 'hear loss',
                                 date_lost_fu ~ 'lost follow-up',
                                 date_death ~ 'lost death'),
                    missing_text = 'NA'
                    ) %>%
        modify_spanning_header(all_stat_cols() ~ header) %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        add_n()
}

# Save tables
dir_save <- str_glue('{here()}/output/descriptive/')
if (!dir.exists(dir_save)) {
    dir.create(dir_save)
}

for (x in seq_along(list_tbl)) {
    fname = list_tbl[[x]]['by'] %>% str_sub(1, -3)
    fname = paste0(dir_save, 'tbl_', fname)
    list_tbl[[x]] %>% as_gt() %>% gt::gtsave(filename = str_glue('{fname}.tex'))
}

# Generate tables grouped by dementia incidence -------------------------------
list_dem <- list_sel[1:3]
list_sel <- append(list_apoll, list_sel[-1:-3])
list_tbl = vector(mode = 'list', length = length(list_dem))
list_header <- c('All cause dementia', 'Alzheimer', 'Vascular dementia')

# List of tables grouped by air pollution quartiles
for (x in seq_along(list_dem)) {
    dummy_sel = list_sel %>% append(list_dem[[x]])
    header = list_header[x]
    list_tbl[[x]] = df %>% 
        select(all_of(dummy_sel)) %>% 
        tbl_summary(
                    by = list_dem[x],
                    statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                          '{min}-{max}'),
                                     all_categorical() ~ '{p}% ({n})'),
                    digits = list(all_continuous() ~ 1,
                                  all_categorical() ~ c(1, 0)),
                    type = list(all_continuous() ~ 'continuous2',
                                all_categorical() ~ 'categorical'),
                    label = list(pm25_q ~ 'pm 2.5 (Q)',
                                 pmcoarse_q ~ 'pm coarse (Q)',
                                 pmabs_q ~ 'pm absorbance (Q)',
                                 pm10_q ~ 'pm 10 (Q)',
                                 no2_q ~ 'NO2 (Q)',
                                 nox_q ~ 'NOx (Q)',
                                 prs_ad ~ 'prs Alzheimer',
                                 tdi_q ~ 'tdi quartiles',
                                 ethnic ~ 'ethnicity',
                                 ea ~ 'education attainment',
                                 pop_density ~ 'population density',
                                 rec_centre ~ 'assessment centre',
                                 alcohol ~ 'alcohol freq.',
                                 smoking ~ 'smoking status',
                                 activity_score ~ 'activity score',
                                 t_out ~ 'time spent outdoor',
                                 noise_poll ~ 'noise pollution',
                                 chol ~ 'cholestol',
#                                 cog_score ~ 'cognitive score',
#                                 date_sleep_disorder ~ 'sleep disorder',
#                                 date_hypet ~ 'hypertension',
#                                 date_diab ~ 'diabetes',
#                                 date_atherosc ~ 'atherosclerosis',
#                                 date_stroke ~ 'stroke',
#                                 date_angina ~ 'angina',
#                                 date_infarct ~ 'infarct',
#                                 date_1cancer ~ 'cancer',
#                                 date_asthma ~ 'asthma',
#                                 date_copd ~ 'copd',
#                                 date_mood_disorder ~ 'mood disorder',
#                                 date_hearloss ~ 'hear loss',
                                 date_lost_fu ~ 'lost follow-up',
                                 date_death ~ 'lost death'),
                    missing_text = 'NA'
                    ) %>%
        modify_spanning_header(all_stat_cols() ~ header) %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        add_n()
}

# Save tables
for (x in seq_along(list_tbl)) {
    fname = list_tbl[[x]]['by'] %>% str_sub(6)
    fname = paste0(dir_save, 'tbl_', fname)
    list_tbl[[x]] %>% as_gt() %>% gt::gtsave(filename = str_glue('{fname}.tex'))
}

# Compute relationships between covariates ------------------------------------ 
# APoll correlation matrix
df[, gsub('_q', '', list_apoll)] %>% 
    cor %>% 
    as.data.frame %>% 
    knitr::kable(., digits =2, format = 'latex') %>% 
    kableExtra::save_kable(file = paste0(dir_save, 'cor_apol.tex'))

# Correlation APoll-TDI
df[, c(gsub('_q', '', list_apoll), 'tdi')] %>% 
    na.omit %>% 
    cor %>% 
    as.data.frame %>%
    .[1:6, 7, drop = F] %>%
    knitr::kable(., digits =2, format = 'latex') %>% 
    kableExtra::save_kable(file = paste0(dir_save, 'cor_apol_tdi.tex'))

# APoll vs population density (rural/urban)
df %>% 
    select(all_of(c(list_apoll,
                    gsub('_q', '_s', list_apoll),
                    list_dem,
                    'pop_density'))) %>%
    tbl_summary(by = pop_density,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(date_dem ~ 'dementia',
                             date_ad ~ 'Alzheimer',
                             date_vd ~ 'vascular dementia',
                             pm25_q ~ 'pm 2.5 (Q)',
                             pmcoarse_q ~ 'pm coarse (Q)',
                             pmabs_q ~ 'pm absorbance (Q)',
                             pm10_q ~ 'pm 10 (Q)',
                             no2_q ~ 'NO2 (Q)',
                             nox_q ~ 'NOx (Q)',
                             pm25_s ~ 'pm 2.5 (IQR)',
                             pmcoarse_s ~ 'pm coarse (IQR)',
                             pmabs_s ~ 'pm absorbance (IQR)',
                             pm10_s ~ 'pm 10 (IQR)',
                             no2_s ~ 'NO2 (IQR)',
                             nox_s ~ 'NOx (IQR)'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Population density') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_apoll_popdensity.tex'))

# tdi_q vs ea, income, ethnicity, population density (sec: socioeconomic covariates)
df %>% 
    select(all_of(c('tdi_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = tdi_q,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Townsend deprivation index') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_tdi_sec.tex'))

#  ea vs tdi_q, income, ethnicity, population density (sec: socioeconomic covariates)
df %>% 
    select(all_of(c('tdi_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = ea,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Education attainment') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_ea_sec.tex'))

#  income vs tdi_q, ea, ethnicity, population density (sec: socioeconomic covariates)
df %>% 
    select(all_of(c('tdi_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = income,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Income') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_income_sec.tex'))

#  population density vs tdi_q, ea, ethnicity, income (sec: socioeconomic covariates)
df %>% 
    select(all_of(c('tdi_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = pop_density,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Population density') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_popdensity_sec.tex'))

#  tdi_q vs multiple deprication index (England)
df %>% 
    select(all_of(c('tdi_q',
                    'mdi_eng'))) %>%
    tbl_summary(by = tdi_q,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             mdi_eng ~ 'deprivation idx (England)'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Townsend deprivation index') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_tdi_mdi.tex'))

#  alcohol consumption vs tdi_q, ea, ethnicity, income, population density (sec: socioeconomic covariates)
df %>% 
    select(all_of(c('alcohol',
                    'tdi_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = alcohol,
                statistic = list(all_continuous() ~ c('{mean} ({sd})', 
                                                       '{min}-{max}'),
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ c(1, 0)),
                type = list(all_continuous() ~ 'continuous2',
                            all_categorical() ~ 'categorical'),
                label = list(tdi_q ~ 'tdi quartiles',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Alcohol consumption freq.') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_alcohol_sec.tex'))


# # Generate plots --------------------------------------------------------------
#ATTENTION: changed date_x date: so either put in new script or run before conversion
# But probably not necessary to make these plots anymore
# # Pollution data histogram
# list_fig = vector(mode = 'list', length = length(list_apoll))
# 
# fig_pm25 <- ggplot(df, aes(x = pm25)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# fig_pmabs <- ggplot(df, aes(x = pmabs)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# fig_pmcoarse <- ggplot(df, aes(x = pmcoarse)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# fig_pm10 <- ggplot(df, aes(x = pm10)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# fig_no2 <- ggplot(df, aes(x = NO2)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# fig_nox <- ggplot(df, aes(x = NOx)) +
#     geom_histogram(bins=100) +
#     labs(x = expression(over(paste(mu, g), m^3)))
# # Outcome histogram
# fig_dem <- df %>%
#     drop_na(date_dem) %>%
#     ggplot(., aes(x = date_dem)) +
#     geom_histogram(bins=100) +
#     labs(x = 'time', y = 'incidence')
# fig_AD <- df %>%
#     drop_na(date_AD) %>%
#     ggplot(., aes(x = date_AD)) +
#     geom_histogram(bins=100) +
#     labs(x = 'time', y = 'incidence')
# fig_vd <- df %>%
#     drop_na(date_vd) %>%
#     ggplot(., aes(x = date_vd)) +
#     geom_histogram(bins=100) +
#     labs(x = 'time', y = 'incidence')
# 
# # Save plots ------------------------------------------------------------------
# path_save <- (paste0(here(), '/output'))
# if (dir.exists(path_save) == FALSE) {
#     dir.create(path_save)
# }
# # Exposures
# ggsave(file = paste0(path_save, '/fig_pm25.eps'), plot = fig_pm25)
# ggsave(file = paste0(path_save, '/fig_pmabs.eps'), plot = fig_pmabs)
# ggsave(file = paste0(path_save, '/fig_pmcoarse.eps'), plot = fig_pmcoarse)
# ggsave(file = paste0(path_save, '/fig_pm10.eps'), plot = fig_pm10)
# ggsave(file = paste0(path_save, '/fig_no2.eps'), plot = fig_no2)
# ggsave(file = paste0(path_save, '/fig_nox.eps'), plot = fig_nox)
# # Outcome
# ggsave(file = paste0(path_save, '/fig_dem.eps'), plot = fig_dem)
# ggsave(file = paste0(path_save, '/fig_AD.eps'), plot = fig_AD)
# ggsave(file = paste0(path_save, '/fig_vd.eps'), plot = fig_vd)











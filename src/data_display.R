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
library(xtable)
library(gt)

# Import helper functions -----------------------------------------------------
source(paste0(here(), '/src/select.R'))

# Data selection parameters ---------------------------------------------------
address_buff <- 1 # buffer time for time at baseline address in years
dem_buff <- 0 # buffer time between baseline and dem incidence in years

# Load and select data --------------------------------------------------------
data_path <- paste0(here(), '/data/processed/df_ukbb.rds')
df_all <- readRDS(file = data_path)
list_covariates <- c('age', # for main model (2)
                     'sex', 
                     'ethnic', 
                     'ea',
                     'income', # causes ~ 70k to drop!
                     'pop_density',
                     'mdi_eng')
dummy <- select_data(df_all, 
                     address_buff = address_buff,
                     dem_buff = dem_buff,
                     list_covariates)
df_sel <- dummy[[1]]

# Recode income
levels(df_sel[['income']]) <- c('$<$18k', 'inc. 18-31k', 'inc. 31-52', 'inc. 52-100k', 'inc. $>$100k')

# Convert dates to occurrance (factor yes or no)
list_dates = c('date_dem', 'date_ad', 'date_vd', 'date_lost_fu', 'date_death')

df_sel[, list_dates] = df_sel %>% select(all_of(list_dates)) %>% sapply(function(x) {
                                                                            mask_na = x %>% is.na
                                                                            x = x %>% as.character
                                                                            x[mask_na] = 'no'
                                                                            x[!mask_na] = 'yes'
                                                                            x
                   })


# Save directory --------------------------------------------------------------
dir_save <- str_glue('{here()}/output/descriptive/')
if (!dir.exists(dir_save)) {
    dir.create(dir_save)
}

# Generate tables grouped by air pollution ------------------------------------
# List variables to include
list_sel = c('date_dem', 'date_ad', 'date_vd', 'age', 'sex', 'ethnic', 'ea', 'income',
             'pop_density', 'mdi_eng_q', 'noise_poll', 'date_lost_fu',
             'date_death')
list_apoll = c('pm25_q', 'pmcoarse_q', 'pmabs_q', 'pm10_q', 'no2_q', 'no_q')
list_tbl = vector(mode = 'list', length = length(list_apoll))

# Make tables
for (x in seq_along(list_apoll)) {
    dummy_sel = list_sel %>% append(list_apoll[[x]]) # selecion for each apoll
    header = list_apoll[[x]] %>% str_sub(1, -3)
    list_tbl[[x]] = df_sel %>% 
        select(all_of(dummy_sel)) %>% 
        tbl_summary(
                    by = list_apoll[x],
                    statistic = list(all_continuous() ~ c('{mean} ({sd})'),
                                     all_categorical() ~ '{p}%'),
                    digits = list(all_continuous() ~ 1,
                                  all_categorical() ~ 1),
                    type = list(all_continuous() ~ 'continuous',
                                all_categorical() ~ 'categorical'),
                    label = list(date_dem ~ 'dementia',
                                 date_ad ~ 'Alzheimer',
                                 date_vd ~ 'vascular dementia',
                                 mdi_eng_q~ 'deprivation (q)',
                                 ethnic ~ 'ethnicity',
                                 ea ~ 'education attainment',
                                 pop_density ~ 'population density',
                                 noise_poll ~ 'noise pollution',
                                 date_lost_fu ~ 'lost follow-up',
                                 date_death ~ 'lost death'),
                    missing_text = 'NA'
                    ) %>%
        modify_spanning_header(all_stat_cols() ~ header) %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        add_n()
}


for (x in seq_along(list_tbl)) {
    fname = list_tbl[[x]]$input['by'] %>% str_sub(1, -3)
    fname = paste0(dir_save, 'tbl_', fname)
    list_tbl[[x]] %>% as_gt() %>% gt::gtsave(filename = str_glue('{fname}.tex'))
}

# Generate tables grouped by dementia incidence -------------------------------
list_dem <- list_sel[1:3]
list_apoll = c('pm25_q', 'pmcoarse_q', 'pmabs_q', 'pm10_q', 'no2_q', 'no_q')
list_sel <- append(list_apoll, list_sel[-1:-3])
list_header <- c('All cause dementia', 'Alzheimer', 'Vascular dementia')
list_tbl = vector(mode = 'list', length = length(list_dem))

# List of tables grouped by dementia incidence
for (x in seq_along(list_dem)) {
    dummy_sel = list_sel %>% append(list_dem[[x]])
    header = list_header[x]
    list_tbl[[x]] = df_sel %>% 
        select(all_of(dummy_sel)) %>% 
        tbl_summary(
                    by = list_dem[x],
                    statistic = list(all_continuous() ~ c('{mean} ({sd})'),
                                     all_categorical() ~ '{p}%'),
                    digits = list(all_continuous() ~ 1,
                                  all_categorical() ~ 1),
                    type = list(all_continuous() ~ 'continuous',
                                all_categorical() ~ 'categorical'),
                    label = list( pm25_q ~ 'pm$_{2.5}$ (q)',
                                 pmcoarse_q ~ 'pm$_{coarse}$ (q)',
                                 pmabs_q ~ 'pm$_{abs}$ (q)',
                                 pm10_q ~ 'pm$_{10}$ (q)',
                                 no2_q ~ 'no$_2$ (q)',
                                 no_q ~ 'no (q)',
                                 ethnic ~ 'ethnicity',
                                 ea ~ 'education',
                                 pop_density ~ 'population density',
                                 mdi_eng_q~ 'deprivation (q)',
                                 noise_poll ~ 'noise pollution',
                                 date_lost_fu ~ 'lost follow-up',
                                 date_death ~ 'lost death'),
                    missing_text = 'nan'
                    ) %>%
        modify_spanning_header(all_stat_cols() ~ header) %>%
        modify_header(label ~ 'Characteristics') %>%
        modify_footnote(all_stat_cols() ~ NA) %>%
        bold_labels()
}

summary_dem <- df_sel %>%
    select(all_of(list_sel)) %>%
    tbl_summary(
                statistic = list(all_continuous() ~ c('{mean} ({sd})'),
                                 all_categorical() ~ '{p}%'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list( pm25_q ~ 'pm$_{2.5}$ (q))',
                             pmcoarse_q ~ 'pm$_{coarse}$ (q)',
                             pmabs_q ~ 'pm$_{abs}$ (q)',
                             pm10_q ~ 'pm$_{10}$ (q)',
                             no2_q ~ 'no$_2$ (q)',
                             no_q ~ 'no (q)',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education',
                             pop_density ~ 'population density',
                             mdi_eng_q~ 'deprivation (q)',
                             noise_poll ~ 'noise pollution',
                             date_lost_fu ~ 'lost follow-up',
                             date_death ~ 'lost death'),
                missing_text = 'nan'
                ) %>%
    modify_spanning_header(all_stat_cols() ~ header) %>%
    modify_header(label ~ 'Characteristics') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    bold_labels()

# Save tables
for (x in seq_along(list_tbl)) { # single dem tables
    fname = list_tbl[[x]]$input['by'] %>% str_sub(6)
    fname = paste0(dir_save, 'tbl_', fname)
    list_tbl[[x]] %>% as_gt() %>% gt::gtsave(filename = str_glue('{fname}.tex'))
}
tbl_merge(tbls = append(list(summary_dem), list_tbl), tab_spanner = c('overall', list_header)) %>% # all dem table
    as_gt() %>% 
    gt::gtsave(filename = str_glue('{dir_save}tbl_dem_all.tex'))

# Compute relationships between covariates ------------------------------------ 
# APoll correlation matrix
df_sel[, gsub('_q', '', list_apoll)] %>% 
    cor %>% 
    as.data.frame %>% 
    knitr::kable(., digits =2, format = 'latex') %>% 
    kableExtra::save_kable(file = paste0(dir_save, 'cor_apol.tex'))

# Correlation APoll-mdi_eng
df_sel[, c(gsub('_q', '', list_apoll), 'mdi_eng')] %>% 
    na.omit %>% 
    cor %>% 
    as.data.frame %>%
    .[1:6, 7, drop = F] %>%
    knitr::kable(., digits =2, format = 'latex') %>% 
    kableExtra::save_kable(file = paste0(dir_save, 'cor_apol_mdi_eng.tex'))

# APoll vs population density (rural/urban)
df_sel %>% 
    select(all_of(c(list_apoll,
                    gsub('_q', '_s', list_apoll),
                    list_dem,
                    'pop_density'))) %>%
    tbl_summary(by = pop_density,
                statistic = list(all_continuous() ~ '{mean} ({sd})', 
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(date_dem ~ 'dementia',
                             date_ad ~ 'Alzheimer',
                             date_vd ~ 'vascular dementia',
                             pm25_q ~ 'pm 2.5 (Q)',
                             pmcoarse_q ~ 'pm coarse (Q)',
                             pmabs_q ~ 'pm absorbance (Q)',
                             pm10_q ~ 'pm 10 (Q)',
                             no2_q ~ 'NO2 (Q)',
                             no_q ~ 'NO (Q)',
                             pm25_s ~ 'pm 2.5 (IQR)',
                             pmcoarse_s ~ 'pm coarse (IQR)',
                             pmabs_s ~ 'pm absorbance (IQR)',
                             pm10_s ~ 'pm 10 (IQR)',
                             no2_s ~ 'NO2 (IQR)',
                             no_s ~ 'NO (IQR)'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Population density') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_apoll_popdensity.tex'))

# mdi_eng_q vs ea, income, ethnicity, population density (sec: socioeconomic covariates)
df_sel %>% 
    select(all_of(c('mdi_eng_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = mdi_eng_q,
                statistic = list(all_continuous() ~ '{mean} ({sd})', 
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(
                             income ~ 'income',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Townsend deprivation index') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_mdi_sec.tex'))

#  ea vs mdi_eng_q, income, ethnicity, population density (sec: socioeconomic covariates)
df_sel %>% 
    select(all_of(c('mdi_eng_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = ea,
                statistic = list(all_continuous() ~ '{mean} ({sd})', 
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(mdi_eng_q ~ 'mdi quartiles',
                             income ~ 'income',
                             ethnic ~ 'ethnicity',
                             pop_density ~ 'population density'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Education attainment') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_ea_sec.tex'))

#  income vs mdi_eng_q, ea, ethnicity, population density (sec: socioeconomic covariates)
df_sel %>% 
    select(all_of(c('mdi_eng_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = income,
                statistic = list(all_continuous() ~ '{mean} ({sd})', 
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(mdi_eng_q ~ 'deprivation (q)',
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

#  population density vs mdi_eng_q, ea, ethnicity, income (sec: socioeconomic covariates)
df_sel %>% 
    select(all_of(c('mdi_eng_q',
                    'ea',
                    'income',
                    'ethnic',
                    'pop_density'))) %>%
    tbl_summary(by = pop_density,
                statistic = list(all_continuous() ~ '{mean} ({sd})', 
                                 all_categorical() ~ '{p}% ({n})'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(mdi_eng_q ~ 'deprivation (q)',
                             income ~ 'income',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment'),
                missing_text = 'NA'
                ) %>% 
    modify_spanning_header(all_stat_cols() ~ 'Population density') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    add_n() %>%
    as_gt() %>% 
    gt::gtsave(filename = paste0(dir_save, 'tbl_popdensity_sec.tex'))

# IQR and quartile values for air pollutants
ap_list_all <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'no')
info_apoll <- map(ap_list_all, function(x) {
                      c(q = quantile(df_sel[, x], probs = 0:4/4, na.rm = T,),
                        IQR = IQR(df_sel[, x], na.rm = T))
                }) %>% setNames(ap_list_all) %>% as.data.frame %>% t %>% as.data.frame
info_apoll %>% xtable %>% print(file = paste0(dir_save, 'apoll_info.tex'))

# Characteristics of participants who lived <5 years at baseline
df_address <- df_all %>% 
    drop_na(pm25) %>%
    mutate(address5 = case_when(time_cuaddress <5 ~ 'less than 5 years',
                                time_cuaddress >= 5 ~ '5 years or more'))

df_address %>% 
    select('address5', 'age', 'sex', 'ethnic', 'ea', 'income', 'pop_density', 'mdi_eng_q') %>%
    tbl_summary(by = address5,
                statistic = list(all_continuous() ~ c('{mean} ({sd})'), 
                                 all_categorical() ~ c('{p}%'),
                digits = list(all_continuous() ~ 1,
                              all_categorical() ~ 1)),
                type = list(all_continuous() ~ 'continuous',
                            all_categorical() ~ 'categorical'),
                label = list(mdi_eng_q~ 'deprivation (q)',
                             ethnic ~ 'ethnicity',
                             ea ~ 'education attainment',
                             pop_density ~ 'population density'),
                missing = 'no') %>%
    modify_spanning_header(c(stat_1, stat_2) ~ 'Time at baseline address') %>%
    modify_header(label ~ 'Characteristics') %>%
    modify_footnote(all_stat_cols() ~ NA) %>%
    bold_labels() %>%
    as_gt() %>%
    tab_stub_indent(rows = everything(), indent = 3) %>%
    gt::gtsave(filename = paste0(dir_save, 'tbl_address5.tex'))


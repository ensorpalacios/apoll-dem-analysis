#' Preprocess data
#'
#' Load, clean and recode data.
#' categorical covariates coded as nominal except alcohol consumption frequency
#' (used as outcome in negative control analysis)
#'
#' @author Ensor Palacios, email{ensorrafael.palacios@bristol.ac.uk}
#' @date 2024-05-30

# Shebang ---------------------------------------------------------------------
#!/usr/loca/bin/Rscript

# Install and import libraries ------------------------------------------------
install.packages('here')
library(data.table)
library(tidyverse)
library(here)

# Load data -------------------------------------------------------------------
path_data <- paste0(here(), '/data/raw/data_participant.csv')
df <- fread(path_data) %>% as.data.frame

# Rename columns --------------------------------------------------------------
var_list = c('eid',  '24006-0.0', '24007-0.0', '24008-0.0', '24005-0.0', '24003-0.0', '24004-0.0', '42018-0.0', '42019-0.0', '42020-0.0', '42021-0.0', '42022-0.0', '42023-0.0',
             '21022-0.0', '31-0.0',
             '53-0.0', '54-0.0', '699-0.0', '129-0.0', '130-0.0',
             '21000-0.0', '6138-0.0', '845-0.0', '738-0.0', '22189-0.0', '20118-0.0',
             '26410-0.0', '26426-0.0', '26427-0.0',
             '20116-0.0', '1558-0.0',
             '22037-0.0', '22038-0.0', '22039-0.0',
             '1060-0.0', '1050-0.0',
             '24024-0.0',
             '1160-0.0', '131060-0.0', '26206-0.0', '21001-0.0', '131286-0.0', '130706-0.0', '130708-0.0', '30690-0.0', '131380-0.0',
             #'131360-0.0', '131362-0.0', '131364-0.0', '131366-0.0', '131368-0.0', '131370-0.0', '131372-0.0', '131374-0.0', '131376-0.0', '131378-0.0',
             '42006-0.0',
             '131296-0.0', '42000-0.0',
             '40005-0.0',
             '42014-0.0', '42016-0.0', # Attention: in new dataset using 42014 (arithmetically defined asthma); old code is '131494-0.0'
             '130892-0.0', '130894-0.0',
             '131258-0.0',
             '20023-0.0', '20018-0.0', '399-0.1', '399-0.2', '399-0.3', '20016-0.0',
             '40000-0.0', '191-0.0', '190-0.0',
             '1349-0.0', '1329-0.0'
)
var_names = c('eid', 'pm25', 'pmabs', 'pmcoarse', 'pm10', 'no2', 'nox', 'date_dem', 'src_dem', 'date_ad', 'src_ad', 'date_vd', 'src_vd',
              'age', 'sex',
              'date_baseline', 'rec_centre','time_cuaddress', 'biloc_north', 'biloc_east',
              'ethnic', 'edu', 'age_edu', 'income', 'tdi', 'pop_density',
              'mdi_eng', 'mdi_wales', 'mdi_scot',
              'smoking', 'alcohol',
              'met_walking', 'met_mod_activity', 'met_vig_activity',
              't_out_winter', 't_out_summer',
              'noise_poll',
              'sleep_time', 'date_sleep_disorder', 'prs_ad', 'bmi', 'date_hypet', 'date_diab_ins', 'date_diab_noins', 'chol', 'date_atherosc',
              #'date_cevd_0', 'date_cevd_1', 'date_cevd_2', 'date_cevd_3', 'date_cevd_4', 'date_cevd_5', 'date_cevd_6', 'date_cevd_7', 'date_cevd_8', 'date_cevd_9',
              'date_stroke',
              'date_angina', 'date_infarct',
              'date_1cancer',
              'date_asthma', 'date_copd',
              'date_bipolar', 'date_depress',
              'date_hearloss',
              'reaction_time', 'prosp_memory', 'pairs_match1', 'pairs_match2', 'pairs_match3', 'fluid_intelligence',
              'date_death', 'date_lost_fu', 'reason_lost_fu',
              'processed_meat', 'fish'
)

name_keys <- data.frame('ukbb' = var_list, 'mine' = var_names)
names(df)[match(name_keys$ukbb, names(df))] <- name_keys$mine

# Helper functions ------------------------------------------------------------
# Clean dates
clean_dates <- function(x) {
    x = case_match(x,
                   c(as.IDate('1901-01-01'), # date before birth
                     as.IDate('1902-02-02'), # date = birth date
                     as.IDate('1903-03-03'), # date same year of birth
                     as.IDate('1909-09-09'), # date in future (1)
                     as.IDate('2037-07-07')  # date in future (2)
                     ) ~ NA,
                   .default = x
    )
    as.Date(x)
}

# Factorise x into quartiles
factorise <- function(x) {
    x_q <- cut(x,
               quantile(x,
                        probs = 0:4/4, 
                        na.rm = T
                        ), 
               include.lowest = T,
               labels = c('1st', '2nd', '3rd', '4th'),
               ordered_result = F)
    x_q
}

# Rescale x by IQR
rescale <- function(x) {
    x_s <- (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
    x_s
}

# Clean/recode data -----------------------------------------------------------
# Extract nitric oxide (no) values - no_x - no_2 (Cyrys et al., 2012)
df <- df %>% mutate(no = nox - no2)

# Transform air pollution into quartiles/rescale by IQR
ap_list_all <- c('pm25', 'pmcoarse', 'pmabs', 'pm10', 'no2', 'no')
for (apoll in ap_list_all) {
         df[, paste0(apoll, '_q')] <- factorise(df[, apoll])
         df[, paste0(apoll, '_s')] <- rescale(df[, apoll])
}

# PCA on apoll_s; retain pc1 only; one with all apoll inlcuded, one with
# pmcoarse and pm10 excluded post-hoc based on results from single pollutant
# models (don't show any association with dementia)
pc_apoll_all <- df %>% # pca with all ap
    select(all_of(ap_list_all)) %>%
    drop_na %>%
    prcomp(center = T, scale. = T)
df['appc_all'] <- df %>% 
    select(all_of(ap_list_all)) %>%
    as.matrix(.) %*% pc_apoll_all[['rotation']][, 'PC1'] %>% .[, 1] * -1

ap_list <- c('pm25_s', 'pmabs_s', 'no2_s', 'no_s')
pc_apoll <- df %>% # pca with selected ap
    select(all_of(ap_list)) %>%
    drop_na %>%
    prcomp(center = T, scale. = T)
df['appc'] <- df %>% 
    select(all_of(ap_list)) %>%
    as.matrix(.) %*% pc_apoll[['rotation']][, 'PC1'] %>% .[, 1] * -1

# Clean time of dementia
df[, 'date_dem'] <- clean_dates(df[, 'date_dem'])
df[, 'date_ad'] <- clean_dates(df[, 'date_ad'])
df[, 'date_vd'] <- clean_dates(df[, 'date_vd'])

# Recode source dementias
df[, 'src_dem'] <- factor(df[, 'src_dem'],
                          levels = c(0, 11, 12, 21, 22),
                          labels = c('self-rep', 'hosp1', 'death1', 'hosp2', 'death2')
)
df[, 'src_ad'] <- factor(df[, 'src_ad'],
                         levels = c(0, 11, 12, 21, 22),
                         labels = c('self-rep', 'hosp1', 'death1', 'hosp2', 'death2')
)
df[, 'src_vd'] <- factor(df[, 'src_vd'],
                         levels = c(0, 11, 12, 21, 22),
                         labels = c('self-rep', 'hosp1', 'death1', 'hosp2', 'death2')
)

# Transform TDI into quartiles
df[, 'tdi_q'] <- factorise(df[, 'tdi'])

# Take log mdi (deprivation index) and combine across Englang, Whales and Scotland
df[, 'mdi_eng_q'] <- factorise(df[, 'mdi_eng'])

# Strata for sex (adult < 60; old >= 60)
df[, 'age_strata'] <- ifelse(df[, 'age'] >= 60, '+60', '<60') %>% factor

# Recode sex
df[, 'sex'] <- factor(df[, 'sex'], levels = c(0, 1), labels = c('female', 'male'))

# Recode ethnicity
df[, 'ethnic'] <- factor(df[, 'ethnic'],
                         levels = c(1, 1001, 2001, 3001, 4001,
                                    2, 1002, 2002, 3002, 4002,
                                    3, 1003, 2003, 3003, 4003,
                                    4, 2004, 3004, 5, 6,
                                    -1 # do not know (assume they are 'other')
                                    ),
                  labels = c('white', 'white', 'other', 'other', 'other',
                             'other', 'white', 'other', 'other', 'other',
                             'other', 'white', 'other', 'other', 'other',
                             'other', 'other', 'other', 'other', 'other',
                             'other'
                             ),
                  exclude = -3 # prefer not to answer
)

# Transform prs_ad into quartiles
df[, 'prs_ad_q'] <- factorise(df[, 'prs_ad'])

# Recode reason lost follow-up
df[, 'reason_lost_fu'] <- factor(df[, 'reason_lost_fu'],
                                 levels = c(1, 2, 3, 4, 5),
                                 labels = c('death', 'lost_fu', 'left_UK', 'left_UK', 'withdrawn')
)

# Create education score (from Okbay et al., Nat Genetics, 2022)
# 1) Compute education attainment (ea) for NVQ/HND/HNC using age left education:
# followed Okbay et al., but additionally:
# capp age_edu to 24 so ea_nvq < ea_college
# substitute 'do not know', 'prefer not to answer' and NA with average ea_nvq
# recode never went to school with NA
mask_nvq <- grepl('5', df[, 'edu'])
ea_nvq <- df[mask_nvq , 'age_edu']
ea_nvq[ea_nvq > 24] <- 24
ea_nvq <- case_match(ea_nvq,
                     c(-1, -3, NA) ~ ea_nvq[ea_nvq > 0] %>%
                         mean(na.rm = T) %>% 
                         round(0), # do not know, prefer not to answer, NA
                     -2 ~ NA,      # never went to school
                     .default = ea_nvq
)
ea_nvq <- ea_nvq - 5 
# 2) Recode qualification with highest ea
ea <- rep(NA, dim(df)[1])          # no response
ea[grepl('-7', df[, 'edu'])] <- 7  # none of the below
ea[grepl('4', df[, 'edu'])]  <- 10 # CSE
ea[grepl('3', df[, 'edu'])]  <- 10 # O levels/GCSEs
ea[grepl('2', df[, 'edu'])]  <- 13 # A levels/AS levels
ea[grepl('6', df[, 'edu'])]  <- 15 # other professional qualifications
ea[grepl('1', df[, 'edu'])]  <- 20 # College/University degree
# 3) Add ea_nvq if it's highest and add to df
ea[grepl('5', df[, 'edu'])] <- pmax(ea[grepl('5', df[, 'edu'])],
                                    ea_nvq,
                                    na.rm = T)
df['ea'] <- ea

# Recode household income before tax
df[, 'income'] <- factor(df[, 'income'],
                         levels = c(1,
                                    2,
                                    3,
                                    4,
                                    5
                                    ),
                  labels = c('<18,000',
                             '18,000-30,999',
                             '31,000-51,999',
                             '52,000-100,000',
                             '>100,000'
                             ),
                  exclude = c(-1, -3), # do not know, prefer not to answer
                  ordered = F
)

# Recode population density into urban and rural areas; followed
# ukbb definitions
df[, 'pop_density'] <- factor(df[, 'pop_density'],
                              levels = c(1, 2, 3, 4, 5, 6,
                                         7, 8, 11, 12, 13,
                                         14, 15, 16, 17, 18
                                         ),
                  labels = c('urban', 'rural', 'rural', 'rural', 'urban', 'rural',
                             'rural', 'rural', 'urban', 'urban', 'rural',
                             'rural', 'rural', 'rural', 'rural', 'rural'
                             ),
                  exclude = 9
)

# Recode recruitment centres
df[, 'rec_centre'] <- factor(df[, 'rec_centre'],
                             levels = c(11012, 11021, 11011, 11008,
                                        11003, 11020, 11005, 11004,
                                        11018, 11010, 11016, 11001,
                                        11017, 11009, 11013, 11002,
                                        11007, 11014, 10003, 11006,
                                        11022, 11023
                                        ),
                  labels = c('Barts', 'Birmingham', 'Bristol', 'Bury',
                             'Cardiff', 'Croydon', 'Edinburgh', 'Glasgow',
                             'Hounslow', 'Leeds', 'Liverpool', 'Manchester',
                             'Middlesborough', 'Newcastel', 'Nottingham', 'Oxford',
                             'Reading', 'Sheffield', 'Stockport', 'Stoke',
                             'Swansea', 'Wrexham'
                  )
)

# Scale noise pollution by interquartile range
df[, paste0('noise_poll', '_s')] <- rescale(df[, 'noise_poll'])

# Categorise alcohol consumption as low,  medium, high; from Wang eClM 2024
df[, 'alcohol'] <- factor(df[, 'alcohol'],
                          levels = c(5, 6,  # special occasions, never
                                     3, 4,  # once/twice a week, one/three times a month
                                     1, 2   # daily/almost daily, three/four times a week
                                     ),
                  labels = c('low', 'low',
                             'medium', 'medium',
                             'high', 'high'
                             ),
                  exclude = -3,             # prefer not to answer
                  ordered = T
)

# Take average time out
df[, 't_out_winter'] <- case_match(df[, 't_out_winter'],
                                   -10 ~ 0,         # < 1 hour a day
                                   c(-1, -3)  ~ NA, # do not know, prefer not to say
                                   .default = df[, 't_out_winter']
)
df[, 't_out_summer'] <- case_match(df[, 't_out_summer'],
                                   -10 ~ 0,         # < 1 hour a day
                                   c(-1, -3)  ~ NA, # do not know, prefer not to say
                                   .default = df[, 't_out_summer']
)
df['t_out'] <- df %>% select(c(t_out_winter, t_out_summer)) %>% rowMeans(na.rm = T)

# Categorise physical activity based on total MET (metabolic equivalent of task):
# low < 600; moderate 600 >= 3000; intense >3000 MET-min/week (IPAQ 2005).
df['activity_score'] <- df %>%
    select(met_walking, met_mod_activity, met_vig_activity) %>%
    mutate(activity_score = rowSums(., na.rm = F),
           .keep = 'none'
    ) 
as_max <- df['activity_score'] %>% max(na.rm = T)
df <- df %>% mutate(activity_score = cut(activity_score, 
                                         breaks = c(0, 600, 3000, as_max), 
                                         include.lowest = T, 
                                         labels = c('low', 'moderate', 'intense'),
                                         ordered_result = F
    )
)


# Encode processed meat intake as ordinal variable
df[, 'processed_meat'] <- factor(df[, 'processed_meat'],
                                 levels = c(0, # never
                                            1, # less than once a week
                                            2, # once a week
                                            3, # 2-4 times a week
                                            4, # 5-6 times a week
                                            5  # one or more times daily
                                            ),
                  labels = c('never',
                             '<1week',
                             '1week',
                             '2-4week',
                             '>5week',
                             '>5week'
#                             '2-4week',
#                             '5-6week',
#                             '>=1day'
                             ),
                  exclude = c(-1, -3), # do not know, prefer not to answer
                  ordered = T
)

# Encode oily fish intake as ordinal variable; order from more often to less
# often intake, to match (assumed) effect of multiple deprivation index on incidence
# of dementia.
df[, 'fish'] <- factor(df[, 'fish'],
                                 levels = c(0, # never
                                            1, # less than once a week
                                            2, # once a week
                                            3, # 2-4 times a week
                                            4, # 5-6 times a week
                                            5  # one or more times daily
                                            ),
                  labels = c('never',
                             '<1week',
                             '1week',
                             '>2week',
                             '>2week',
                             '>2week'
#                             '2-4week',
#                             '5-6week',
#                             '>=1day'
                             ),
                  exclude = c(-1, -3), # do not know, prefer not to answer
                  ordered = T
) %>% fct_rev

# Clean sleep duration variable
df[, 'sleep_time'] <- case_match(df[, 'sleep_time'],
                                 c(-1, -3)  ~ NA, # do not know, prefer not to say
                                 .default = df[, 't_out_summer']
)

# Clean date sleep disorder
df[, 'date_sleep_disorder'] <- clean_dates(df[, 'date_sleep_disorder'])

# Categorise obesety as bmi >= 30
df <- within(df, obese <- bmi >= 30)
df[, 'obese'] <- factor(df[, 'obese'],
                        levels = c(F, T),
                        labels = c('non-obese', 'obese')
)

# Combine diabetes columns; retain earliest occurance
df <- df %>% mutate(date_diab = pmin(date_diab_ins,
                                     date_diab_noins,
                                     na.rm = T
                                     )
)

# Clean date diabetes
df[, 'date_diab'] <- clean_dates(df[, 'date_diab'])

# Clean date hypertension
df[, 'date_hypet'] <- clean_dates(df[, 'date_hypet'])

# Clean date atherosc
df[, 'date_atherosc'] <- clean_dates(df[, 'date_atherosc'])

# Categorise high cholesterol as >= 6.5 (Anstey et al., JAD, 2019)
df[, 'chol'] <- df[, 'chol'] >= 6.5
df[, 'chol'] <- factor(df[, 'chol'],
                       levels = c(F, T),
                       labels = c('non-high', 'high'),
                       ordered = F
)

# Clean date stroke
df[, 'date_stroke'] <- clean_dates(df[, 'date_stroke'])

# Clean date angina
df[, 'date_angina'] <- clean_dates(df[, 'date_angina'])

# Clean date infarct
df[, 'date_infarct'] <- clean_dates(df[, 'date_infarct'])


# Clean date asthma
df[, 'date_asthma'] <- clean_dates(df[, 'date_asthma'])

# Clear date copd
df[, 'date_copd'] <- clean_dates(df[, 'date_copd'])
 
# Combine mood disorders
df <- df %>% mutate(date_mood_disorder = pmin(date_depress,
                                              date_bipolar,
                                              na.rm = T
                                              )
)

# Clear date mood disorder
df[, 'date_mood_disorder'] <- clean_dates(df[, 'date_mood_disorder'])

# Clear date hear loss
df[, 'date_hearloss'] <- clean_dates(df[, 'date_hearloss'])

# Create cognitive score:
# log transform reaction times;
# +1 %>% log transform n incorrect matches with 6 pairs;
# Recode prospective memory;
# PCA with log reaction_time, +1 log pairs_match2, prosp_memory and fluid_intellgence;
# from Lyall et at., PlosOne, 2016
df[, 'reaction_time'] <- log(df[, 'reaction_time'])
df[, 'pairs_match2'] <- log(df[, 'pairs_match2'] + 1)
df[, 'prosp_memory'] <- case_match(df[, 'prosp_memory'],
                                   0 ~ -1, # Instruction not recalled, either skipped or incorrect
                                   2 ~ 0,  # correct recall on second attempt
                                   1 ~ 1,   # correct_recall on first attempt
                                   NA ~ NA
)
cog_measures <- c('reaction_time', 'pairs_match2', 'prosp_memory', 'fluid_intelligence')
pc_cog <- df %>%
    select(all_of(cog_measures)) %>%
    drop_na %>%
    prcomp(center = T, scale. = T)

df['cog_score'] <- df %>%
    select(all_of(cog_measures)) %>%
    as.matrix(.) %*% pc_cog[['rotation']][, 'PC1']  %>% .[, 1]
    

# Recode smoking status
df[, 'smoking'] <- factor(df[, 'smoking'],
                          levels = c(0, 1, 2),
                          labels = c('never', 'former', 'current'),
                          exclude = -3 # prefer not to answer
)

# Clean time at baseline address (address at assessment period)
df[, 'time_cuaddress'] <- case_match(df[, 'time_cuaddress'],
                                     c(-1, -3, -10) ~ NA, # do not know, prefer not to say, < 1 year
                                     .default = df[, 'time_cuaddress']
)

# Clean birth location east
df[, 'biloc_east'] <- case_match(df[, 'biloc_east'],
                                 -1 ~ NA, # location could not be mapped
                                 .default = df[, 'biloc_east']
)

# Clean birth location north
df[, 'biloc_north'] <- case_match(df[, 'biloc_north'],
                                  -1 ~ NA, # location could not be mapped
                                  .default = df[, 'biloc_north']
)

#df['time_cuaddress'] <- as.data.frame(lapply(df['time_cuaddress'], function(x){
#                                  gsub('Less than a year|Do not know|Prefer not to answer|^$', '0', x) %>%
#                                      as.integer()
#                                  }))

# Convert date_baseline to Date
df[, 'date_baseline'] <- as.Date(df[, 'date_baseline'])

# Add end of follow-up per country (based on rec_centre)
df[, 'end_study'] <- factor(df[, 'rec_centre'],
                            levels = c('Barts', 'Birmingham', 'Bristol', 'Bury',
                                     'Cardiff', 'Croydon', 'Edinburgh', 'Glasgow',
                                     'Hounslow', 'Leeds', 'Liverpool', 'Manchester',
                                     'Middlesborough', 'Newcastel', 'Nottingham', 'Oxford',
                                     'Reading', 'Sheffield', 'Stockport', 'Stoke',
                                     'Swansea', 'Wrexham'),
                  labels = c('2022-10-31', '2022-10-31', '2022-10-31', '2022-10-31',
                             '2022-05-31', '2022-10-31', '2022-08-31', '2022-08-31',
                             '2022-10-31', '2022-10-31', '2022-10-31', '2022-10-31',
                             '2022-10-31', '2022-10-31', '2022-10-31', '2022-10-31',
                             '2022-10-31', '2022-10-31', '2022-10-31', '2022-10-31',
                             '2022-05-31','2022-05-31')
) %>% ymd
 
# Save data -------------------------------------------------------------------
saveRDS(df, file = paste0(here(), '/data/processed/df_ukbb.rds'))
# saveRDS(df, file = 'df_ukbb.rda')
# write.csv(df, paste0(here(), '/data/processed/df_ukbb.csv'), row.names=F)

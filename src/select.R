#' Select data
#'
#' Select ukbb data; drop participants who lived less than address_buff years at baseline address (time_cuaddress);
#' drop participants with all cause dementia incidence before time baseline + dem_buff; drop individual with
#' missing covariates.
#' Additionally: if copd = T (positive control analysis), ignore time of dementia incidence and drop participants 
#' based on copd time (not dem) compared to baseline + dem_buf (default 0); if fish = T (negative control 
#' analysis), ignore time of dementia incidence and drop participants based on presence of response to oily fish 
#' consumption frequency at baseline, plus keep only people with time_cuaddress within address_buff boundaries.
#' Attention: dropping participants that did not know or did not want to say how many years lived at baseline.
#'
#' @param df_path: path of df_ukbb.csv
#' @param address_buff: time at baseline address
#' @param dem_buff: buffer time between baseline assessment and dementia incidence
#' @param copd: if true, keep based on copd times, not dem
#' @param fish: if true, drop if na and drop participants if lived more than
#' address_buff (attention: dem/copd its less than address_buff)
#' @param ...: covariates
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#' @date 2024-03-27

select_data <- function(df_all, 
                        ...,
                        address_buff = 1, 
                        dem_buff = 0, 
                        copd = F,
                        fish = F){
    # Initialise
    dem_buff <- dem_buff * 365 # convert years in days
    list_cov <- list(...)
    `%>%` <- dplyr::`%>%`

    # Save initial number of participants
    n_all <- df_all %>% dim %>% .[1]

    # Drop participants without pm data
    # attention: no2/x has ~30k more data; here dropped for consistency with other apoll!
    mask_pm <- df_all['pm25'] %>% complete.cases
    n_apoll_keep <- sum(mask_pm)
    n_apoll_drop <- sum(!mask_pm)
    eid_drop <- df_all[!mask_pm, 'eid'] # save electronic id of dropped
    df <- df_all[mask_pm, ]

    # Drop based on dementia incidence time (don't if copd/fish=T)
    if (copd == F && fish == F) {
        mask_db <- df['date_dem'] >= (df['date_baseline'] + dem_buff)
    } else { # don't drop for copd analysis
        mask_db <- df[, 'date_dem'] >= ymd('1900-01-01')
    }
    n_dem <- sum(mask_db, na.rm = T) # save n dem cases
    n_dem_drop <- sum(!mask_db, na.rm = T) # save n dropped
    n_dem_keep <- n_apoll_keep - n_dem_drop # save n keep
    mask_db[is.na(mask_db)] <- T # include participants without dementia
    df <- df[mask_db, ]

    # Drop based on copd incidence time
    n_copd = NA
    n_copd_drop = 0
    n_copd_keep = n_dem_keep 
    if (copd == T) {
        mask_cb <- df['date_copd'] >= (df['date_baseline'])
        n_copd <- sum(mask_cb, na.rm = T) # save n copd cases
        n_copd_drop <- sum(!mask_cb, na.rm = T) # save n dropped
        n_copd_keep <- n_copd_keep - n_copd_drop # save n dropped
        mask_cb[is.na(mask_cb)] <- T # include participants without copd
        df <- df[mask_cb, ]
    }

    # Drop based on oily fish response at baseline
    n_fish = NA
    n_fish_drop = 0
    n_fish_keep = n_dem_keep
    if (fish == T) {
        mask_fish <- complete.cases(df['fish'])
        n_fish <- sum(mask_fish, na.rm = T) # save n fish cases
        n_fish_drop <- sum(!mask_fish, na.rm = T) # save n dropped
        n_fish_keep <- n_fish_keep  - n_fish_drop # save n dropped
        df <- df[mask_fish, ]
    }

    # Drop based on covariate missing values
    n_cov_keep_all = df %>% dim %>% .[1]
    n_cov_drop_all = 0
    n_cov_drop = list()
    if (!!length(list_cov)) { # if list_cov has length > 0
        for (cov_ in list_cov[[1]]) {
            mask_cov = complete.cases(df[cov_])
            n_cov_drop_all = n_cov_drop_all + sum(!mask_cov)
            n_cov_drop[cov_] = sum(!mask_cov)
            if (cov_ == 'mdi_eng') {
                n_mdi_scot <- df[['mdi_scot']] %>% complete.cases %>% sum
                n_mdi_wales <- df[['mdi_wales']] %>% complete.cases %>% sum
                x100mdi_eng <- 100 - ((n_mdi_scot + n_mdi_wales) / 
                                      sum(mask_cov) * 100)
                
            }
            df = df[mask_cov, ]
        }
    }
    n_cov_keep_all = n_cov_keep_all - n_cov_drop_all  

    # Drop based on time at address (default 1 year); if fish = T,
    # then keep people within address_buff boundaries instead
    if (fish == T) {
        mask_ab <- df['time_cuaddress'] >= address_buff[1] & df['time_cuaddress'] <= address_buff[2]
    } else {
        mask_ab <- df['time_cuaddress'] >= address_buff
    }
    mask_ab[is.na(mask_ab)] <- F # exclude codes -10, -1, -3
    n_address_keep <- sum(mask_ab) # save n remaining
    n_address_drop <- sum(!mask_ab, na.rm = T) # save n dropped
    df <- df[mask_ab, ]

    # Print n participats kept and dropped
    cat('apoll keep', n_apoll_keep, '\n', 'apoll drop', n_apoll_drop,
        '\n', 'n dem:', n_dem, '\n', 'dem drop:', n_dem_drop,
        '\n', 'n copd:', n_copd, '\n', 'copd drop:', n_copd_drop, 
        '\n', 'n fish :', n_fish, '\n', 'fish drop:', n_fish_drop, 
        '\n', 'cov keep:', n_cov_keep_all, 'cov drop:', n_cov_drop_all,
        '\n', 'address keep', n_address_keep, '\n', 'address drop', n_address_drop)
    map(n_cov_drop, 1) %>% print

    # Return output
    return(list(df = df, 
                eid_drop = eid_drop,
                n_all = n_all,
                apoll_keep = n_apoll_keep,
                apoll_drop = n_apoll_drop,
                n_dem = n_dem,
                dem_keep = n_dem_keep,
                dem_drop = n_dem_drop,
                n_copd = n_copd,
                copd_keep = n_copd_keep,
                copd_drop = n_copd_drop,
                n_fish = n_fish,
                fish_keep = n_fish_keep,
                fish_drop = n_fish_drop,
                cov_keep_all = n_cov_keep_all,
                cov_drop_all = n_cov_drop_all,
                cov_drop = n_cov_drop,
                x100mdi_eng = x100mdi_eng,
                address_keep = n_address_keep,
                address_drop = n_address_drop
                )
    )
}

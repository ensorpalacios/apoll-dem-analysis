#' Select data
#'
#' Select ukbb data; drop participants who lived less that address_buff years at baseline address (time_cuaddress);
#' drop participants with all cause dementia incidence before time baseline + dem_buff; drop individual with
#' missing covariates.
#' Attention: dropping participants that did not know or did not want to say how many years lived at baseline.
#'
#' @param df_path: path of df_ukbb.csv
#' @param address_buff: time at baseline address
#' @param dem_buff: buffer time between baseline assessment and dementia incidence
#' @param ...: covariates
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#' @date 2024-03-27

select_data <- function(df_all, 
                        address_buff = 1, 
                        dem_buff = 0, 
                        drop_nowhite = T,
                        ...){
    # Initialise
    dem_buff <- dem_buff * 365 # convert years in days
    list_cov <- list(...)
    `%>%` <- dplyr::`%>%`

    # Drop participants without pm data
    # attention: no2/x has ~30k more data; here dropped!
    mask_pm <- df_all['pm25'] %>% complete.cases
    n_apoll_drop <- sum(!mask_pm)
    eid_drop <- df_all[!mask_pm, 'eid']
    df <- df_all[mask_pm, ]

    # Drop based on time at address (default 1 year)
    mask_ab <- df['time_cuaddress'] >= address_buff
    mask_ab[is.na(mask_ab)] <- F # exclude codes -10, -1, -3
    n_address_keep <- sum(mask_ab) # save n remaining
    n_address_drop <- sum(!mask_ab, na.rm = T) # save n dropped
    df <- df[mask_ab, ]

    # Drop based on dementia incidence time
    mask_db <- df['date_dem'] >= (df['date_baseline'] + dem_buff)
    n_dem <- sum(mask_db, na.rm = T) # save n dem cases
    n_dem_drop <- sum(! mask_db, na.rm = T) # save n dropped
    mask_db[is.na(mask_db)] <- T # include participants without dementia
    df <- df[mask_db, ]

    # Drop based on ethnicity
    if (drop_nowhite) {
        mask_w <- df['ethnic'] == 'white'
        mask_w[is.na(mask_w)] <- F # drop NA
        n_ethnic_keep <- sum(mask_w)
        n_ethnic_drop <- sum(!mask_w)
        df <- df[mask_w, ]
    } else {
        n_ethnic_keep <- 'all'
        n_ethnic_drop <- 'none'
    }

    # Drop based on covariate missing values
    n_cov_drop = 0
    for (cov_ in list_cov) {
        mask_cov = complete.cases(df[cov_])
        n_cov_drop = n_cov_drop + sum(!mask_cov)
        df = df[complete.cases(df[cov_]), ]
    }

    # Print n participats kept and dropped
    cat('apoll drop', n_apoll_drop,
        '\n', 'address keep', n_address_keep, '\n', 'address drop', n_address_drop,
        '\n', 'n dem:', n_dem, '\n', 'dem drop:', n_dem_drop,
        '\n', 'ethnic keep', n_ethnic_keep, '\n', 'ethnic drop', n_ethnic_drop,
        '\n', 'cov drop:', n_cov_drop)

    return(list(df = df, 
                eid_drop = eid_drop
                )
    )
}

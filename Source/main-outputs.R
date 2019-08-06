source("haiti_ocv_JHUmodel/Source/utils_clean.R")
source("haiti_ocv_JHUmodel/Source/utils_eval.R")
reload_source()

source("haiti_ocv_JHUmodel/Source/src_generate_fc_params.R")
out_dir <- paste0("haiti_ocv_JHUmodel/GeneratedData/", "final1_summaries/")
dir.create(out_dir, showWarnings = FALSE)


## -------------------------------------- ##
## script settings ##
ts_module <- TRUE; pr_elim_module <- TRUE; t_elim_module <- TRUE; numvacc_module <- TRUE


## -------------------------------------- ##
fix_fc_settings_otf <- function(modcode, scencode){ ## on-the-fly
  fc_seed <- 20192104
  fc_nsims <- 25
  fc_horizon <- 52*11 ## 11 years forecast

  fc_nd = NA; fc_nw = NA; fc_c2 = NA; fc_c1 = NA; fc_vescen = NA

  ## grouped by deployment strategy
  if (scencode %in% paste0("id", seq(2, 34, by = 4))){ ## 2 dept campaigns
    fc_nd = 2; fc_nw = 52
  } else if (scencode %in% paste0("id", seq(4, 36, by = 4))){ ## 3 dept campaigns
    fc_nd = 3; fc_nw = 33
  } else if (scencode %in% paste0("id", seq(3, 35, by = 4))){ ## slow national campaigns
    fc_nd = 10; fc_nw = 26
  } else if (scencode %in% paste0("id", seq(1, 33, by = 4))){ ## fast national campaigns
    fc_nd = 10; fc_nw = 10
  } 

  ## grouped by coverage level
  if (scencode %in% paste0("id", 1:12)){
    fc_c2 = 0.7; fc_c1 = 0.1
  } else if (scencode %in% paste0("id", 13:24)){
    fc_c2 = 0.4; fc_c1 = 0.2
  } else if (scencode %in% paste0("id", 25:36)){
    fc_c2 = 0.95; fc_c1 = 0.0167
  }

  ## grouped by ve scenario
  if (scencode %in% paste0("id", c(1:4, 13:16, 25:28))){
    fc_vescen = "ve_s1"
  } else if (scencode %in% paste0("id", c(5:8, 17:20, 29:32))){
    fc_vescen = "ve_s2"
  } else if (scencode %in% paste0("id", c(9:12, 21:24, 33:36))){
    fc_vescen = "ve_s3"
  }

  fc_str_if <- paste0(scencode, "_nsims", fc_nsims, "_seed", fc_seed)
  fc_str_if2 <- paste0(scencode, "_nsims", fc_nsims)

  return(list(mcode=modcode, scode=scencode, seed=fc_seed, horizon = fc_horizon, nsims=fc_nsims, nd = fc_nd, nw = fc_nw, c2 = fc_c2, c1 = fc_c1, vescen = fc_vescen, str_if = fc_str_if, str_if2 = fc_str_if2))
}

## -------------------------------------- ##
## get ts data for vacmodels ##
if(ts_module){

  scodes <- c("id1", "id2", "id3", "id4", "id25", "novac")
  ts_outputs <- map_dfr(scodes, function(sc){
    mcode <- "final1"; scode = sc

    episettings <- fix_epi_settings()
    endsettings <- fix_end_settings()
    fcsettings <- fix_fc_settings_otf(mcode, scode)
    out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
    out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
    out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/")

    ## -------------------------------------- ##
    ## data import and setup ##
    pop.haiti <- get.haiti.pop()

    haiti.dat <- get.mspp.agg.data() 
    haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
    haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
      dplyr::mutate(week_end = seq_along(week)-1)

    covartab.all <- make.covartab(0, nrow(haiti.dat)+1, byt=1, degree=6, nbasis=6, per=52.14)
    covartab <- covartab.all[(episettings$nweeks+1):nrow(covartab.all),] %>%
      dplyr::mutate(time_end = seq_along(time)-1)
    num_betas <- ncol(covartab.all)-1


    fn_epi <- list.files(out_dir_epi, "fitstates_if_")
    fn_epi_full <- paste0(out_dir_epi, fn_epi)
    fn_end <- list.files(out_dir_end, "fitstates_if_")
    fn_end_full <- paste0(out_dir_end, fn_end)
    fn_fc_elim <- list.files(out_dir_fc, "fcstates_stoch_")
    fn_fc_elim_full <- paste0(out_dir_fc, fn_fc_elim)
    fn_fc_plt <- list.files(out_dir_fc, "stochfcs_") ## summary across all paramsets & sims
    fn_fc_plt_full <- paste0(out_dir_fc, fn_fc_plt)

    ## -------------------------------------- ##
    ## ts data processing for a single model ##
    fit_epidat <- process_fepi_plot(fn_epi_full, haiti.dat.epi) %>%
      dplyr::mutate(period = "epidemic")
    fit_enddat <- process_fend_plot(fn_end_full, haiti.dat.end) %>% 
      dplyr::filter(week > episettings$nweeks) %>%
      dplyr::mutate(period = "endemic")
    fitdat <- bind_rows(fit_epidat, fit_enddat)

    fcdat <- process_ffc_plot(fn_fc_plt_full) %>%
      dplyr::filter(week > max(haiti.dat$week)) %>%
      dplyr::mutate(period = "forecast")
    ffcdat <- bind_rows(fitdat, fcdat)

    date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
    datedat <- data.frame(weekdate = seq(date.start$date_sat, length.out = nrow(ffcdat), by = 7), week = seq_len(nrow(ffcdat)))

    plotdat_ffc <- full_join(ffcdat, datedat, by = c("week")) %>%
      dplyr::mutate(scenario = scode)

    return(plotdat_ffc)

  }) %>% 
    dplyr::mutate(scenario = 
                    case_when(scenario == "id1" ~ "1",
                              scenario == "id2" ~ "2",
                              scenario == "id3" ~ "3",
                              scenario == "id4" ~ "4",
                              scenario == "id25" ~ "25",
                              scenario == "novac" ~ "0")) %>%
    dplyr::mutate(team = "JHU") %>%
    dplyr::rename(obs_cases_med = cases_med,
                  obs_cases_low = cases_lo,
                  obs_cases_hi = cases_hi,
                  true_inf_med = incid_med,
                  true_inf_low = incid_lo,
                  true_inf_hi = incid_hi,
                  population_size = N_med)
    
  vacts_output <- process_cumvac_output() %>%
    dplyr::select(weekdate, scenario, vacc_med) 

  ## ts_[team name]_[output date in the format yyyy-mm-dd].csv
  ## team scenario  saturday_date obs_cases_med obs_cases_low obs_cases_hi  true_inf_med  true_inf_low  true_inf_hi population_size vacc_med  vacc_low  vacc_hi
  ts_output_full <- full_join(ts_outputs, vacts_output, by = c("weekdate", "scenario")) %>%
    dplyr::rename(saturday_date = weekdate) %>%
    dplyr::select(team, scenario, saturday_date, obs_cases_med, obs_cases_low, obs_cases_hi, true_inf_med, true_inf_low, true_inf_hi, population_size, vacc_med) %>%
    dplyr::mutate(vacc_low = NA, vacc_hi = NA)

  write_csv(ts_output_full, paste0(out_dir, "ts_JHU_", Sys.Date(), ".csv"))

}


# ## -------------------------------------- #
## save elimination ##
## pr_elim_[team name]_[output date in the format yyyy-mm-dd].csv
## team scenario  time_frame  pr_elim_thresh1 pr_resurg_thresh1 pr_elim_thresh2 pr_resurg_thresh2 cum_inf_med cum_inf_low cum_inf_hi

if(pr_elim_module){

  scodes <- c("id1", "id2", "id3", "id4", "id25", "novac")
  prelim_outputs <- map_dfr(scodes, function(sc){
    mcode <- "final1"; scode = sc

    episettings <- fix_epi_settings()
    endsettings <- fix_end_settings()
    fcsettings <- fix_fc_settings_otf(mcode, scode)
    out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", mcode, "_", scode, "/")
    fcelim <- readRDS(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds")) %>%
      dplyr::mutate(scenario = scode)

    pElim_dat <- find_probElim(fcelim) %>% 
      dplyr::mutate(scenario = scode) %>%
      dplyr::select(scenario, time_frame, pr_elim_thresh1, pr_resurg_thresh1, pr_elim_thresh2, pr_resurg_thresh2)

    return(pElim_dat)
  }) %>% 
    dplyr::mutate(scenario = 
                    case_when(scenario == "id1" ~ 1,
                              scenario == "id2" ~ 2,
                              scenario == "id3" ~ 3,
                              scenario == "id4" ~ 4,
                              scenario == "id25" ~ 25,
                              scenario == "novac" ~ 0))

  ## identify relevant time points
  vacstart_wk <- max(get.mspp.agg.data()$week)+4
  year3_wk <- round(vacstart_wk + 3*52.14)
  year5_wk <- round(vacstart_wk + 5*52.14)
  year10_wk <- round(vacstart_wk + 10*52.14)

  vacstart_dt <- get_date_df()[which(get_date_df()$week == vacstart_wk),]$saturday_date
  year3_dt <- get_date_df()[which(get_date_df()$week == year3_wk),]$saturday_date
  year5_dt <- get_date_df()[which(get_date_df()$week == year5_wk),]$saturday_date
  year10_dt <- get_date_df()[which(get_date_df()$week == year10_wk),]$saturday_date

  ts_output_full <- read_csv(paste0(out_dir, "ts_JHU_", Sys.Date(), ".csv")) %>%
    dplyr::filter(saturday_date >= vacstart_dt) %>%
    group_by(scenario) %>%
    dplyr::mutate(cum_inf_med = cumsum(true_inf_med), cum_inf_low = cumsum(true_inf_low), cum_inf_hi = cumsum(true_inf_hi)) %>%
    ungroup %>%
    dplyr::filter(saturday_date %in% c(year3_dt, year5_dt, year10_dt)) %>%
    dplyr::mutate(time_frame = case_when(saturday_date == as.Date(year3_dt) ~ 3,
                                        saturday_date == as.Date(year5_dt) ~ 5,
                                        saturday_date == as.Date(year10_dt) ~ 10)) %>%
    dplyr::select(team, scenario, time_frame, contains("cum_"))

  pr_elim_full <- full_join(prelim_outputs, ts_output_full, by = c("scenario", "time_frame")) %>%
    dplyr::select(team, scenario, time_frame, pr_elim_thresh1, pr_resurg_thresh1, pr_elim_thresh2, pr_resurg_thresh2, cum_inf_med, cum_inf_low, cum_inf_hi)

  write_csv(pr_elim_full, paste0(out_dir, "pr_elim_JHU_", Sys.Date(), ".csv"))

}


## -------------------------------------- ##
## save time to elimination ##
## tte_[team name]_[output date in the format yyyy-mm-dd].csv
## team scenario  t_elim_thresh1_med  t_elim_thresh1_low  t_elim_thresh1_hi t_elim_thresh2_med  t_elim_thresh2_low  t_elim_thresh2_hi vacc_startdate

if(t_elim_module){

  scodes <- c("id1", "id2", "id3", "id4", "id25", "novac")
  telim_outputs <- map_dfr(scodes, function(sc){
    mcode <- "final1"; scode = sc

    episettings <- fix_epi_settings()
    endsettings <- fix_end_settings()
    fcsettings <- fix_fc_settings_otf(mcode, scode)
    out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", mcode, "_", scode, "/")
    fcelim <- readRDS(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds")) %>%
      dplyr::mutate(scenario = scode)

    tElim_dat <- find_timeToElim(fcelim) %>% 
      dplyr::mutate(scenario = scode) %>%
      dplyr::select(scenario, t_elim_thresh1_med, t_elim_thresh1_low, t_elim_thresh1_hi, t_elim_thresh2_med,  t_elim_thresh2_low, t_elim_thresh2_hi)

    return(tElim_dat)
  }) %>% 
    dplyr::mutate(scenario = 
                    case_when(scenario == "id1" ~ 1,
                              scenario == "id2" ~ 2,
                              scenario == "id3" ~ 3,
                              scenario == "id4" ~ 4,
                              scenario == "id25" ~ 25,
                              scenario == "novac" ~ 0)) %>%
    dplyr::mutate(team = "JHU", vacc_startdate = vacstart_wk <- max(get.mspp.agg.data()$week)+4) %>%
    dplyr::select(team, scenario, t_elim_thresh1_med, t_elim_thresh1_low, t_elim_thresh1_hi, t_elim_thresh2_med,  t_elim_thresh2_low, t_elim_thresh2_hi)

  write_csv(telim_outputs, paste0(out_dir, "tte_JHU_", Sys.Date(), ".csv"))

}

if(numvacc_module){
  ## approximately how many doses of vaccine are needed for each campaign?

  scencodes <- paste0("id", c(1,2,3,4,25))
  cov_df <- data.frame(scode = scencodes, cov2d = c(rep(.7, 4), .95), cov1d = c(rep(.1, 4), .0167))

  pop_df <- map_dfr(scencodes, function(sc){
    dummy <- data.frame(scode = sc, ocv_order = 1:10, dept = c("Centre", "Artibonite", "Ouest", "Nord Ouest", "Nord", "Sud", "Nippes", "Nord Est", "Sud Est", "Grand'Anse"), pop = c(746236, 1727524, 4029705, 728807, 1067177, 774976, 342525, 393967, 632601, 468301))
    if(sc == "id2"){
      df <- dummy %>% dplyr::filter(ocv_order <= 2)
    } else if(sc == "id4"){
      df <- dummy %>% dplyr::filter(ocv_order <= 3)
    } else{
      df <- dummy
    }
    return(df)
  })
  
  full_df <- pop_df %>%
    group_by(scode) %>%
    summarise(targetpop = sum(pop)) %>%
    left_join(cov_df, by = c("scode")) %>%
    dplyr::mutate(vaccines_needed = cov2d*2*targetpop + cov1d*1*targetpop)
  

}

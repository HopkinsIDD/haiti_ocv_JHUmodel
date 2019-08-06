
#' @description model projection settings 
fix_fc_settings <- function(){
  modcode <- "final"
  scencode <- "id1" ## should be a string "id#", where # represents the code of a different vaccination deployment strategy (1 = fast national, 2 = 2-department, 3 = slow national, 4 = 3-department, 25 = fast national high coverage, and all other scenarios represent combinations of each)
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


#' @description model fit settings for epidemic model that is feeding into the projection model
fix_epi_settings <- function(){
  epi_nweeks <- 232 ## 180 end of 3/2014; 219 end of 2014; 232 end of 3/2015
  epi_seed <- 20190415
  epi_nsamps <- NA; epi_nprofs <- 3
  epi_nmif <- 50
  epi_nparticles <- 100
  epi_passTm <- "true"
  modcode <- "final"
  usepriors <- FALSE

  ## create strings for filenames
  epi_str_tm <- paste0("epi", epi_nweeks, ifelse(is.na(epi_nsamps), paste0("_nprofs", epi_nprofs), paste0("_nsamps", epi_nsamps)), "_seed", epi_seed)
  epi_str_tm2 <- paste0("epi", epi_nweeks, ifelse(is.na(epi_nsamps), paste0("_nprofs", epi_nprofs), paste0("_nsamps", epi_nsamps)))
  epi_str_if <- paste0("epi", epi_nweeks, ifelse(is.na(epi_nsamps), paste0("_nprofs", epi_nprofs), paste0("_nsamps", epi_nsamps)), "_nmif", epi_nmif, "_nparticles", epi_nparticles, "_passTm", epi_passTm, "_seed", epi_seed)
  epi_str_if2 <- paste0("epi", epi_nweeks, ifelse(is.na(epi_nsamps), paste0("_nprofs", epi_nprofs), paste0("_nsamps", epi_nsamps)), "_nmif", epi_nmif, "_nparticles", epi_nparticles, "_passTm", epi_passTm)

  return(list(nweeks=epi_nweeks, seed=epi_seed, nsamps=epi_nsamps, nmif=epi_nmif, nparticles=epi_nparticles, mcode=modcode, passTm=epi_passTm, priors=usepriors, nprofs=epi_nprofs, str_tm=epi_str_tm, str_tm2=epi_str_tm2, str_if=epi_str_if, str_if2=epi_str_if2))
}


#' @description model fit settings for endemic model that is feeding into the projection model
fix_end_settings <- function(){

  episettings <- fix_epi_settings()

  end_seed <- 4172019
  end_nmif <- 50
  end_nparticles <- 100
  end_passTm <- "false"
  modcode <- paste0(episettings$mcode, "1")

  end_str_tm <- paste0("end", episettings$nweeks, "_seed", end_seed)
  end_str_tm2 <- paste0("end", episettings$nweeks)
  end_str_if <- paste0("end", episettings$nweeks, "_nmif", end_nmif, "_nparticles", end_nparticles, "_passTm", end_passTm, "_seed", end_seed)
  end_str_if2 <- paste0("end", episettings$nweeks, "_nmif", end_nmif, "_nparticles", end_nparticles, "_passTm", end_passTm)

  return(list(seed=end_seed, nmif=end_nmif, nparticles=end_nparticles, mcode=modcode, passTm=end_passTm, str_tm = end_str_tm, str_tm2 = end_str_tm2, str_if = end_str_if, str_if2 = end_str_if2))
}


#' @description set parameters for projection model
generate_fc_params <- function(out_dir_epi, out_dir_end, out_dir_fc){
  fcsettings <- fix_fc_settings()
  endsettings <- fix_end_settings()
  episettings <- fix_epi_settings()

  set.seed(fcsettings$seed)

  kap <- 0.95

  est_epi <- readRDS(paste0(out_dir_epi, "iffits_", episettings$str_if, ".rds")) %>%
    dplyr::select(parid, loglik) %>%
    dplyr::rename(loglik_epi = loglik)
  est_end <- readRDS(paste0(out_dir_end, "iffits_", endsettings$str_if, ".rds")) %>%
    dplyr::filter(loglik >-3000) %>%
    dplyr::mutate(mcode = endsettings$mcode) %>%
    dplyr::mutate(parid = as.integer(parid)) %>%
    left_join(est_epi, by = c("parid")) %>%
    dplyr::mutate(loglik_f = loglik + loglik_epi) %>%
    dplyr::select(mcode, parid, loglik_f, rho, theta, beta1, beta2, beta3, beta4, beta5, beta6, nu, gamma, sigma, theta0, alpha, mu, delta) %>%
    dplyr::mutate(kappa = kap)

  incl_parids <- unlist(est_end$parid)
  fitstate_fns <- paste0("fitstates_if_", endsettings$str_if2, "_parid", incl_parids, "_seed", endsettings$seed, ".csv")
  initconds <- map_dfr(1:length(fitstate_fns), function(i){
    par_id <- incl_parids[i]
    read_csv(paste0(out_dir_end, fitstate_fns[i])) %>%
      dplyr::filter(week == max(week)) %>%
      dplyr::mutate(S.0 = S_med/N_med, E.0 = E_med/N_med, I.0 = I_med/N_med, A.0 = A_med/N_med, incid.0 = incid_med) %>%
      dplyr::mutate(R.0 = 1-(S.0+E.0+I.0+A.0), parid = as.numeric(par_id), N0 = N_med) %>%
      dplyr::mutate(E.0 = ifelse(E.0 == 0, 1E-7, E.0), I.0 = ifelse(I.0 == 0, 1E-7, I.0)) %>%
      dplyr::select(parid, S.0, E.0, I.0, A.0, R.0, incid.0, N0)
  })

  starts_fc <- inner_join(est_end, initconds, by = c("parid")) %>%
    dplyr::filter(!is.na(S.0)) %>%
    dplyr::arrange(desc(loglik_f))
  write_csv(starts_fc, paste0(out_dir_fc, "starts_", fcsettings$mcode, "_epi", episettings$nweeks, ".csv"))

  return(starts_fc)
}
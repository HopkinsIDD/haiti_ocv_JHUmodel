
#' @description model fit settings for epidemic model feeding into the given endemic model
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


#' @description model fit settings for the given endemic model
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


#' @description set parameters for endemic fit model
generate_end_params <- function(out_dir_epi, out_dir_end){

  episettings <- fix_epi_settings()
  endsettings <- fix_end_settings()

  set.seed(endsettings$seed)
  est_epi <- readRDS(paste0(out_dir_epi, "/iffits_", episettings$str_if, ".rds")) 

  if(file.exists(paste0(out_dir_end, "starts_", endsettings$mcode, "_epi", episettings$nweeks, ".csv"))){
    starts_end <- read_csv(paste0(out_dir_end, "starts_", endsettings$mcode, "_epi", episettings$nweeks, ".csv"))

  } else{
    starts_end1 <- est_epi %>%
    dplyr::filter(nfail.max == 0) %>%
    dplyr::mutate(epi_mcode = episettings$mcode) %>%
    dplyr::mutate(parid = as.integer(parid)) %>%
    dplyr::select(epi_mcode, parid, loglik, rho, theta, beta1, beta2, beta3, beta4, beta5, beta6, nu, gamma, sigma, theta0, alpha, mu, delta)

    fitstate_fns <- list.files(out_dir_epi, paste0("fitstates_if_", episettings$str_if2))
    initconds <- map_dfr(1:length(fitstate_fns), function(i){
      par_id <- gsub(paste0("_seed", episettings$seed, ".csv"), "", gsub(paste0("fitstates_if_", episettings$str_if2, "_parid"), "", fitstate_fns[i]))
      read_csv(paste0(out_dir_epi, fitstate_fns[i])) %>%
        dplyr::filter(week == max(week)) %>%
        dplyr::mutate(S.0 = S_med/N_med, E.0 = E_med/N_med, I.0 = I_med/N_med, A.0 = A_med/N_med, incid.0 = incid_med) %>%
        dplyr::mutate(R.0 = 1-(S.0+E.0+I.0+A.0), parid = as.numeric(par_id), N0 = N_med) %>%
        dplyr::select(parid, S.0, E.0, I.0, A.0, R.0, incid.0, N0)
    }) %>%
      arrange(parid)

    starts_end <- inner_join(starts_end1, initconds, by = c("parid")) %>%
      dplyr::filter(!is.na(S.0)) %>%
      dplyr::filter(nu > 0.9 & beta1 < 100) %>% ## filter reasonable parameter ranges
      dplyr::arrange(desc(loglik))
    write_csv(starts_end, paste0(out_dir_end, "starts_", endsettings$mcode, "_epi", episettings$nweeks, ".csv")) 
  }

  return(starts_end)
}


#' @description explore coverage of starting endemic parameters
explore_end_inits_coverage <- function(fig_dir){

  epi_mcode <- fix_epi_settings()$mcode
  end_mcode <- fix_end_settings()$mcode
  nweeks <- fix_epi_settings()$nweeks
  out_dir_epi <- paste0("haiti_mass_ocv/GeneratedData/", epi_mcode, "/")
  out_dir_end <- paste0("haiti_mass_ocv/GeneratedData/", end_mcode, "/")
  starts_end <- generate_end_params(out_dir_epi, out_dir_end)

  pdf(file = paste0(fig_dir, "loglik_", fix_epi_settings()$str_if, ".pdf"), width=4, height=4)
  hist(starts_end$loglik)
  dev.off()

  pltpars <- starts_end  %>% 
    dplyr::select(rho, theta, beta1, nu, loglik)
  plt <- ggpairs(pltpars)
  ggsave(paste0(fig_dir, "start_coverage_", end_mcode, "_epi", nweeks, ".pdf"), plt, width = 6, height = 6, units = "in")
  return(plt)
}

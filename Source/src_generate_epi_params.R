
#' @description model fit settings for epidemic fit model
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


#' @description set parameters for epidemic fit model
generate_epi_params <- function(out_dir, priors){

  parseed <- fix_epi_settings()$seed
  nweeks <- fix_epi_settings()$nweeks
  nsamps <- fix_epi_settings()$nsamps
  mcode <- fix_epi_settings()$mcode
  nprofs <- fix_epi_settings()$nprofs

  if(file.exists(paste0(out_dir, "starts_", mcode, "_epi", nweeks, ifelse(is.na(nsamps), paste0("_nprofs", nprofs), paste0("_nsamps", nsamps)), "_", parseed, ".csv"))){
    starts_epi <- read_csv(paste0(out_dir, "starts_", mcode, "_epi", nweeks, ifelse(is.na(nsamps), paste0("_nprofs", nprofs), paste0("_nsamps", nsamps)), "_", parseed, ".csv"))
  
  } else{
    set.seed(parseed)

    pop <- get.haiti.pop()

    ## Fixed parameters
    gamma <- 7/2 ## 2 day infectious period
    sigma <- 7/1.4 ## 1.4 day latent period
    theta0 <- 0 ## 0% going to asymptomatic in base model
    alpha <- 7/2920 ## 8 year mean duration of natural immunity
    mu <- ((1+22.6/1000)^(1/52.14))-1 ## 22.6/1000 average annual birth rate, adjusted for compounding by week
    delta <- ((1+7.5/1000)^(1/52.14))-1 ## 7.5/1000 average annual death rate, adjusted for compounding by week

    ## Starting states
    E0 <- 10/pop ## rpois(nsamps, 10)/pop
    I0 <- 10/pop ## rpois(nsamps, 10)/pop
    A0 <- 0.0/pop
    R0 <- 0.000
    S0 <- 1-R0-I0-E0-A0

    ## beta parameter settings
    blo <- 1E-9; bup <- 10 ## uniform beta settings
    bmn <- 4.5; bse <- 0.5 ## median R0 among initial values: median(exp(rnorm(1000, log(bmn), bse)))*2/7
    ## nu parameter settings
    nlo <- 0.95; nup <- 1

    ## Fitted parameters
    if(priors){
      ## create data frame of initial starting parameters for epidemic fitting procedure
      theta <- runif(nsamps, 1, 10)
      rho_t <- rnorm(nsamps,log(.25/.75), 1.5)
      rho <- exp(rho_t)/(1+exp(rho_t)) ## effective range of rho is roughly 0.03 to 0.85 with central tendency at 0.25
      
      ## uniform dist of betas
      beta1 <- runif(nsamps, blo, bup)
      beta2 <- runif(nsamps, blo, bup)
      beta3 <- runif(nsamps, blo, bup)
      beta4 <- runif(nsamps, blo, bup)
      beta5 <- runif(nsamps, blo, bup)
      beta6 <- runif(nsamps, blo, bup)

      ## uniform dist for nu
      nu <- runif(nsamps, nlo, nup)

      starts_epi <- data.frame(parid=1:nsamps, rho=rho, theta=theta, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4, beta5=beta5, beta6=beta6, nu=nu, gamma=gamma, sigma=sigma, theta0=theta0, alpha=alpha, mu=mu, delta=delta, S.0=S0, E.0=E0, I.0=I0, A.0=A0, R.0=R0)

    } else{
      rhosamps <- profileDesign(rho = seq(1E-8, 1, length = 30),
                    upper = c(theta=20, beta1=bup, beta2=bup, beta3=bup, beta4=bup, beta5=bup, beta6=bup, nu=nup),
                    lower = c(theta=1, beta1=blo, beta2=blo, beta3=blo, beta4=blo, beta5=blo, beta6=blo, nu=nlo),
                    nprof = nprofs)
      thetasamps <- profileDesign(theta = seq(1, 20, length = 30),
                    upper = c(rho=1, beta1=bup, beta2=bup, beta3=bup, beta4=bup, beta5=bup, beta6=bup, nu=nup),
                    lower = c(rho=1E-8, beta1=blo, beta2=blo, beta3=blo, beta4=blo, beta5=blo, beta6=blo, nu=nlo),
                    nprof = nprofs)
      betasamps <- profileDesign(beta1 = seq(blo, bup, length = 30),
                    upper = c(rho=1, theta=20, beta2=bup, beta3=bup, beta4=bup, beta5=bup, beta6=bup, nu=nup),
                    lower = c(rho=1E-8, theta=1, beta2=blo, beta3=blo, beta4=blo, beta5=blo, beta6=blo, nu=nlo),
                    nprof = nprofs)
      nusamps <- profileDesign(nu = seq(nlo, nup, length = 10),
                    upper = c(rho=1, theta=20, beta1=bup, beta2=bup, beta3=bup, beta4=bup, beta5=bup, beta6=bup),
                    lower = c(rho=1E-8, theta=1, beta1=blo, beta2=blo, beta3=blo, beta4=blo, beta5=blo, beta6=blo),
                    nprof = nprofs)
      starts_epi <- bind_rows(rhosamps, thetasamps, betasamps, nusamps) %>%
        dplyr::mutate(parid = seq_along(rho)) %>%
        dplyr::mutate(gamma=gamma, sigma=sigma, theta0=theta0, alpha=alpha, mu=mu, delta=delta, nu=nu, S.0=S0, E.0=E0, I.0=I0, A.0=A0, R.0=R0) %>%
        dplyr::select(parid, rho, theta, beta1, beta2, beta3, beta4, beta5, beta6, nu, gamma, sigma, theta0, alpha, mu, delta, S.0, E.0, I.0, A.0, R.0)
    }
    
    write_csv(starts_epi, paste0(out_dir, "starts_", mcode, "_epi", nweeks, ifelse(is.na(nsamps), paste0("_nprofs", nprofs), paste0("_nsamps", nsamps)), "_", parseed, ".csv"))
  }  

  return(starts_epi)
}


#' @description explore coverage of starting epidemic parameters
explore_epi_param_coverage <- function(fig_dir){

  parseed <- fix_epi_settings()$epi_seed
  nweeks <- fix_epi_settings()$nweeks
  mcode <- fix_epi_settings()$mcode
  par_dir <- paste0("haiti_ocv_JHUmodel/GeneratedData/", mcode, "/")
  starts_epi <- generate_epi_params(out_dir=par_dir)

  p1 <- qplot(theta, rho, data = starts_epi)
  p2 <- qplot(theta, beta1, data = starts_epi)
  p3 <- qplot(rho, beta1, data = starts_epi)
  plt <- grid.arrange(p1, p2, p3, nrow = 1)
  ggsave(paste0(fig_dir, "start_coverage_", mcode, "_epi", nweeks, "_seed", parseed, ".png"), plt, width = 6, height = 3, units = "in")
  return(plt)
}
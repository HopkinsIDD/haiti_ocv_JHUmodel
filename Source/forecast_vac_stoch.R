rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_fc_params.R")


## -------------------------------------- ##
fcsettings <- fix_fc_settings()
endsettings <- fix_end_settings()
episettings <- fix_epi_settings()
scode <- fcsettings$scode

fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", fcsettings$mcode, "_", scode, "/")
out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", fcsettings$mcode, "_", scode, "/")
out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_fc, showWarnings = FALSE)

# ## -------------------------------------- ##
# ## testing area
#  topn <- 5

## -------------------------------------- ##
## control flow for pomp model and department-vacc initial conditions
if (scode %in% paste0("id", seq(2, 34, by = 4))) { ## 2 dept campaigns
    source("haiti_ocv_JHUmodel/Source/pomp_fc_2dept.R")
    starts <- generate_fc_params(out_dir_epi, out_dir_end, out_dir_fc) %>%
      dplyr::mutate(S1.0 = 0.0, E1.0 = 0.0, I1.0 = 0.0, A1.0 = 0.0, R1.0 = 0.0,
                    S2.0 = 0.0, E2.0 = 0.0, I2.0 = 0.0, A2.0 = 0.0, R2.0 = 0.0) #%>% top_n(topn, wt = loglik_f)
    allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","kappa",
      "S.0","E.0","I.0","A.0","R.0",
      "S1.0","E1.0","I1.0","A1.0","R1.0",
      "S2.0","E2.0","I2.0","A2.0","R2.0",
      "incid.0","N0")
} else if (scode %in% paste0("id", seq(4, 36, by = 4))) { ## 3 dept campaigns
    source("haiti_ocv_JHUmodel/Source/pomp_fc_3dept.R") 
    starts <- generate_fc_params(out_dir_epi, out_dir_end, out_dir_fc) %>%
      dplyr::mutate(S1.0 = 0.0, E1.0 = 0.0, I1.0 = 0.0, A1.0 = 0.0, R1.0 = 0.0,
                    S2.0 = 0.0, E2.0 = 0.0, I2.0 = 0.0, A2.0 = 0.0, R2.0 = 0.0,
                    S3.0 = 0.0, E3.0 = 0.0, I3.0 = 0.0, A3.0 = 0.0, R3.0 = 0.0)
    allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","kappa",
      "S.0","E.0","I.0","A.0","R.0",
      "S1.0","E1.0","I1.0","A1.0","R1.0",
      "S2.0","E2.0","I2.0","A2.0","R2.0",
      "S3.0","E3.0","I3.0","A3.0","R3.0",
      "incid.0","N0")
} else if (scode %in% paste0("id", c(seq(3, 35, by = 4), seq(1, 33, by = 4)))) { ## slow national & fast national campaigns
    source("haiti_ocv_JHUmodel/Source/pomp_fc_national.R")
    starts <- generate_fc_params(out_dir_epi, out_dir_end, out_dir_fc) %>%
      dplyr::mutate(S1.0 = 0.0, E1.0 = 0.0, I1.0 = 0.0, A1.0 = 0.0, R1.0 = 0.0,
                    S2.0 = 0.0, E2.0 = 0.0, I2.0 = 0.0, A2.0 = 0.0, R2.0 = 0.0,
                    S3.0 = 0.0, E3.0 = 0.0, I3.0 = 0.0, A3.0 = 0.0, R3.0 = 0.0,
                    S4.0 = 0.0, E4.0 = 0.0, I4.0 = 0.0, A4.0 = 0.0, R4.0 = 0.0,
                    S5.0 = 0.0, E5.0 = 0.0, I5.0 = 0.0, A5.0 = 0.0, R5.0 = 0.0,
                    S6.0 = 0.0, E6.0 = 0.0, I6.0 = 0.0, A6.0 = 0.0, R6.0 = 0.0,
                    S7.0 = 0.0, E7.0 = 0.0, I7.0 = 0.0, A7.0 = 0.0, R7.0 = 0.0,
                    S8.0 = 0.0, E8.0 = 0.0, I8.0 = 0.0, A8.0 = 0.0, R8.0 = 0.0,
                    S9.0 = 0.0, E9.0 = 0.0, I9.0 = 0.0, A9.0 = 0.0, R9.0 = 0.0,
                    S10.0 = 0.0, E10.0 = 0.0, I10.0 = 0.0, A10.0 = 0.0, R10.0 = 0.0)
    allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","kappa",
      "S.0","E.0","I.0","A.0","R.0",
      "S1.0","E1.0","I1.0","A1.0","R1.0",
      "S2.0","E2.0","I2.0","A2.0","R2.0",
      "S3.0","E3.0","I3.0","A3.0","R3.0",
      "S4.0","E4.0","I4.0","A4.0","R4.0",
      "S5.0","E5.0","I5.0","A5.0","R5.0",
      "S6.0","E6.0","I6.0","A6.0","R6.0",
      "S7.0","E7.0","I7.0","A7.0","R7.0",
      "S8.0","E8.0","I8.0","A8.0","R8.0",
      "S9.0","E9.0","I9.0","A9.0","R9.0",
      "S10.0","E10.0","I10.0","A10.0","R10.0",
      "incid.0","N0")
}


## -------------------------------------- ##
# import settings and data
pop.haiti <- get.haiti.pop()
haiti.dat <- get.mspp.agg.data() 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1) %>%
  dplyr::select(week_end, cases)

## vac covariates
vactab <- make.vactab(t0 = 0, 
                      tmax = nrow(haiti.dat)+fcsettings$horizon+1, 
                      ndept = fcsettings$nd,
                      nweeks = fcsettings$nw,
                      coverage_2dose = fcsettings$c2, 
                      coverage_1dose = fcsettings$c1,
                      first_vac_t = nrow(haiti.dat)+4, 
                      ve_scen = fcsettings$vescen)

# ## -------------------------------------- ##
# ## testing area
# vactab <- vactab %>% dplyr::mutate(num_vacc=0, ve_d1=0, ve_d2=0)

## spline covariates
covartab.all <- make.covartab(0, nrow(haiti.dat)+fcsettings$horizon+1, byt=1, degree=6, nbasis=6, per=52.14) %>%
  dplyr::mutate(time_endfc = ifelse(time < episettings$nweeks, NA, time)) %>%
  dplyr::mutate(time_fc = ifelse(time < max(haiti.dat$week), NA, time)) %>%
  dplyr::mutate(time_end = ifelse(time <= max(haiti.dat$week), time, NA))
covartab.fc <- covartab.all %>%
  dplyr::filter(!is.na(time_fc)) %>%
  dplyr::select(time_fc, contains("seas")) 
covartab.all2 <- covartab.all %>%
  dplyr::select(time, contains("seas")) %>%
  full_join(vactab, by = c("time"))

## -------------------------------------- ##
## models
haiti.mod.fc <- build.fc.mod(pop = pop.haiti,
                      dat = haiti.dat,
                      my.times = "week",
                      covar.times = "time",
                      my.t0 = 0,
                      covar = covartab.all2)
time(haiti.mod.fc) <- covartab.fc$time_fc
timezero(haiti.mod.fc) <- min(covartab.fc$time_fc)

models <- list(haiti.mod.fc)
dim(models) <- c(1,1)
dimnames(models) <- list(c("stoch.forecast"), c("haiti"))

## -------------------------------------- ##
## stochastic forecasts 
bake(file=paste0(out_dir_fc, "/stochfcs_", fcsettings$str_if, ".rds"), seed = fcsettings$seed, {

  foreach(start = iter(starts, by = 'row'),
          .inorder = FALSE, .combine = rbind,
          .packages = c("pomp", "magrittr"),
          .errorhandling = c('remove'),
          .noexport = c(),
          .verbose = TRUE) %dopar%
    {

    tic <- Sys.time()
    print(paste("forecast stoch params vac", start$parid))

    M2 <- models["stoch.forecast","haiti"][[1]]
    coef(M2) <- unlist(start[which(names(start) %in% allpars)])
    
    base_sims <- simulate(M2, nsim = fcsettings$nsims, as.data.frame=TRUE) 

    if (scode %in% paste0("id", seq(2, 34, by = 4))) { ## 2 dept
      fc_states_sims <- base_sims %>%
        dplyr::mutate(Nnv = S+E+I+A+R, 
                      N1 = S1+E1+I1+A1+R1, 
                      N2 = S2+E2+I2+A2+R2, 
                      N = Nnv+N1+N2, 
                      loglik = start$loglik_f) %>%
        dplyr::rename(week = time) %>%
        dplyr::select(sim, week, cases, 
                      S, E, I, A, R, Nnv, 
                      S1, E1, I1, A1, R1, N1, 
                      S2, E2, I2, A2, R2, N2, 
                      N, incid, incidU, incidV, asymV, newV, loglik, foival, Str0, Sout, Sin)
    } else if (scode %in% paste0("id", seq(4, 36, by = 4))) { ## 3 dept
      fc_states_sims <- base_sims %>%
        dplyr::mutate(Nnv = S+E+I+A+R, 
                      N1 = S1+E1+I1+A1+R1, 
                      N2 = S2+E2+I2+A2+R2, 
                      N3 = S3+E3+I3+A3+R3,
                      N = Nnv+N1+N2+N3, 
                      loglik = start$loglik_f) %>%
        dplyr::rename(week = time) %>%
        dplyr::select(sim, week, cases, 
                      S, E, I, A, R, Nnv, 
                      S1, E1, I1, A1, R1, N1, 
                      S2, E2, I2, A2, R2, N2,
                      S3, E3, I3, A3, R3, N3, 
                      N, incid, incidU, incidV, asymV, newV, loglik, foival, Str0, Sout, Sin)
    } else if (scode %in% paste0("id", c(seq(3, 35, by = 4), seq(1, 33, by = 4)))) { ## all dept
      fc_states_sims <- base_sims %>%
        dplyr::mutate(Nnv = S+E+I+A+R, 
                      N1 = S1+E1+I1+A1+R1, 
                      N2 = S2+E2+I2+A2+R2, 
                      N3 = S3+E3+I3+A3+R3, 
                      N4 = S4+E4+I4+A4+R4,
                      N5 = S5+E5+I5+A5+R5, 
                      N6 = S6+E6+I6+A6+R6,
                      N7 = S7+E7+I7+A7+R7, 
                      N8 = S8+E8+I8+A8+R8,
                      N9 = S9+E9+I9+A9+R9, 
                      N10 = S10+E10+I10+A10+R10,
                      N = Nnv+N1+N2+N3+N4+N5+N6+N7+N8+N9+N10, 
                      loglik = start$loglik_f) %>%
        dplyr::rename(week = time) %>%
        dplyr::select(sim, week, cases, 
                      S, E, I, A, R, Nnv, 
                      S1, E1, I1, A1, R1, N1, 
                      S2, E2, I2, A2, R2, N2,
                      S3, E3, I3, A3, R3, N3, 
                      S4, E4, I4, A4, R4, N4, 
                      S5, E5, I5, A5, R5, N5, 
                      S6, E6, I6, A6, R6, N6, 
                      S7, E7, I7, A7, R7, N7, 
                      S8, E8, I8, A8, R8, N8, 
                      S9, E9, I9, A9, R9, N9, 
                      S10, E10, I10, A10, R10, N10,  
                      N, incid, incidU, incidV, asymV, newV, loglik, foival, Str0, Sout, Sin)
    }

    write_csv(fc_states_sims, paste0(out_dir_fc, "fcstates_stoch_", fcsettings$str_if2, "_parid", start$parid, "_seed", fcsettings$seed, ".csv"))

    fc_states_sims 

    } %>% group_by(week) %>%
      summarise(S_med = median(S), E_med = median(E), I_med = median(I), A_med = median(A), R_med = median(R), 
                incid_med = median(incid), N_med = median(N), cases_med = median(cases),
                S_lo = quantile(S, probs = c(.025)), E_lo = quantile(E, probs = c(.025)), I_lo = quantile(I, probs = c(.025)), A_lo = quantile(A, probs = c(.025)), R_lo = quantile(R, probs = c(.025)), 
                incid_lo = quantile(incid, probs = c(.025)), N_lo = quantile(N, probs = c(.025)), cases_lo = quantile(cases, probs = c(.025)),
                S_hi = quantile(S, probs = c(.975)), E_hi = quantile(E, probs = c(.975)), I_hi = quantile(I, probs = c(.975)), A_hi = quantile(A, probs = c(.975)), R_hi = quantile(R, probs = c(.975)), 
                incid_hi = quantile(incid, probs = c(.975)), N_hi = quantile(N, probs = c(.975)), cases_hi = quantile(cases, probs = c(.975)))

})  -> fc_stoch

plot(fc_stoch$week, fc_stoch$cases_med)


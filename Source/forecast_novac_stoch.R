rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_fc_novac_params.R")
source("haiti_ocv_JHUmodel/Source/pomp_fc_novac.R") ## novac

## -------------------------------------- ##
fcsettings <- fix_fc_settings()
endsettings <- fix_end_settings()
episettings <- fix_epi_settings()

fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", fcsettings$mcode, "_", fcsettings$scode, "/")
out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/")
out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_fc, showWarnings = FALSE)

# ## -------------------------------------- ##
# ## testing area
# topn <- 5

## -------------------------------------- ##
# import settings and data
starts <- generate_fc_params(out_dir_epi, out_dir_end, out_dir_fc) #%>% top_n(topn, wt = loglik_f)
pop.haiti <- get.haiti.pop()
haiti.dat <- get.mspp.agg.data() 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1) %>%
  dplyr::select(week_end, cases)

covartab.all <- make.covartab(0, nrow(haiti.dat)+fcsettings$horizon+1, byt=1, degree=6, nbasis=6, per=52.14) %>%
  dplyr::mutate(time_endfc = ifelse(time < episettings$nweeks, NA, time)) %>%
  dplyr::mutate(time_fc = ifelse(time < max(haiti.dat$week), NA, time)) %>%
  dplyr::mutate(time_end = ifelse(time <= max(haiti.dat$week), time, NA))
covartab.fit <- covartab.all %>%
  dplyr::filter(!is.na(time_end)) %>%
  dplyr::select(time_end, contains("seas"))
covartab.fc <- covartab.all %>%
  dplyr::filter(!is.na(time_fc)) %>%
  dplyr::select(time_fc, contains("seas")) 
covartab.all2 <- covartab.all %>%
  dplyr::select(time, contains("seas"))
num_betas <- ncol(covartab.all2)-1

## -------------------------------------- ##
## models

haiti.mod.fit <- build.fc.mod(pop = pop.haiti,
                      dat = haiti.dat,
                      my.times = "week",
                      covar.times = "time",
                      my.t0 = 0,
                      covar = covartab.all2)
haiti.mod.fc <- build.fc.mod(pop = pop.haiti,
                      dat = haiti.dat,
                      my.times = "week",
                      covar.times = "time",
                      my.t0 = 0,
                      covar = covartab.all2)
time(haiti.mod.fc) <- covartab.fc$time_fc
timezero(haiti.mod.fc) <- min(covartab.fc$time_fc)

models <- list(haiti.mod.fit, haiti.mod.fc)
dim(models) <- c(2,1)
dimnames(models) <- list(c("fit", "stoch.forecast"), c("haiti"))

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
    print(paste("forecast stoch params novac", start$parid))

    M2 <- models["stoch.forecast","haiti"][[1]]
    allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","S.0","E.0","I.0","A.0","R.0","incid.0","N0")
    coef(M2) <- unlist(start[which(names(start) %in% allpars)])
    
    fc_states_sims <- simulate(M2, nsim = fcsettings$nsims, as.data.frame=TRUE) %>%
      dplyr::mutate(N = S+E+I+A+R, loglik = start$loglik_f) %>%
      dplyr::rename(week = time) %>%
      dplyr::select(sim, week, cases, S, E, I, A, R, N, incid, loglik, foival, Str0, Sout, Sin) 

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




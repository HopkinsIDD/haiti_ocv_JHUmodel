## rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_end_params.R")
source("haiti_ocv_JHUmodel/Source/pomp_fit_end.R")

# registerDoParallel(cores=2)

## -------------------------------------- ##
episettings <- fix_epi_settings()
endsettings <- fix_end_settings()
fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", endsettings$mcode, "/")
out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_end, showWarnings = FALSE)

nmif <- endsettings$nmif; nparticles <- endsettings$nparticles

## -------------------------------------- ##
## import settings and data
starts <- generate_end_params(out_dir_epi, out_dir_end)
pop.haiti <- get.haiti.pop()

haiti.dat <- get.mspp.agg.data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1)

covartab.all <- make.covartab(0, nrow(haiti.dat)+1, byt=1, degree=6, nbasis=6, per=52.14)
covartab <- covartab.all[(episettings$nweeks+1):nrow(covartab.all),] %>%
  dplyr::mutate(time_end = seq_along(time)-1) %>%
  dplyr::select(time_end, contains("seas"))
num_betas <- ncol(covartab.all)-1

haiti.mod.end <- build.end.mod(pop = pop.haiti, 
                              dat = haiti.dat.end,
                              my.times = "week_end",
                              covar.times = "time_end",
                              covar = covartab)

## -------------------------------------- ##
## perform end trajectory iterated filtering
bake(file = paste0(out_dir_end, "iffits_", endsettings$str_if, ".rds"), seed = endsettings$seed, {

    if(endsettings$passTm == "true"){
      prof_tm <- readRDS(paste0(out_dir_end, "/tmfits_", endsettings$str_tm, ".rds")) 
      starts_if <- prof_tm %>%
        dplyr::filter(is.finite(loglik)) %>%
        # dplyr::filter(parid == par) %>%
        # dplyr::filter(is.finite(loglik) & conv %in% c(0,1)) %>%
        # dplyr::sample_n(ix) %>%
        dplyr::select(-model, -loglik, -conv, -etime)
    } else{
      parplt <- explore_end_inits_coverage(fig_dir = fig_dir)
      starts_if <- starts
    }
    

    foreach(start = iter(starts_if, by = 'row'),
            .combine = rbind, .inorder = FALSE,
            .packages = c("pomp", "magrittr"),
            .errorhandling = c('remove'),
            .export = c("haiti.dat.epi", "num_betas", "pop.haiti", "covartab"),
            .noexport = c(),
            .verbose = TRUE) %dopar% 
    {
      tic <- Sys.time()  
      print(paste(start$parid, "iterated filtering"))

      po <- haiti.mod.end
      allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","S.0","E.0","I.0","A.0","R.0","incid.0","N0")
      coef(po) <- unlist(start[which(names(start) %in% allpars)])
      est_params <- c(paste0("beta", 1:num_betas), "rho", "theta", "nu")
      
      pfilt_init <- replicate(n=20, logLik(pfilter(po, Np=10)))
      print(paste(start$parid, "init loglik", mean(pfilt_init)))

      if (coef(po,"E.0")==0) coef(po,"E.0") <- 1e-9
      if (coef(po,"I.0")==0) coef(po,"I.0") <- 1e-9

      mf.mod <- mif2(po, Nmif=nmif,
                rw.sd = rw.sd(
                          beta1=.05,
                          beta2=.05,
                          beta3=.05,
                          beta4=.05,
                          beta5=.05,
                          beta6=.05,
                          theta=.005,
                          rho = .01,
                          nu=.04),
                # ivps = c("E.0","I.0"),
                Np = nparticles,
                cooling.type = "hyperbolic",
                cooling.fraction.50 = 0.05,
                transform = TRUE,
                verbose = FALSE)
      print(paste("mif loglik after", nmif, "iter", logLik(pfilter(mf.mod, Np=25))))
      mf.mod <- continue(mf.mod, Nmif = nmif)
      print(paste("mif loglik after", nmif+nmif, "iter", logLik(pfilter(mf.mod, Np=25))))
        
      est_rho <- coef(mf.mod, "rho")
      est_theta <- coef(mf.mod, "theta")
      fit_states_sims <- simulate(po, params = coef(mf.mod), nsim = 50, as.data.frame = TRUE, transform = TRUE) %>%
        dplyr::select(-week) %>%
        dplyr::mutate(parid = start$parid, N = S+E+I+A+R) %>%
        dplyr::rename(week = time) %>%
        dplyr::select(parid, sim, week, S, E, I, A, R, incid, N, cases)      

      fit_states <- fit_states_sims %>%
        group_by(week) %>%
        summarise(S_med = median(S), E_med = median(E), I_med = median(I), A_med = median(A), R_med = median(R), incid_med = median(incid), N_med = median(N), cases_med = median(cases),
                  S_lo = quantile(S, probs = c(.025)), E_lo = quantile(E, probs = c(.025)), I_lo = quantile(I, probs = c(.025)), A_lo = quantile(A, probs = c(.025)), R_lo = quantile(R, probs = c(.025)), incid_lo = quantile(incid, probs = c(.025)), N_lo = quantile(N, probs = c(.025)), cases_lo = quantile(cases, probs = c(.025)),
                  S_hi = quantile(S, probs = c(.975)), E_hi = quantile(E, probs = c(.975)), I_hi = quantile(I, probs = c(.975)), A_hi = quantile(A, probs = c(.975)), R_hi = quantile(R, probs = c(.975)), incid_hi = quantile(incid, probs = c(.975)), N_hi = quantile(N, probs = c(.975)), cases_hi = quantile(cases, probs = c(.975)))
        
      write_csv(fit_states, paste0(out_dir_end, "fitstates_if_", endsettings$str_if2, "_parid", start$parid, "_seed", endsettings$seed, ".csv"))

      pf.lik <- replicate(10, pfilter(mf.mod, Np=25, max.fail=Inf))
      ll <- sapply(pf.lik, logLik)
      ll <- logmeanexp(ll, se = TRUE)
      nfail <- sapply(pf.lik, getElement, "nfail")
      
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      dummy <- data.frame(model = "if", 
                parid = start$parid,
                as.list(coef(mf.mod)),
                loglik = ll[1],
                loglik.se = ll[2], 
                nfail.min = min(nfail),
                nfail.max = max(nfail), ## best if nfail.max==0
                etime = as.numeric(etime))
      write_csv(dummy, paste0(out_dir_end, "parest_if_", endsettings$str_if2, "_parid", start$parid, "_seed", endsettings$seed, ".csv"))

      rm(fit_states, fit_states_sims, mf.mod)
      gc()
      
      dummy
    
    } 

})  -> prof_if


## -------------------------------------- #### -------------------------------------- ##
## get a sample index for subsequent checks
conv_if <- prof_if %>% filter(nfail.max==0)
print(paste(nrow(conv_if), "valid iterated filtering estimates"))
ismp <- sample(1:nrow(conv_if), size=1)
## -------------------------------------- ##
## check beta and rough R0 estimate
betas <- conv_if[ismp,] %>% dplyr::select(contains("beta")) %>% as.matrix
bspline <- covartab %>% dplyr::select(contains("seas")) %>% as.matrix %>% t
mybeta <- betas %*% bspline %>% as.vector
plot(mybeta, type = "l")
## adjustment for mixing coefficient
I0s <- conv_if[ismp,] %>% dplyr::select(I.0) %>% unlist %>% median
nus <- conv_if[ismp,] %>% dplyr::select(nu) %>% unlist %>% median
print(paste("R0 ranges (init, min, max):", (I0s^nus)/I0s*mybeta[1]*2/7, (I0s^nus)/I0s*min(mybeta)*2/7, (I0s^nus)/I0s*max(mybeta)*2/7))
## -------------------------------------- ##
## check one parameter set2
coefnames <- names(starts)[-which(names(starts)=="parid")]
params <- parmat(unlist(conv_if[ismp, coefnames]), nrep=1)
simdat <- simulate(haiti.mod.epi,
          params=params, nsim=50,
          as.data.frame = TRUE,
          transform = TRUE) %>%
      dplyr::rename(week = time) %>%
      group_by(week) %>%
      dplyr::summarise(cases = mean(cases), incid = mean(incid)) %>%
      dplyr::left_join(haiti.dat.epi %>% dplyr::rename(orig=cases), by = c("week"))      

fitplt <- ggplot(simdat, aes(x = week)) +
      geom_line(aes(y = cases), colour = "blue", alpha = 0.35) +
      geom_point(aes(y = orig), size = 1) +
      annotate("text", x = 100, y = 20000, label = unlist(conv_if[ismp, "loglik"])) +
      theme_bw() +
      scale_y_continuous("cases", limits = c(0, 30000), oob = scales::squish)
print(fitplt)    
## -------------------------------------- ##
## check range of important params
plot(conv_if$theta, conv_if$rho)
plot(conv_if$beta1, conv_if$rho)
plot(conv_if$beta1, conv_if$rho)
plot(log(conv_if$E.0+conv_if$I.0), conv_if$rho)

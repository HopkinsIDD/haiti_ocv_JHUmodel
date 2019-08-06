# rm(list=ls())
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

## -------------------------------------- ##
## import settings and data
starts <- generate_end_params(out_dir_epi, out_dir_end)
parplt <- explore_end_inits_coverage(fig_dir = fig_dir)

pop.haiti <- get.haiti.pop()

haiti.dat <- get.mspp.agg.data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1)

covartab.all <- make.covartab(0, nrow(haiti.dat)+1, byt=1, degree=6, nbasis=6, per=52.14)
covartab <- covartab.all[(episettings$nweeks+1):nrow(covartab.all),] %>%
  dplyr::mutate(time_end = seq_along(time)-1)
num_betas <- ncol(covartab.all)-1

haiti.mod.end <- build.end.mod(pop = pop.haiti, 
                              dat = haiti.dat.end,
                              my.times = "week_end",
                              covar.times = "time_end",
                              covar = covartab)

## -------------------------------------- ##
## perform end trajectory matching fits
bake(file = paste0(out_dir_end, "/tmfits_", endsettings$str_tm, ".rds"), seed = endsettings$seed, {

    foreach(start = iter(starts, by = 'row'),
            .combine = rbind, .inorder = FALSE,
            .packages = c("pomp", "magrittr"),
            .errorhandling = c('remove'),
            .export = c("haiti.dat.end", "num_betas", "pop.haiti", "covartab", "endsettings"),
            .noexport = c(),
            .verbose = TRUE) %dopar% 
    {
      tic <- Sys.time()
      print(paste(start$parid, "trajmatch"))

      tm <- haiti.mod.end
      allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","S.0","E.0","I.0","A.0","R.0","incid.0","N0")
      coef(tm) <- unlist(start[which(names(start) %in% allpars)])
      est_params <- c(paste0("beta", 1:num_betas), "rho", "theta", "nu")
      
      pfilt_init <- replicate(n=50, logLik(pfilter(tm, N=10)))
      print(paste(start$parid, "init loglik", pfilt_init))

      tm.mod <- traj.match(tm, est = est_params, transform = TRUE)
      if (coef(tm.mod,"E.0")==0) coef(tm.mod,"E.0") <- 1e-9
      if (coef(tm.mod,"I.0")==0) coef(tm.mod,"I.0") <- 1e-9
      tm.mod <- traj.match(tm.mod, est = est_params, transform = TRUE, maxit = 1E5)
      
      toc <- Sys.time()
      etime <- toc-tic
      units(etime) <- "hours"

      # export fitted states
      est_rho <- coef(tm.mod, "rho")
      est_theta <- coef(tm.mod, "theta")
      fit_states <- tbl_df(t(tm.mod@states)) %>% 
        dplyr::mutate(week = seq_along(S), parid = start$parid) %>%
        dplyr::select(parid, week, S, E, I, A, R, incid) %>%
        dplyr::mutate(N = S+E+I+A+R) %>%
        rowwise %>%
        dplyr::mutate(est_cases = rnbinom(1, mu = est_rho*incid, size = est_theta)) %>% ungroup
      write_csv(fit_states, paste0(out_dir_end, "fitstates_tm_", endsettings$str_tm2, "_parid", start$parid, "_seed", episettings$seed, ".csv"))

      dummy <- data.frame(model = "tm",
                parid = start$parid,
                as.list(coef(tm.mod)),
                loglik = logLik(tm.mod),
                conv = tm.mod$convergence,
                etime = as.numeric(etime))
      write_csv(dummy, paste0(out_dir_end, "parest_tm_", endsettings$str_tm2, "_parid", start$parid, "_seed", episettings$seed, ".csv"))

      rm(fit_states, tm.mod)
      gc()
      
      dummy
    
    } %>% dplyr::mutate(sum = S.0+E.0+I.0+A.0+R.0,
                S.0 = round(N0*S.0/sum),
                E.0 = round(N0*E.0/sum),
                I.0 = round(N0*I.0/sum),
                A.0 = round(N0*A.0/sum),
                R.0 = round(N0*R.0/sum)) %>%
          dplyr::select(-sum) %>%
          unique() 

})  -> prof_tm


## -------------------------------------- #### -------------------------------------- ##
## get a sample index for subsequent checks
conv_tm <- prof_tm %>% filter(conv %in% c(0,1))
print(paste(nrow(conv_tm), "potential valid estimates"))
ismp <- sample(1:nrow(conv_tm), 1)
## -------------------------------------- ##
## check beta and rough R0 estimate
betas <- conv_tm[ismp,] %>% dplyr::select(contains("beta")) %>% as.matrix
bspline <- covartab %>% dplyr::select(contains("seas")) %>% as.matrix %>% t
mybeta <- betas %*% bspline %>% as.vector
plot(mybeta, type = "l")
## adjustment for mixing coefficient
I0s <- conv_tm[ismp,] %>% dplyr::select(I.0) %>% unlist %>% median
nus <- conv_tm[ismp,] %>% dplyr::select(nu) %>% unlist %>% median
print(paste("R0 ranges (init, min, max):", (I0s^nus)/I0s*mybeta[1]*2/7, (I0s^nus)/I0s*min(mybeta)*2/7, (I0s^nus)/I0s*max(mybeta)*2/7))
## -------------------------------------- ##
## check one parameter set2
coefnames <- names(starts)[-which(names(starts) %in% c("parid", "epi_mcode", "loglik"))]
params <- parmat(unlist(conv_tm[ismp, coefnames]), nrep=1)
simdat <- simulate(haiti.mod.end,
          params=params, nsim=50,
          as.data.frame = TRUE,
          transform = TRUE) %>%
      dplyr::rename(week = time) %>%
      # group_by(week) %>%
      # dplyr::summarise(cases = mean(cases), incid = mean(incid)) %>%
      dplyr::left_join(haiti.dat.end %>% dplyr::rename(orig=cases), by = c("week"))      

fitplt <- ggplot(simdat, aes(x = week)) +
      geom_line(aes(y = cases), colour = "blue", alpha = 0.35) +
      geom_point(aes(y = orig), size = 1) +
      annotate("text", x = 100, y = 20000, label = unlist(conv_tm[ismp, "loglik"])) +
      theme_bw() +
      scale_y_continuous("cases", limits = c(0, 30000), oob = scales::squish)
print(fitplt)    
## -------------------------------------- ##
## check range of important params
plot(conv_tm$theta, conv_tm$rho)
plot(conv_tm$beta1, conv_tm$rho)
plot(conv_tm$beta1, conv_tm$rho)
plot(log(conv_tm$E.0+conv_tm$I.0), conv_tm$rho)

# rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_epi_params.R")
source("haiti_ocv_JHUmodel/Source/pomp_fit_epi.R")

# registerDoParallel(cores=2)

## -------------------------------------- ##
episettings <- fix_epi_settings()
fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", episettings$mcode, "/")
out_dir <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

## -------------------------------------- ##
## import settings and data
starts <- generate_epi_params(out_dir = out_dir, priors = episettings$priors)
parplt <- explore_epi_param_coverage(fig_dir = fig_dir)

pop.haiti <- get.haiti.pop()

haiti.dat <- get.mspp.agg.data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 

covartab <- make.covartab(0, nrow(haiti.dat.epi)+1, byt=1, degree=6, nbasis=6, per=52.14)
num_betas <- ncol(covartab)-1

haiti.mod.epi <- build.epi.mod(pop = pop.haiti, 
                              dat = haiti.dat.epi,
                              covar = covartab)

## -------------------------------------- ##
## filters for testing
# starts <- starts[which(starts$parid==84),]
# starts <- data.frame(parid=1001, rho=.98, theta=7.23, beta1=1.67, beta2=5.77, beta3=.717, beta4=5.96, beta5=1.58, beta6=5.7, gamma=3.5, sigma=5.0, theta0=0, alpha=.0024, mu=.0002, delta=.0001, S.0=.999, E.0=.0005, I.0=.000000094, A.0=0, R.0=0)
# starts <- data.frame(parid=1001, rho=.31, theta=3.79, beta1=5.51, beta2=2.95, beta3=4.07, beta4=2.59, beta5=5.7, beta6=1.73, gamma=3.5, sigma=5.0, theta0=0, alpha=.0024, mu=.0002, delta=.0001, S.0=.999, E.0=1.39E-14, I.0=5.79E-4, A.0=0, R.0=0)
# starts <- data.frame(parid=1001, rho=.31, theta=3.79, beta1=4.08, beta2=4.19, beta3=2.88, beta4=3.72, beta5=4.32, beta6=3.36, gamma=3.5, sigma=5.0, theta0=0, alpha=.002397, mu=.0002, delta=.0001, S.0=.99936, E.0=3.899E-12, I.0=6.4E-4, A.0=0, R.0=0)


## -------------------------------------- ##
## perform epi trajectory matching fits
bake(file = paste0(out_dir, "/tmfits_", episettings$str_tm, ".rds"), seed = episettings$seed, {

    foreach(start = iter(starts, by = 'row'),
            .combine = rbind, .inorder = FALSE,
            .packages = c("pomp", "magrittr"),
            .errorhandling = c('remove'),
            .export = c("haiti.dat.epi", "num_betas", "pop.haiti", "covartab", "episettings"),
            .noexport = c(),
            .verbose = TRUE) %dopar% 
    {
      tic <- Sys.time()
      # parid=2
      # start <- starts[parid,]
      print(paste(start$parid, "trajmatch"))
      tm <- haiti.mod.epi
      allpars <- c("rho","theta","beta1","beta2","beta3","beta4","beta5","beta6","gamma","sigma","theta0","alpha","mu","delta","nu","S.0","E.0","I.0","A.0","R.0")
      coef(tm) <- unlist(start[which(names(start) %in% allpars)])
      est_params <- c(paste0("beta", 1:num_betas), "rho", "theta", "nu", "I.0", "E.0")
      
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
      write_csv(fit_states, paste0(out_dir, "fitstates_tm_", episettings$str_tm2, "_parid", start$parid, "_seed", episettings$seed, ".csv"))

      dummy <- data.frame(model = "tm",
                parid = start$parid,
                as.list(coef(tm.mod)),
                loglik = logLik(tm.mod),
                conv = tm.mod$convergence,
                etime = as.numeric(etime))
      write_csv(dummy, paste0(out_dir, "parest_tm_", episettings$str_tm2, "_parid", start$parid, "_seed", episettings$seed, ".csv"))
      
      rm(fit_states, tm.mod)
      gc()
      
      dummy
    
    } %>% dplyr::mutate(sum = S.0+E.0+I.0+A.0+R.0,
                S.0 = round(pop.haiti*S.0/sum),
                E.0 = round(pop.haiti*E.0/sum),
                I.0 = round(pop.haiti*I.0/sum),
                A.0 = round(pop.haiti*A.0/sum),
                R.0 = round(pop.haiti*R.0/sum)) %>%
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
coefnames <- names(starts)[-which(names(starts)=="parid")]
params <- parmat(unlist(conv_tm[ismp, coefnames]), nrep=1)
simdat <- simulate(haiti.mod.epi,
          params=params, nsim=50,
          as.data.frame = TRUE,
          transform = TRUE) %>%
      dplyr::rename(week = time) %>%
      # group_by(week) %>%
      # dplyr::summarise(cases = mean(cases), incid = mean(incid)) %>%
      dplyr::left_join(haiti.dat.epi %>% dplyr::rename(orig=cases), by = c("week"))      

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

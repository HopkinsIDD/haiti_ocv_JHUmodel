## functions to prepare data and build pomp models - shared across vacc implementations
## addresses compilation errors
library(devtools)
library(data.table)
find_rtools()
has_devel()

#### data utils ###########################################

## -------------------------------- ##
## A few project specific functions ##
get.mspp.dept.data <- function(){

  ## Saturday dates represent last day of week cumulatively
  allDat <- read_csv("haiti_ocv_JHUmodel/Data/haiti-data-from-2010-10-to-2019-01.csv", skip = 1, col_names = c("date_sat_orig", "report", "Artibonite", "Centre", "Grand_Anse", "Nippes", "Nord", "Nord_Est", "Nord_Ouest", "Ouest", "Sud", "Sud_Est"), col_types = "cciiiiiiiiii")

  splitDate <- strsplit(allDat$date_sat_orig, "-")
  setattr(splitDate[[1]], 'names', c("year", "month", "day"))
  dateDf <- tbl_df(as.data.frame(do.call(rbind, splitDate))) %>%
    dplyr::mutate(month = as.character(month)) %>%
    dplyr::mutate(day = as.character(day)) %>%
    dplyr::mutate(year = as.character(year)) %>%
    dplyr::mutate(month = ifelse(nchar(month)==1, paste0("0", month), month)) %>%
    dplyr::mutate(day = ifelse(nchar(day)==1, paste0("0", day), day)) %>%
    dplyr::mutate(date_sat = as.Date(paste(year, month, day, sep = "-")))
  
  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by=7))
    
  cleanDat <- allDat %>% 
    dplyr::mutate(date_sat = as.Date(dateDf$date_sat)) %>%
    dplyr::select(-date_sat_orig) %>%
    full_join(fullDateVec, by = c("date_sat")) %>%
    dplyr::arrange(date_sat) %>%
    dplyr::mutate(week = seq_along(date_sat)) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)

  return(cleanDat)
} 

get.mspp.dept.agg.data <- function(departmentDat, tscale){

  cleanDat <- departmentDat %>%
    dplyr::group_by(!!tscale) %>% ## "day" or "week"
    dplyr::summarise(cases = sum(cases, na.rm=TRUE))

  return(cleanDat)
}


## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.5,0.975)) {
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=2)
  rval$y
}

## --------------------------------- ##
## construct seasonal B-spline terms ##
## --------------------------------- ##

make.covartab <- function(t0, tmax, byt, nbasis=3, degree=3, per=52.18){

    tbasis <- seq(from=t0,to=tmax,by=byt)

    covartab <-  data.frame(cbind(time=tbasis,
                                  periodic.bspline.basis(
                                      x=tbasis,
                                      nbasis=nbasis,
                                      degree=degree,
                                      period=per,
                                      names="seas%d"
                                      )))

    return(covartab)
}


make.vactab <- function(t0=0, tmax, byt=1, ndept=10, nweeks=10, coverage_2dose, coverage_1dose, first_vac_t, ve_scen){

  ## deployment scenario 1 (fast national): one departmental campaign every 10 weeks for all depts

  tbasis <- seq(from=t0,to=tmax,by=byt)
  time_check <- c()
  for(i in 1:ndept) {
    time_check <- c(time_check, rep(0, nweeks-2), rep(i, 2))
  } ## 10 week vacc campaigns (nweeks)

  ## number of vaccines per week by department
  pop_dept <- data.frame(ocv_order = 1:10, dept = c("Centre", "Artibonite", "Ouest", "Nord Ouest", "Nord", "Sud", "Nippes", "Nord Est", "Sud Est", "Grand'Anse"), pop = c(746236, 1727524, 4029705, 728807, 1067177, 774976, 342525, 393967, 632601, 468301)) %>%
    dplyr::mutate(num_vacc = (coverage_2dose+coverage_1dose)*pop/1) ## pulse vaccinees in last 1 week of campaign

  ## create dataframe with number vaccinated for each campaign
  vactab <- data.frame(time = tbasis, vac_tcheck = 0)
  vactab[which(vactab$time %in% first_vac_t:(first_vac_t+length(time_check)-1)),]$vac_tcheck <- time_check
  vactab2 <- left_join(vactab, pop_dept %>% dplyr::select(ocv_order, num_vacc), by = c("vac_tcheck"="ocv_order")) %>%
    dplyr::mutate(num_vacc = ifelse(is.na(num_vacc), 0, round(num_vacc))) #%>%
    ## for some reason it works when you put all vaccinees into 2 consecutive time points (seemingly would double the number vaccinated, but it doesn't) in the covariates
    ## tried adding a very small value to get num_vacc working correctly -- number vaccinated was not corrected
    # group_by(vac_tcheck) %>%
    # dplyr::mutate(num_vacc = ifelse(time==min(time), 1E-6, num_vacc)) %>%
    # ungroup

  ## ve decay after x weeks
  veDecay_mo <- read_csv("haiti_ocv_JHUmodel/Data/ve_decay_bymonth.csv")
  veDecay_adult <- bind_rows(veDecay_mo, veDecay_mo, veDecay_mo, veDecay_mo) %>% 
    arrange(month) %>%
    dplyr::select(!!ve_scen) %>% 
    unlist %>% unname
  ## adjust for lower VE in U5 population (U5 VE is 0.4688*adult VE; roughly 11% of the population is 0-4 years old according to UN World Population Prospects 2017)
  ## adjust ve for age (population VE = adult VE * (1-(1-0.4688)*proportion under-5))
  veDecay <- veDecay_adult * (1-(1-0.4688)*0.11)
  ## adjust for one-dose decay after 52 weeks
  veDecay[53:length(veDecay)] <- coverage_2dose/(coverage_2dose+coverage_1dose)*veDecay[53:length(veDecay)]

  ## add vaccine immunity decay, time shifted for each campaign
  decay_times <- vactab2 %>%
    dplyr::filter(vac_tcheck > 0) %>%
    group_by(vac_tcheck) %>%
    dplyr::filter(time == max(time)) %>% 
    dplyr::select(-num_vacc) %>%
    ungroup
  for(i in 1:ndept) {
    decay_start <- decay_times %>%
      dplyr::filter(vac_tcheck == i) %>% 
      dplyr::select(time) %>% unlist %>% unname
    decay_df <- tbl_df(data.frame(time = (seq_along(veDecay)+decay_start))) %>%
      dplyr::mutate(!!paste0("ve_d", i) := veDecay)
    vactab2 <- left_join(vactab2, decay_df, by = c("time")) 
  }

  vactab3 <- vactab2 %>%
    dplyr::mutate_at(vars(contains("ve_")), funs(ifelse(is.na(.), 0, .)))

  return(vactab3)
}


plot.covartab <- function(covartab_dat, paramVec){

  pltDat <- covartab_dat %>%
    gather(varname, value, contains("seas"))

  pltBsplines <- ggplot(pltDat, aes(x = time, y = value, group = varname)) +
    geom_line(aes(colour = varname)) +
    ggtitle("B splines over time")
  plot(pltBsplines)

  subVec <- paramVec[which(grepl("beta", names(paramVec)))]
  params <- data.frame(paramname = names(subVec), varname = paste0("seas", 1:length(subVec)), param = subVec)
  
  pltDat2 <- left_join(pltDat, params, by = c("varname")) %>%
    dplyr::mutate(intermediate = param * value) %>%
    dplyr::select(time, intermediate) %>%
    group_by(time) %>%
    summarise(beta_t = sum(intermediate))

  pltBeta <- ggplot(pltDat2, aes(x = time, y = beta_t)) +
    geom_line()
  print(pltBeta)

  return(pltBeta)  
}

calc.foi <- function(covartab_dat, params_df, fitdat, pop){

  ## fitdat should be a summary of I&A across fits
  pltDat <- covartab_dat %>%
    gather(varname, value, contains("seas"))

  pltDat2 <- map_dfr(1:nrow(params_df), function(i){

    paramVec <- as.vector(params_df[i,])
    subVec <- unlist(paramVec[which(grepl("beta", names(paramVec)))])
    params <- data.frame(varname = paste0("seas", 1:length(subVec)), param = subVec, stringsAsFactors = FALSE)
    
    left_join(pltDat, params, by = c("varname")) %>%
      dplyr::mutate(intermediate = param * value) %>%
      dplyr::select(time, intermediate) %>%
      group_by(time) %>%
      summarise(beta_t = sum(intermediate)) %>%
      dplyr::mutate(pid = i, iota = paramVec$iota) %>%
      full_join(fitdat, by = c("time")) %>%
      dplyr::filter(time<=404) %>%
      dplyr::mutate(foi = iota + (I_med+A_med)*beta_t/pop) %>%
      dplyr::mutate(foi_up = iota + (I_up+A_up)*beta_t/pop) %>%
      dplyr::mutate(foi_low = iota + (I_low+A_low)*beta_t/pop) 
  }) 
  
  return(pltDat2)

}

## -------------------------------------- ##
## tools to find good starting parameters ##
## -------------------------------------- ##
match_trajectories <- function(haiti_mod_data, param_vec, lhs_traj, rdm_seed = 0, lhs_row = 0){
  
  num_betas <- length(param_vec)-11
  est_params <- c(paste0("beta", 1:num_betas), "rho", "iota", "I.0", "E.0")
  tm.haiti <- traj.match(haiti_mod_data,
                                start = param_vec,
                                est = est_params,
                                method = "Nelder-Mead",
                                reltol = 1e-8,
                                maxit = 15000,
                                transform = TRUE
                                )

        

    if(lhs_traj){

      print(summary(tm.haiti))
      print(paste("orig guess lik", logLik(pfilter(haiti_mod_data, params = param_vec, Np=10000))))
      ll_tm <- logLik(pfilter(haiti_mod_data, params=coef(tm.haiti), Np=10000))
      print(paste("traj lik", ll_tm))

      sim.haiti.tm <- simulate(tm.haiti,
                                params=coef(tm.haiti),
                                nsim=300,
                                transform=TRUE)
      
      png(paste0("haiti_ocv_JHUmodel/figures/find_starting_params/trajmatch_seed", rdm_seed, "_lhs", lhs_row, ".png"), width = 600, height = 600)
        plot(haiti.dat, main = paste("traj match", row, "loglik", tm))
        for (i in 1:300) {
            lines(sim.haiti.tm[[i]]@data[1,], lty=2, col=addalpha(3,.05)) # 
        }
        dev.off()
    }
        

  return(tm.haiti)
}

#### summarizing utils ###########################################
get_sim_quants <- function(pomp.sim.output){
  sim.quants <- apply(sapply(pomp.sim.output,function(x) x@data[1,]),1,function(y) quantile(y,c(.025,.5,.975)))
  return(sim.quants)
}

get_sim_mean <- function(pomp.sim.output){
  sim.mean <- apply(sapply(pomp.sim.output,function(x) x@data[1,]),1,mean)
  return(sim.mean)
}

calculate_cases_averted <- function(pomp.sim.output.vacc, pomp.sim.output.novacc){
  cases.mx.vacc <- sapply(pomp.sim.output.vacc,function(x) x@data[1,])
  cases.mx.novacc <- sapply(pomp.sim.output.novacc,function(x) x@data[1,])
  
  ## cumulative sum of cases for each simulation
  cases.cumsum.vacc <- apply(cases.mx.vacc, 2, sum)
  cases.cumsum.novacc <- apply(cases.mx.novacc, 2, sum)
  cases.averted <- cases.cumsum.novacc-cases.cumsum.vacc
  
  return(list(mn.casesAverted = mean(cases.averted), quant.casesAverted = mean(cases.averted)))
}

prepare_tsPlot_data <- function(pomp.sim.output.vacc, pomp.sim.output.novacc, haiti.dat, t.vac){
  mn.vacc <- get_sim_mean(pomp.sim.output.vacc)
  mn.novacc <- get_sim_mean(pomp.sim.output.novacc)
  quant.vacc <- get_sim_quants(pomp.sim.output.vacc)
  quant.novacc <- get_sim_quants(pomp.sim.output.novacc)

  vacc_df <- data.frame(tstep = (t.vac+1):(length(mn.vacc)+t.vac), scenario = "vacc", mn = mn.vacc, q_025 = quant.vacc[1,], q_5 = quant.vacc[2,], q_975 = quant.vacc[3,])
  novacc_df <- data.frame(tstep = (t.vac+1):(length(mn.vacc)+t.vac), scenario = "no-vacc", mn = mn.novacc, q_025 = quant.novacc[1,], q_5 = quant.novacc[2,], q_975 = quant.novacc[3,])
  haiti_df <- data.frame(tstep = haiti.dat$week, scenario = "reported", mn = haiti.dat$cases, q_025 = NA, q_5 = NA, q_975 = NA)

  df <- bind_rows(vacc_df, novacc_df, haiti_df)
  return(df)
}

#### model diagnostics ###########################################
estimate_R0 <- function(no_vac_fc, pop){
  ## Rt is approx 1 at equilibrium point
  S_prop_mn <- no_vac_fc %>% 
    summarise(S_prop = mean(S)/pop) %>%
    unlist
  return(1/S_prop_mn)
}

## calculate VE from simulation data
estimate_VE <- function(scen_fc, novac_fc, num_weeks, pop){

  ## time after vaccination
  time0 <- scen_fc %>%
    dplyr::filter(vac_tcheck != 0) %>%
    dplyr::filter(time == max(time)) %>%
    distinct(time) %>%
    dplyr::mutate(time = time+1) %>% unlist

  ## vac scenario data
  scen_fc_vacpts <- transmute_all(scen_fc, funs(ifelse(is.na(.), 0, .))) %>%
    dplyr::filter(time %in% (time0:(time0+num_weeks))) %>%
    dplyr::mutate(incidU = incid-incidV) %>%
    group_by(time) %>%
    summarise_all(funs(mean(., na.rm=TRUE))) %>%
    ungroup %>%
    dplyr::mutate(allV = S1+S2+E1+E2+I1+I2+A1+A2+R1+R2+S3+E3+I3+A3+R3+S4+E4+I4+A4+R4+S5+E5+I5+A5+R5+S6+E6+I6+A6+R6+S7+E7+I7+A7+R7+S8+E8+I8+A8+R8+S9+E9+I9+A9+R9+S10+E10+I10+A10+R10,
                  allU = S+E+I+A+R) %>%
    summarise(incid = sum(incid), incidU = sum(incidU), incidV = sum(incidV), asymV = sum(asymV), allV = first(allV), allU = first(allU)) %>%
    dplyr::mutate(all = allV + allU) %>%
    dplyr::mutate(pV = incidV/allV, pU = incidU/allU, pVscen = incid/all)
  
  ## novac scenario data
  novac_fc_vacpts <- novac_fc %>%
    dplyr::filter(time %in% (time0:(time0+num_weeks))) %>%
    group_by(time) %>%
    summarise_all(funs(mean(., na.rm=TRUE))) %>%
    ungroup %>%
    dplyr::mutate(all = S+E+I+A+R) %>%
    summarise(incid = sum(incid), all = first(all)) %>%
    dplyr::mutate(pUscen = incid/all)

  
  ## calculate direct, indirect, total and overall ve
  ve_dir <- 1 - (scen_fc_vacpts$pV/scen_fc_vacpts$pU)
  ve_indir <- 1 - (scen_fc_vacpts$pU/novac_fc_vacpts$pUscen)
  ve_tot <- 1 - (scen_fc_vacpts$pV/novac_fc_vacpts$pUscen)
  ve_ov <- 1 - (scen_fc_vacpts$pVscen/novac_fc_vacpts$pUscen)
  
  return(list(ve_dir = ve_dir, ve_indir = ve_indir, ve_tot = ve_tot, ve_ov = ve_ov))
} 

compare_I <- function(scen_fc, novac_fc, num_weeks, pop){
  
  ## time after vaccination
  time0 <- scen_fc %>%
    dplyr::filter(vac_tcheck != 0) %>%
    dplyr::filter(time == max(time)) %>%
    distinct(time) %>%
    dplyr::mutate(time = time+1) %>% unlist

  ## vac scenario data
  scen_fc_vacpts <- transmute_all(scen_fc, funs(ifelse(is.na(.), 0, .))) %>%
    dplyr::filter(time %in% (time0:(time0+num_weeks))) %>%
    dplyr::mutate(incidU = incid-incidV) %>%
    group_by(time) %>%
    summarise_all(funs(mean(., na.rm=TRUE))) %>%
    ungroup %>%
    dplyr::mutate(E_sc = E+E1+E2+E3+E4+E5+E6+E7+E8+E9+E10, 
                  I_sc = I+I1+I2+I3+I4+I5+I6+I7+I8+I9+I10
                  ) %>%
    dplyr::select(time, E, E1, E2, I, I1, I2, incid, incidU, incidV, E_sc, I_sc, E3, I3, E4, I4, E5, I5, E6, I6, E7, I7, E8, I8, E9, I9, E10, I10
      ) %>%
    dplyr::rename(incid_sc = incid, E0 = E, I0 = I) 

  ## novac scenario data
  novac_fc_vacpts <- novac_fc %>%
    dplyr::filter(time %in% (time0:(time0+num_weeks))) %>%
    group_by(time) %>%
    summarise_all(funs(mean(., na.rm=TRUE))) %>%
    ungroup %>%
    dplyr::select(time, E, I, incid)

  fc <- full_join(scen_fc_vacpts, novac_fc_vacpts, by = c("time")) %>%
    dplyr::mutate(diff_I = I-I_sc, diff_incid = incid-incid_sc, flag_I = ifelse(I_sc > I, TRUE, FALSE), flag_incid = ifelse(incid_sc > incid, TRUE, FALSE))
  return(fc)

}


find_elimPeriod <- function(simdata, vacstart, pop){
  ## should work with sims_tm and sims_if
  newSimdata <- simdata %>% 
    filter(time >= vacstart) %>%
    dplyr::mutate(
      elimA = ifelse(new_cases < 1/10000*pop, TRUE, FALSE),
      elimB = ifelse(new_cases < 1/10000000*pop, TRUE, FALSE)) ## change new_cases<1
  uqsims <- unique(newSimdata$uqsim)

  returnData <- map_dfr(uqsims, function(sim){
    rleSimdata <- newSimdata %>% dplyr::filter(uqsim == sim)
    dummyA <- rep(F, nrow(rleSimdata))
    rle.resultsA <- rle(rleSimdata$elimA)
    if(sum(rle.resultsA$values & rle.resultsA$lengths > 52 & rle.resultsA$lengths == max(rle.resultsA$lengths))){ ## max elim duration that is at least over 52 weeks long
      
      pre.indexA = min(which(rle.resultsA$values & rle.resultsA$lengths > 52 & rle.resultsA$lengths == max(rle.resultsA$lengths)))
      post.indexA = min(which(rle.resultsA$values & rle.resultsA$lengths > 52 & rle.resultsA$lengths == max(rle.resultsA$lengths))) + max(rle.resultsA$lengths) - 1
      dummyA[pre.indexA:post.indexA] <- T
      print(paste("elimA in sim", sim))
    } 

    dummyB <- rep(F, nrow(rleSimdata))
    rle.resultsB <- rle(rleSimdata$elimB)
    if(sum(rle.resultsB$values & rle.resultsB$lengths > 52 & rle.resultsB$lengths == max(rle.resultsB$lengths))){ ## max elim duration that is at least over 52 weeks long
      
      pre.indexB = min(which(rle.resultsB$values & rle.resultsB$lengths > 52 & rle.resultsB$lengths == max(rle.resultsB$lengths)))
      post.indexB = min(which(rle.resultsB$values & rle.resultsB$lengths > 52 & rle.resultsB$lengths == max(rle.resultsB$lengths))) + max(rle.resultsB$lengths) - 1
      dummyB[pre.indexB:post.indexB] <- T
      print(paste("elimB in sim", sim))
    } 

    rleSimdata %>% dplyr::mutate(elimPeriodA = dummyA, elimPeriodB = dummyB)
  })

  return(returnData)
}

#### plotting utils ###########################################

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

#### plots ###########################################
plot_fit <- function(plt_df){
  fit_plt <- ggplot(plt_df, aes(x = week, y = cases_med, group = period)) +
    geom_ribbon(aes(ymin = cases_low, ymax = cases_up, fill = period), alpha = 0.5) +
    geom_line(aes(colour = period), size = 1) +
    geom_point(aes(y = data), colour = "black", alpha = 0.5, size = 1) +
    theme_bw() +
    theme(legend.position = "none", legend.title = element_blank()) + 
    scale_colour_manual(values = colours_two, aesthetics = c("fill", "colour")) +
    scale_x_date("week") +
    scale_y_continuous("cases", limit = c(0,30000)) +
    guides(fill = "none", colour = "none")
  return(fit_plt)
}

## main-forecast-...
plot_forecast <- function(plt_df){
  fc_plt <- ggplot(plt_df, aes(x = time, y = cases_med, group = period)) +
    geom_ribbon(aes(ymin = cases_low, ymax = cases_up, fill = period), alpha = 0.5) +
    geom_line(aes(colour = period), size = 1) +
    geom_point(aes(y = data), colour = "black", alpha = 0.5, size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    scale_colour_manual(values = colours_two, aesthetics = c("fill", "colour")) +
    scale_x_continuous("week") +
    scale_y_continuous("cases", limit = c(0,30000)) +
    guides(fill = "none")
  return(fc_plt)
}

plot_forecast_zm <- function(plt_df){
  fc_plt <- ggplot(plt_df, aes(x = time, y = cases_med, group = period)) +
    geom_ribbon(aes(ymin = cases_low, ymax = cases_up, fill = period), alpha = 0.5) +
    geom_line(aes(colour = period), size = 1) +
    geom_point(aes(y = data), colour = "black", alpha = 0.5, size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    scale_colour_manual(values = colours_two, aesthetics = c("fill", "colour")) +
    scale_x_continuous("week", limits = c(350,NA)) +
    scale_y_continuous("cases", limit = c(0,500)) +
    guides(fill = "none")
  return(fc_plt)
}
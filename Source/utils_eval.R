
## -------------------------------------- ##
## PLOTTING
breakinterval <- "1 year"

## observed fit
plot_of <- function(fdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(fdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = cases_lo, ymax = cases_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = cases_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = breakinterval) +
    labs(x = "Week", y = "Reported, Symptomatic Incidence") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

## observed forecast
plot_ofc <- function(fcdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(fcdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = cases_lo, ymax = cases_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = cases_med, colour = period), size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = breakinterval) +
    labs(x = "Week", y = "Reported, Symptomatic Incidence") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

## observed fit & forecast
plot_offc_full <- function(ffcdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(ffcdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = cases_lo, ymax = cases_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = cases_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = breakinterval) +
    labs(x = "Week", y = "Reported, Symptomatic Incidence") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

# true forecast
plot_tfc <- function(fcdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(fcdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = incid_lo, ymax = incid_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = incid_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = breakinterval) +
    labs(x = "Week", y = "True Cases") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

# true fit & forecast
plot_tffc_full <- function(ffcdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a", "forecast" = "#8d2300")
  plt <- ggplot(ffcdata, aes(x = weekdate)) +
    geom_ribbon(aes(ymin = incid_lo, ymax = incid_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = incid_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black", alpha = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size=11), axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = breakinterval) +
    labs(x = "Week", y = "True Cases") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}

## -------------------------------------- ##
## process data for elimination

process_fc_elim <- function(fcstate_fns){
  fcsettings <- fix_fc_settings()
  ## fcstate_fns include path
  ## sim, week, cases, S, E, I, A, R, N, incid, loglik
  map_dfr(1:length(fcstate_fns), function(i){
    par_id <- gsub(paste0("_seed", fcsettings$seed, ".csv"), "", gsub(paste0("haiti_mass_ocv/GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/fcstates_stoch_", fcsettings$str_if2, "_parid"), "", fcstate_fns[i]))
    print(paste("***** par id", par_id, "*****"))
    dummy <- read_csv(fcstate_fns[i]) %>%
      dplyr::select(sim, week, cases, incid, asymV, N)
    output <- find_elimPeriod(dummy) %>%
      dplyr::mutate(parid = as.numeric(par_id))
    rm(dummy)
    gc()
    return(output)
  }) 
}

find_elimPeriod <- function(simdata){
  vacstart <- min(simdata$week)+4
  ## should work with sims_tm and sims_if
  simdata2 <- simdata %>% 
    dplyr::filter(week >= vacstart) %>%
    dplyr::mutate(true_incid = incid + asymV) %>%
    ungroup 
  uqsims <- unique(simdata2$sim)

  returnData <- map_dfr(uqsims, function(uqsim){
    rleSimdata <- simdata2 %>% dplyr::filter(sim == uqsim) %>%
      dplyr::mutate(window_marker = seq_along(sim))
    print(paste("uqsim", uqsim, "**********"))

    ## 1/10000 incidence threshold
    rleinit <- rleSimdata %>% dplyr::filter(window_marker %in% 1:52)
    while(sum(rleinit$true_incid)>=(1/10000*mean(rleinit$N)) & nrow(rleinit)>=52){
      ## use the while loop to find and filter the minimum (52 weeks) elimination dataset
      new_window_start <- min(rleinit$window_marker)+1
      new_window_end <- new_window_start+51
      rleinit <- rleSimdata %>% dplyr::filter(window_marker %in% new_window_start:new_window_end)
    }
    if (nrow(rleinit) < 52){ ## no elimination period dataset
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- NA
      print(paste("no elim @ 1/10K threshold"))
    } else{ ## find true elimination period dataset
      elim_testresurg_start <- max(rleinit$window_marker)+1
      testresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_testresurg_start:(elim_testresurg_start+52))
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- FALSE
      print(paste("elim achieved @ 1/10K threshold"))
      while(sum(testresurg$true_incid)<(1/10000*mean(testresurg$N)) & max(testresurg$window_marker)<max(rleSimdata$window_marker)){
        new_window_start <- min(testresurg$window_marker)+1
        new_window_end <- new_window_start+51
        testresurg <- rleSimdata %>% dplyr::filter(window_marker %in% new_window_start:new_window_end)
      }
      if(sum(testresurg$true_incid)>=(1/10000*mean(testresurg$N))){ ## resurgence occurs
        ## find entire duration of elimination period (prior to resurgence)
        elim_beforeresurg_start <- min(rleinit$window_marker)
        elim_beforeresurg_end <- min(testresurg$window_marker)-1
        beforeresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_beforeresurg_start:elim_beforeresurg_end)
        dummyA[which(rleSimdata$window_marker %in% beforeresurg$window_marker)] <- TRUE
        resA <- TRUE
        print(paste("resurgence occurs"))
      } else{ ## resurgence does not occur
        elim_start <- min(rleinit$window_marker)
        elim_end <- max(testresurg$window_marker)
        noresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_start:elim_end) 
        dummyA[which(rleSimdata$window_marker %in% noresurg$window_marker)] <- TRUE
        resA <- FALSE
        print(paste("resurgence does not occur"))
      }
      
    }

    ## 1 case threshold
    print("processing one case threshold")
    dummyB <- rep(F, nrow(rleSimdata))
    resB <- NA
    rle.resultsB <- rle(rleSimdata$true_incid)
    if(rle.resultsB$values[length(rle.resultsB$values)]==0 & rle.resultsB$lengths[length(rle.resultsB$values)] > 52){ 
      ## elimination period that is at least 52 weeks long and endures through the end of the 10 year simulation
      first.rle.indexB = length(rle.resultsB$values)
      pre.indexB = ifelse(first.rle.indexB>1, sum(rle.resultsB$lengths[1:(first.rle.indexB-1)])+1, 1)
      post.indexB = pre.indexB + rle.resultsB$lengths[first.rle.indexB]-1
      dummyB[pre.indexB:post.indexB] <- T
      resB <- NA ## 6/11/2019 resurgence is not possible in this new definition of elimination
      print(paste("elim @ 1 case threshold in sim", uqsim))
    } else{
      print(paste("no elim @ 1 case threshold"))
    }
    
    rleSimdata %>% 
      dplyr::mutate(elimPeriodA = dummyA, elimPeriodB = dummyB, 
                    resurgA = resA, 
                    resurgB = resB) %>% ## is there resurgence to non-elimination levels for a given simulation (same value for all time steps)
      dplyr::select(sim, week, true_incid, N, elimPeriodA, elimPeriodB, resurgA, resurgB)
  })

  return(returnData)
}

find_probElim <- function(elimdata){
  year3 <- round(min(elimdata$week) + 3*52.14)
  year5 <- round(min(elimdata$week) + 5*52.14)
  year10 <- round(min(elimdata$week) + 10*52.14)
  
  y3 <- elimdata %>%
    dplyr::filter(week == year3) %>%
    summarise(nElimA = sum(elimPeriodA), nElimB = sum(elimPeriodB), nsims = nrow(.), resurgA = sum(resurgA, na.rm=TRUE), resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, pr_elim_thresh2 = nElimB/nsims, pr_resurg_thresh1 = resurgA/nElimA, pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 3)
  y5 <- elimdata %>%
    dplyr::filter(week == year5) %>%
    summarise(nElimA = sum(elimPeriodA), nElimB = sum(elimPeriodB), nsims = nrow(.), resurgA = sum(resurgA, na.rm=TRUE), resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, pr_elim_thresh2 = nElimB/nsims, pr_resurg_thresh1 = resurgA/nElimA, pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 5) 
  y10 <- elimdata %>%
    dplyr::filter(week == year10) %>%
    summarise(nElimA = sum(elimPeriodA), nElimB = sum(elimPeriodB), nsims = nrow(.), resurgA = sum(resurgA, na.rm=TRUE), resurgB = sum(resurgB, na.rm=TRUE)) %>%
    dplyr::mutate(pr_elim_thresh1 = nElimA/nsims, pr_elim_thresh2 = nElimB/nsims, pr_resurg_thresh1 = resurgA/nElimA, pr_resurg_thresh2 = resurgB/nElimB) %>%
    dplyr::mutate(time_frame = 10) 
  
  bind_rows(y3, y5, y10)
}

find_timeToElim <- function(elimdata){
  vacstart <- min(elimdata$week)
  periodA <- elimdata %>%
    dplyr::filter(elimPeriodA) %>%
    group_by(parid, sim) %>% 
    dplyr::filter(week == min(week)) %>%
    ungroup %>%
    dplyr::mutate(tElim_wks = week-vacstart) %>%
    summarise(t_elim_thresh1_med = median(tElim_wks), t_elim_thresh1_low = quantile(tElim_wks, probs = c(.025)), t_elim_thresh1_hi = quantile(tElim_wks, probs = c(.975)))

  periodB <- elimdata %>%
    dplyr::filter(elimPeriodB) %>%
    group_by(parid, sim) %>% 
    dplyr::filter(week == min(week)) %>%
    ungroup %>%
    dplyr::mutate(tElim_wks = week-vacstart) %>%
    summarise(t_elim_thresh2_med = median(tElim_wks), t_elim_thresh2_low = quantile(tElim_wks, probs = c(.025)), t_elim_thresh2_hi = quantile(tElim_wks, probs = c(.975))) 

  bind_cols(periodA, periodB)
}


process_fc_elim_novac <- function(fcstate_fns){
  fcsettings <- fix_fc_settings()
  ## fcstate_fns include path
  ## sim, week, cases, S, E, I, A, R, N, incid, loglik
  map_dfr(1:length(fcstate_fns), function(i){
    par_id <- gsub(paste0("_seed", fcsettings$seed, ".csv"), "", gsub(paste0("haiti_mass_ocv/GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/fcstates_stoch_", fcsettings$str_if2, "_parid"), "", fcstate_fns[i]))
    print(paste("***** par id", par_id, "*****"))
    dummy <- read_csv(fcstate_fns[i]) %>%
      dplyr::select(sim, week, cases, incid, N)
    output <- find_elimPeriod_novac(dummy) %>%
      dplyr::mutate(parid = as.numeric(par_id))
    rm(dummy)
    gc()
    return(output)
  }) 
}

find_elimPeriod_novac <- function(simdata){
  vacstart <- min(simdata$week)+4
  ## should work with sims_tm and sims_if
  simdata2 <- simdata %>% 
    dplyr::filter(week >= vacstart) %>%
    dplyr::mutate(true_incid = incid) %>%
    ungroup
  uqsims <- unique(simdata2$sim)

  returnData <- map_dfr(uqsims, function(uqsim){
    rleSimdata <- simdata2 %>% dplyr::filter(sim == uqsim) %>%
      dplyr::mutate(window_marker = seq_along(sim))
    print(paste("uqsim", uqsim, "****************"))

    ## 1/10000 incidence threshold
    rleinit <- rleSimdata %>% dplyr::filter(window_marker %in% 1:52)
    while(sum(rleinit$true_incid)>=(1/10000*mean(rleinit$N)) & nrow(rleinit)>=52){
      ## use the while loop to find and filter the minimum (52 weeks) elimination dataset
      new_window_start <- min(rleinit$window_marker)+1
      new_window_end <- new_window_start+51
      rleinit <- rleSimdata %>% dplyr::filter(window_marker %in% new_window_start:new_window_end)
    }
    if (nrow(rleinit) < 52){ ## no elimination period dataset
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- NA
      print(paste("no elim @ 1/10K threshold"))
    } else{ ## find true elimination period dataset
      # browser()
      elim_testresurg_start <- max(rleinit$window_marker)+1
      testresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_testresurg_start:(elim_testresurg_start+52))
      dummyA <- rep(F, nrow(rleSimdata))
      resA <- FALSE
      print(paste("elim @ 1/10K threshold"))
      while(sum(testresurg$true_incid)<(1/10000*mean(testresurg$N)) & max(testresurg$window_marker)<max(rleSimdata$window_marker)){
        new_window_start <- min(testresurg$window_marker)+1
        new_window_end <- new_window_start+51
        testresurg <- rleSimdata %>% dplyr::filter(window_marker %in% new_window_start:new_window_end)
      }
      if(sum(testresurg$true_incid)>=(1/10000*mean(testresurg$N))){ ## resurgence occurs
        print(paste(sum(testresurg$true_incid), (1/10000*mean(testresurg$N))))
        ## find entire duration of elimination period (prior to resurgence)
        elim_beforeresurg_start <- min(rleinit$window_marker)
        elim_beforeresurg_end <- min(testresurg$window_marker)-1
        beforeresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_beforeresurg_start:elim_beforeresurg_end)
        dummyA[which(rleSimdata$window_marker %in% beforeresurg$window_marker)] <- TRUE
        resA <- TRUE
        print(paste("resurgence occurs"))
      } else{ ## resurgence does not occur
        print(paste(sum(testresurg$true_incid), (1/10000*mean(testresurg$N))))
        elim_start <- min(rleinit$window_marker)
        elim_end <- max(testresurg$window_marker)
        noresurg <- rleSimdata %>% dplyr::filter(window_marker %in% elim_start:elim_end) 
        dummyA[which(rleSimdata$window_marker %in% noresurg$window_marker)] <- TRUE
        resA <- FALSE
        print(paste("resurgence does not occur"))
      }
      
    }

    ## 1 case threshold
    print("processing one case threshold")
    dummyB <- rep(F, nrow(rleSimdata))
    resB <- NA
    rle.resultsB <- rle(rleSimdata$true_incid)
    if(rle.resultsB$values[length(rle.resultsB$values)]==0 & rle.resultsB$lengths[length(rle.resultsB$values)] > 52){ 
      ## elimination period that is at least 52 weeks long and endures through the end of the 10 year simulation
      first.rle.indexB = length(rle.resultsB$values)
      pre.indexB = ifelse(first.rle.indexB>1, sum(rle.resultsB$lengths[1:(first.rle.indexB-1)])+1, 1)
      post.indexB = pre.indexB + rle.resultsB$lengths[first.rle.indexB]-1
      dummyB[pre.indexB:post.indexB] <- T
      resB <- NA ## 6/11/2019 resurgence is not possible in this new definition of elimination
      print(paste("elim @ 1 case threshold in sim", uqsim))
    } else{
      print(paste("no elim @ 1 case threshold"))
    }
    
    rleSimdata %>% 
      dplyr::mutate(elimPeriodA = dummyA, elimPeriodB = dummyB, 
                    resurgA = resA, 
                    resurgB = resB) %>% ## is there resurgence to non-elimination levels for a given simulation (same value for all time steps)
      dplyr::select(sim, week, true_incid, N, elimPeriodA, elimPeriodB, resurgA, resurgB)
  })

  return(returnData)
}

## -------------------------------------- ##
## process data for plotting

process_fepi_plot <- function(fitstate_fns, hti_data){
  ## fitstate_fns include path
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i]) 
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::rename(cmed = cases_med, imed = incid_med) %>%
    dplyr::summarise(cases_med = median(cmed), cases_lo = quantile(cmed, probs = c(.025)), cases_hi = quantile(cmed, probs = c(.975)), incid_med = median(imed), incid_lo = quantile(imed, probs = c(.025)), incid_hi = quantile(imed, probs = c(.975)), N_med = median(N_med)) %>%
    dplyr::full_join(hti_data %>% dplyr::rename(orig_cases = cases), by = c("week"))
}

process_fend_plot <- function(fitstate_fns, hti_data, cname){
  episettings <- fix_epi_settings()
  ## fitstate_fns include path
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i]) %>% 
      dplyr::mutate(week = week + episettings$nweeks) 
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::rename(cmed = cases_med, imed = incid_med) %>%
    dplyr::summarise(cases_med = median(cmed), cases_lo = quantile(cmed, probs = c(.025)), cases_hi = quantile(cmed, probs = c(.975)), incid_med = median(imed), incid_lo = quantile(imed, probs = c(.025)), incid_hi = quantile(imed, probs = c(.975)), N_med = median(N_med)) %>%
    dplyr::full_join(hti_data %>% dplyr::select(-week_end) %>% dplyr::rename(orig_cases = cases), by = c("week")) %>%
    dplyr::arrange(week)
}

process_ffc_plot <- function(fcstate_fns){
  ## fcstate_fns include path
  ## sim, week, cases, S, E, I, A, R, N, incid, loglik
  map_dfr(1:length(fcstate_fns), function(i){
    readRDS(fcstate_fns[i])
    # read_csv(fcstate_fns[i])
  }) %>%
    dplyr::mutate(orig_cases = NA) %>%
    dplyr::select(week, cases_med, cases_lo, cases_hi, incid_med, incid_lo, incid_hi, orig_cases, N_med)
}


## -------------------------------------- ##
## process cumulative vaccination plot

get_vac_settings <- function(scencode, ndata){
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

  vactab <- make.vactab(t0 = 0, 
                      tmax = ndata+53*6, 
                      ndept = fc_nd,
                      nweeks = fc_nw,
                      coverage_2dose = fc_c2, 
                      coverage_1dose = fc_c1,
                      first_vac_t = ndata+4, 
                      ve_scen = fc_vescen) %>%
    dplyr::select(time, num_vacc) %>% 
    dplyr::filter(time > ndata)

  return(vactab)
}

get_date_df <- function(){
  ndp <- nrow(get.mspp.agg.data())

  date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
  datedf <- data.frame(saturday_date = seq(date.start$date_sat, length.out = ndp+52*11, by = 7), week = seq_len(ndp+52*11))
  return(datedf)
}

process_cumvac_output <- function(){
  ndp <- nrow(get.mspp.agg.data())

  date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
  datedf <- data.frame(weekdate = seq(date.start$date_sat, length.out = ndp+52*11, by = 7), week = seq_len(ndp+52*11))

  s1v <- get_vac_settings("id1", ndp) %>% dplyr::mutate(scenario = 1)
  s2v <- get_vac_settings("id2", ndp) %>% dplyr::mutate(scenario = 2)
  s3v <- get_vac_settings("id3", ndp) %>% dplyr::mutate(scenario = 3)
  s4v <- get_vac_settings("id4", ndp) %>% dplyr::mutate(scenario = 4)

  vacc_scens <- bind_rows(s1v, s2v, s3v, s4v) %>%
    rename(week = time) %>%
    dplyr::mutate(scenario = as.character(scenario)) %>%
    left_join(datedf, by = c("week")) %>%
    dplyr::rename(vacc_med = num_vacc)

  return(vacc_scens)
}

process_cumvac_plot <- function(datedf){
  ndp <- nrow(get.mspp.agg.data())
  s1v <- get_vac_settings("id1", ndp) %>% dplyr::rename(s1 = num_vacc)
  s2v <- get_vac_settings("id2", ndp) %>% dplyr::rename(s2 = num_vacc)
  s3v <- get_vac_settings("id3", ndp) %>% dplyr::rename(s3 = num_vacc)
  s4v <- get_vac_settings("id4", ndp) %>% dplyr::rename(s4 = num_vacc)

  vacc_scens <- full_join(s1v, s2v, by = c("time")) %>%
    full_join(s3v, by = c("time")) %>%
    full_join(s4v, by = c("time")) %>%
    rename(week = time) %>%
    left_join(datedf, by = c("week")) 

  return(vacc_scens)
}

plot_cumvac <- function(vacc_scens){

  lastvactime <- vacc_scens %>% dplyr::mutate(any_vacc = s1 + s2 + s3 + s4) %>%
    dplyr::filter(any_vacc > 0) %>%
    dplyr::filter(week == max(week)) %>% dplyr::select(week) %>% unlist %>% unname

  deploylabs <- data.frame(scenario = c("s2", "s4", "s3", "s1"), pltlab = c("2 dept", "3 dept", "slow national", "fast national"))

  pltdat <- tbl_df(vacc_scens) %>%
    dplyr::filter(week <= lastvactime+1) %>% 
    dplyr::mutate_at(vars(c("s1", "s2", "s3", "s4")), cumsum) %>%
    dplyr::select(weekdate, week, contains("s")) %>%
    gather(scenario, num_vacc, s1:s4) %>%
    dplyr::mutate(scenario = factor(scenario, levels = deploylabs$scenario, labels = deploylabs$pltlab))


  plt <- ggplot(pltdat, aes(x = weekdate, y = num_vacc, group = scenario)) +
    geom_line(aes(colour = scenario), size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_viridis_d("deployment strategy", option = "C") +
    scale_y_continuous("Vaccinated People") +
    scale_x_date(date_breaks = "6 months")

  return(plt)

}

plot_cumvac_frac <- function(vacc_scens, plotdf_fc){

  lastvactime <- vacc_scens %>% dplyr::mutate(any_vacc = s1 + s2 + s3 + s4) %>%
    dplyr::filter(any_vacc > 0) %>%
    dplyr::filter(week == max(week)) %>% dplyr::select(week) %>% unlist %>% unname

  deploylabs <- data.frame(scenario = c("s2", "s4", "s3", "s1"), pltlab = c("2 dept", "3 dept", "slow national", "fast national"))

  pop_df <- plotdf_fc %>% dplyr::select(week, N_med)

  pltdat <- tbl_df(vacc_scens) %>%
    dplyr::filter(week <= lastvactime+1) %>% 
    dplyr::mutate_at(vars(c("s1", "s2", "s3", "s4")), cumsum) %>%
    dplyr::select(weekdate, week, contains("s")) %>%
    gather(scenario, num_vacc, s1:s4) %>%
    left_join(pop_df, by = c("week")) %>%
    dplyr::mutate(scenario = factor(scenario, levels = deploylabs$scenario, labels = deploylabs$pltlab),
                  frac_vacc = num_vacc/N_med)


  plt <- ggplot(pltdat, aes(x = weekdate, y = frac_vacc, group = scenario)) +
    geom_line(aes(colour = scenario), size = 1) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + 
    scale_colour_viridis_d("deployment strategy", option = "C") +
    scale_y_continuous("Fraction of Population Vaccinated", limits = c(0,1)) +
    scale_x_date(date_breaks = "6 months")

  return(plt)

}
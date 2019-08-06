
reload_source <- function(){
  options(warn=1,error=traceback)
  library(plyr)
  library(tidyverse)
  library(pomp)
  library(devtools)
  library(data.table)
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(iterators)
  library(pkgbuild)
  library(gridExtra)
  library(GGally)

  # checkpoint("2019-04-29", project = "haiti_ocv_JHUmodel/") ## original R version was 3.4.4
  find_rtools()
  has_devel()
  source("haiti_ocv_JHUmodel/Source/utils_clean.R")
  # source("haiti_ocv_JHUmodel/Source/src_generate_epi_params.R")

}

## -------------------------------------- ##
## get haiti population
get.haiti.pop <- function(){
  pop <- sum(c(746236, 1727524, 4029705, 728807, 1067177, 774976, 342525, 393967, 632601, 468301)) ## use sum of population in departments 
  return(pop)
}

## -------------------------------------- ##
## load mspp data aggregated to country scale
get.mspp.agg.data <- function(){

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
    dplyr::mutate(date_sat = as.Date(paste(year, month, day, sep = "-"), origin = "1900-01-01"))
  
  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by=7))
    
  cleanDat <- allDat %>% 
    dplyr::mutate(date_sat = as.Date(dateDf$date_sat, origin = "1900-01-01")) %>%
    dplyr::select(-date_sat_orig) %>%
    full_join(fullDateVec, by = c("date_sat")) %>%
    dplyr::arrange(date_sat) %>%
    dplyr::mutate(week = seq_along(date_sat)) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)

  aggDat <- cleanDat %>%
    dplyr::group_by(week) %>% ## "day" or "week"
    dplyr::summarise(cases = sum(cases, na.rm=TRUE))

  return(aggDat)
} 

## -------------------------------------- ##
## load mspp data aggregated to department scale
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
    dplyr::mutate(date_sat = as.Date(paste(year, month, day, sep = "-"), origin = "1900-01-01"))
  
  fullDateVec <- data.frame(date_sat = seq(min(dateDf$date_sat), max(dateDf$date_sat), by=7))
    
  cleanDat <- allDat %>% 
    dplyr::mutate(date_sat = as.Date(dateDf$date_sat, origin = "1900-01-01")) %>%
    dplyr::select(-date_sat_orig) %>%
    full_join(fullDateVec, by = c("date_sat")) %>%
    dplyr::arrange(date_sat) %>%
    dplyr::mutate(week = seq_along(date_sat)) %>%
    tidyr::gather(department, cases, Artibonite:Sud_Est)

  return(cleanDat)
} 

## -------------------------------------- ##
## construct seasonal b-spline terms
make.covartab <- function(t0, tmax, byt=1, nbasis=6, degree=6, per=52.14){

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

## -------------------------------------- ##
## construct vaccination covariate table
make.vactab <- function(t0=0, tmax, byt=1, ndept=10, nweeks=10, coverage_2dose, coverage_1dose, first_vac_t, ve_scen){

  ## deployment scenario 1 (fast national): one departmental campaign every 10 weeks for all depts

  tbasis <- seq(from=t0,to=tmax,by=byt)
  time_check <- c()
  for(i in 1:ndept) {
    time_check <- c(time_check, rep(0, nweeks-1), rep(i, 1))
  } ## 10 week vacc campaigns (nweeks)

  ## number of vaccines per week by department
  pop_dept <- data.frame(ocv_order = 1:10, dept = c("Centre", "Artibonite", "Ouest", "Nord Ouest", "Nord", "Sud", "Nippes", "Nord Est", "Sud Est", "Grand'Anse"), pop = c(746236, 1727524, 4029705, 728807, 1067177, 774976, 342525, 393967, 632601, 468301)) %>%
    dplyr::mutate(num_vacc = (coverage_2dose+coverage_1dose)*pop/1) ## pulse vaccinees in last 1 week of campaign
  # browser()
  ## create dataframe with number vaccinated for each campaign
  vactab <- data.frame(time = tbasis, vac_tcheck = 0)
  vactab[which(vactab$time %in% first_vac_t:(first_vac_t+length(time_check)-1)),]$vac_tcheck <- time_check
  vactab2 <- left_join(vactab, pop_dept %>% dplyr::select(ocv_order, num_vacc), by = c("vac_tcheck"="ocv_order")) %>%
    dplyr::mutate(num_vacc = ifelse(is.na(num_vacc), 0, round(num_vacc))) 

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
    dplyr::mutate_at(vars(contains("ve_")), list(~ifelse(is.na(.), 0, .)))

  return(vactab3)
}
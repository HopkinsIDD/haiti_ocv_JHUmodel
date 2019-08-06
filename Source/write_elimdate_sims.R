## write elimination 

source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()

outdir <- "haiti_ocv_JHUmodel/GeneratedData/final1_summaries/"

#### read Haiti data ####
deptdat <- get.mspp.dept.data() 
fittimes <- deptdat %>%
  distinct(week, date_sat) 
htdat <- get.mspp.agg.data() %>%
  left_join(fittimes, by = c("week"))

fcweeks <- (max(fittimes$week)+1):(max(fittimes$week)+521)
fcdates <- seq.Date(max(fittimes$date_sat)+7, along.with = fcweeks, by = "week")
fctimes <- tbl_df(data.frame(week = fcweeks, date_sat = fcdates))

#### read JHU data ####
s1 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_id1/fcelim_if_id1_nsims25_seed20192104.rds") %>% 
  dplyr::mutate(scenario = 1)
s2 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_id2/fcelim_if_id2_nsims25_seed20192104.rds") %>% 
  dplyr::mutate(scenario = 2)
s3 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_id3/fcelim_if_id4_nsims25_seed20192104.rds") %>% 
  dplyr::mutate(scenario = 3)
s4 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_id4/fcelim_if_id4_nsims25_seed20192104.rds") %>% 
  dplyr::mutate(scenario = 4)
s25 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_id25/fcelim_if_id25_nsims25_seed20192104.rds") %>% 
  dplyr::mutate(scenario = 25)
s0 <- readRDS("haiti_ocv_JHUmodel/GeneratedData/final1_novac/fcelim_if_novac_nsims25_seed20192104.rds") %>%
  dplyr::mutate(scenario = 0)

process_elimweek <- function(df){
  left_join(
    df %>% distinct(scenario, sim, parid), 
    df %>% 
      group_by(sim, parid) %>% 
      dplyr::filter(elimPeriodB) %>%
      dplyr::filter(week == min(week)) %>%
      dplyr::select(scenario, sim, parid, week), 
    by = c("scenario","sim","parid")
  )
}

output <- bind_rows(process_elimweek(s0), process_elimweek(s1), process_elimweek(s2), process_elimweek(s3), process_elimweek(s4), process_elimweek(s25)) %>%
  left_join(fctimes, by = c("week")) %>%
  dplyr::rename(elimination_date = date_sat) %>%
  dplyr::mutate(team = "JHU") 

plt <- ggplot(output, aes(x=week)) + 
  geom_bar(aes(fill=factor(scenario)), position="dodge") +
  theme(legend.position = "bottom")

print(plt)

write_csv(output %>% dplyr::select(team, scenario, elimination_date), paste0(outdir, "elimdate_JHU_2019-06-26.csv"))
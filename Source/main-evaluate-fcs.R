source("haiti_ocv_JHUmodel/Source/utils_clean.R")
source("haiti_ocv_JHUmodel/Source/utils_eval.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_fc_params.R")


## -------------------------------------- ##
## script settings ##
process_modules <- TRUE; plot_modules <- FALSE; cumvac_modules <- TRUE

## -------------------------------------- ##
## model settings ##
episettings <- fix_epi_settings()
endsettings <- fix_end_settings()
fcsettings <- fix_fc_settings()
fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", fcsettings$mcode, "_", fcsettings$scode, "/")
gen_fig_dir <- paste0("haiti_ocv_JHUmodel/figures/supp_052019/")
out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
out_dir_fc <- paste0("haiti_ocv_JHUmodel/GeneratedData/", fcsettings$mcode, "_", fcsettings$scode, "/")
dir.create(gen_fig_dir, showWarnings = FALSE)
dir.create(out_dir_fc, showWarnings = FALSE)

## -------------------------------------- ##
## data import and setup ##
pop.haiti <- get.haiti.pop()

haiti.dat <- get.mspp.agg.data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks) 
haiti.dat.end <- haiti.dat %>% dplyr::filter(week >= episettings$nweeks) %>%
  dplyr::mutate(week_end = seq_along(week)-1)

covartab.all <- make.covartab(0, nrow(haiti.dat)+1, byt=1, degree=6, nbasis=6, per=52.14)
covartab <- covartab.all[(episettings$nweeks+1):nrow(covartab.all),] %>%
  dplyr::mutate(time_end = seq_along(time)-1)
num_betas <- ncol(covartab.all)-1


fn_epi <- list.files(out_dir_epi, "fitstates_if_")
fn_epi_full <- paste0(out_dir_epi, fn_epi)
fn_end <- list.files(out_dir_end, "fitstates_if_")
fn_end_full <- paste0(out_dir_end, fn_end)
fn_fc_elim <- list.files(out_dir_fc, "fcstates_stoch_")
fn_fc_elim_full <- paste0(out_dir_fc, fn_fc_elim)
fn_fc_plt <- list.files(out_dir_fc, "stochfcs_") ## summary across all paramsets & sims
fn_fc_plt_full <- paste0(out_dir_fc, fn_fc_plt)

## -------------------------------------- ##
## data preparation -- plotting ##
if(process_modules){
  fit_epidat <- process_fepi_plot(fn_epi_full, haiti.dat.epi) %>%
    dplyr::mutate(period = "epidemic")
  fit_enddat <- process_fend_plot(fn_end_full, haiti.dat.end) %>% 
    dplyr::filter(week > episettings$nweeks) %>%
    dplyr::mutate(period = "endemic")
  fitdat <- bind_rows(fit_epidat, fit_enddat)

  fcdat <- process_ffc_plot(fn_fc_plt_full) %>%
    dplyr::filter(week > max(haiti.dat$week)) %>%
    dplyr::mutate(period = "forecast")
  ffcdat <- bind_rows(fitdat, fcdat)

  date.start <- get.mspp.dept.data() %>% dplyr::filter(date_sat == min(date_sat)) %>% dplyr::distinct(date_sat)
  datedat <- data.frame(weekdate = seq(date.start$date_sat, length.out = nrow(ffcdat), by = 7), week = seq_len(nrow(ffcdat)))

  plotdat_ffc <- full_join(ffcdat, datedat, by = c("week"))
  plotdat_f <- plotdat_ffc %>% dplyr::filter(week <= max(haiti.dat$week))
  plotdat_fc <- plotdat_ffc %>% dplyr::filter(week > max(haiti.dat$week))

}

## -------------------------------------- ##
## plot clean fits and forecasts ##

if(plot_modules){
  plt_of <- plot_of(plotdat_f)
  ggsave(paste0(fig_dir, "o_ffull_if_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_of, width = 6, height = 4, units = "in")
  print(plt_of)

  plt_ofc <- plot_ofc(plotdat_fc)
  ggsave(paste0(fig_dir, "o_fcfull_if_", fcsettings$str_if, "_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_ofc, width = 6, height = 4, units = "in")
  print(plt_ofc)

  plt_offc <- plot_offc_full(plotdat_ffc)
  ggsave(paste0(fig_dir, "o_ffcfull_if_", fcsettings$str_if, "_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_offc, width = 6, height = 4, units = "in")
  print(plt_offc
)
  plt_tfc <- plot_tfc(plotdat_fc)
  ggsave(paste0(fig_dir, "t_fcfull_if_", fcsettings$str_if, "_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_tfc, width = 6, height = 4, units = "in")
  print(plt_tfc)

  plt_tffc <- plot_tffc_full(plotdat_ffc)
  ggsave(paste0(fig_dir, "t_ffcfull_if_", fcsettings$str_if, "_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_tffc, width = 6, height = 4, units = "in")
  print(plt_tffc)

}

## -------------------------------------- ##
## plot cumulative vaccination across scenarios ##

if(cumvac_modules){
  cumvac_df <- process_cumvac_plot(datedat)
  plt_cumvac <- plot_cumvac(cumvac_df)
  ggsave(paste0(gen_fig_dir, "cumvacc_deployment.pdf"), plt_cumvac, width = 6, height = 4, units = "in")
  plt_cumvac_frac <- plot_cumvac_frac(cumvac_df, plotdat_fc)
  ggsave(paste0(gen_fig_dir, "cumvacc_deployment_frac.pdf"), plt_cumvac_frac, width = 6, height = 4, units = "in")
  print(plt_cumvac)
  print(plt_cumvac_frac)
}

## -------------------------------------- ##
## data preparation -- elimination ##
if(!file.exists(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))){
  fcelim <- process_fc_elim(fn_fc_elim_full)
  saveRDS(fcelim, paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
  } else{
    fcelim <- readRDS(paste0(out_dir_fc, "fcelim_if_", fcsettings$str_if, ".rds"))
  }

## -------------------------------------- ##
## print elimination ##

pElim_dat <- find_probElim(fcelim)
print(paste("probability of elimination", fcsettings$str_if))
print(pElim_dat)

tElim_dat <- find_timeToElim(fcelim)
print(paste("time to elimination (weeks)"))
print(tElim_dat)

sink(paste0(out_dir_fc, "fcelimprob.txt"))
print.table(pElim_dat)
print(paste("*********"))
print.table(tElim_dat)
sink()
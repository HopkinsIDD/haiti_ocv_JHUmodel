## rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_end_params.R")
source("haiti_ocv_JHUmodel/Source/pomp_fit_end.R")

## -------------------------------------- ##
## model settings ##
episettings <- fix_epi_settings()
endsettings <- fix_end_settings()
fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", endsettings$mcode, "/")
out_dir_epi <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
out_dir_end <- paste0("haiti_ocv_JHUmodel/GeneratedData/", endsettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir_end, showWarnings = FALSE)

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

starts <- generate_end_params(out_dir_epi, out_dir_end)

## -------------------------------------- ##
## custom functions ##
plot_fit <- function(fitdata){
  plt <- ggplot(fitdata, aes(x = week)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_hi), fill = "#458b74", alpha = 0.5) +
    geom_line(aes(y = est_med), size = 1, colour = "#458b74") +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw()
  return(plt)
}
plot_fit_full <- function(fitdata){
  colvec <- c("epidemic" = "#008d23", "endemic" = "#8d006a")
  plt <- ggplot(fitdata, aes(x = week)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_hi, fill = period), alpha = 0.5) +
    geom_line(aes(y = est_med, colour = period), size = 1) +
    geom_point(aes(y = orig_cases), colour = "black") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_fill_manual("", values = colvec, aesthetics = c("colour", "fill"))
  return(plt)
}
process_fitstates <- function(fitstate_fns, hti_data, cname){
  ## fitstate_fns include path
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i]) %>% 
      dplyr::mutate(week = week + episettings$nweeks) %>%
      dplyr::rename(est_cases = cname)
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::summarise(est_med = median(est_cases), est_lo = quantile(est_cases, probs = c(.025)), est_hi = quantile(est_cases, probs = c(.975))) %>%
    dplyr::full_join(hti_data %>% dplyr::select(-week_end) %>% dplyr::rename(orig_cases = cases), by = c("week")) %>%
    dplyr::arrange(week)
}
process_fitstates_epi <- function(fitstate_fns, hti_data){
  ## fitstate_fns include path
  map_dfr(1:length(fitstate_fns), function(i){
    read_csv(fitstate_fns[i]) 
  }) %>%
    ungroup %>%
    group_by(week) %>%
    dplyr::summarise(est_med = median(cases_med), est_lo = quantile(cases_med, probs = c(.025)), est_hi = quantile(cases_med, probs = c(.975))) %>%
    dplyr::full_join(hti_data %>% dplyr::rename(orig_cases = cases), by = c("week"))
}

## -------------------------------------- ##
## import endemic fits from all tms
fn_tm <- list.files(out_dir_end, "fitstates_tm_")
fn_tm_full <- paste0(out_dir_end, fn_tm)
fit_tm_summ <- process_fitstates(fn_tm_full, haiti.dat.end, "est_cases")

plt_tm_summ <- plot_fit(fit_tm_summ)
ggsave(paste0(fig_dir, "fit_tm_", endsettings$str_tm, ".pdf"), plt_tm_summ, width = 4, height = 4, units = "in")
print(plt_tm_summ)

plt_tm_summ_zm <- plot_fit(fit_tm_summ) +
  scale_x_continuous(limits = c(max(fit_tm_summ$week)-40, max(fit_tm_summ$week))) +
  scale_y_continuous(limits = c(0, 4000))
ggsave(paste0(fig_dir, "fitzm_tm_", endsettings$str_tm, ".pdf"), plt_tm_summ_zm, width = 4, height = 4, units = "in")
print(plt_tm_summ_zm)

## -------------------------------------- ##
## import endemic fits only from all converged tm
prof_tm <- readRDS(paste0(out_dir_end, "tmfits_", endsettings$str_tm, ".rds")) %>%
  arrange(desc(loglik))
parid_ls <- prof_tm %>% dplyr::filter(conv==0) %>% dplyr::select(parid) %>% unlist %>% unname
fn_tm_conv <- paste0(out_dir_end, "fitstates_tm_", endsettings$str_tm2, "_parid", parid_ls, "_seed", endsettings$seed, ".csv")
fit_tm_lim <- process_fitstates(fn_tm_conv, haiti.dat.end, "est_cases")

plt_tm_lim <- plot_fit(fit_tm_lim)
ggsave(paste0(fig_dir, "fitconv_tm_", endsettings$str_tm, ".pdf"), plt_tm_lim, width = 4, height = 4, units = "in")
print(plt_tm_lim)

## -------------------------------------- ##
## import endemic fits from all ifs
fn_if <- list.files(out_dir_end, "fitstates_if_")
fn_if_full <- paste0(out_dir_end, fn_if)
fit_if_summ <- process_fitstates(fn_if_full, haiti.dat.end, "cases_med")

plt_if_summ <- plot_fit(fit_if_summ)
ggsave(paste0(fig_dir, "fit_if_", endsettings$str_if, ".pdf"), plt_if_summ, width = 4, height = 4, units = "in")
print(plt_if_summ)

plt_if_summ_zm <- plot_fit(fit_if_summ) +
  scale_x_continuous(limits = c(max(fit_if_summ$week)-40,max(fit_if_summ$week))) +
  scale_y_continuous(limits = c(0, 2500))
ggsave(paste0(fig_dir, "fitzm_if_", endsettings$str_if, ".pdf"), plt_if_summ_zm, width = 4, height = 4, units = "in")
print(plt_if_summ_zm)

## -------------------------------------- ##
## import endemic fits from all converged ifs
prof_if <- readRDS(paste0(out_dir_end, "iffits_", endsettings$str_if, ".rds"))

parid_conv_ls <- prof_if %>% filter(nfail.max == 0) %>% select(parid) %>% unlist %>% unname
fn_if_conv <- paste0(out_dir_end, "fitstates_if_", endsettings$str_if2, "_parid", parid_conv_ls, "_seed", endsettings$seed, ".csv")
fit_if_lim <- process_fitstates(fn_if_conv, haiti.dat.end, "cases_med")

plt_if_lim <- plot_fit(fit_if_lim)
ggsave(paste0(fig_dir, "fitconv_if_", endsettings$str_if, ".pdf"), plt_if_lim, width = 4, height = 4, units = "in")
print(plt_if_lim) ## at least one filtering failure

## -------------------------------------- ##
## explore endemic estimates and likelihood
pltcols <-  c("rho", "theta", "beta1", "nu", "loglik")

all_if <- prof_if %>% dplyr::select(pltcols)
plt_all_if <- ggpairs(all_if)
ggsave(paste0(fig_dir, "est_coverage_", endsettings$str_if, ".pdf"), plt_all_if, width = 6, height = 6, units = "in")
print(plt_all_if)

lim_if <- prof_if %>% filter(loglik>-3000) %>% dplyr::select(pltcols)
plt_lim_if <- ggpairs(lim_if)
ggsave(paste0(fig_dir, "estlim_coverage_", endsettings$str_if, ".pdf"), plt_lim_if, width = 6, height = 6, units = "in")
print(plt_lim_if)

## -------------------------------------- ##
## explore full epidemic/endemic iterated filtering fits
fn_epi <- list.files(out_dir_epi, "fitstates_if_")
fn_epi_full <- paste0(out_dir_epi, fn_epi)
fn_end <- list.files(out_dir_end, "fitstates_if_")
fn_end_full <- paste0(out_dir_end, fn_if)

fit_epidat <- process_fitstates_epi(fn_epi_full, haiti.dat.epi) %>%
  dplyr::mutate(period = "epidemic")
fit_enddat <- process_fitstates(fn_end_full, haiti.dat.end, "cases_med") %>% 
  dplyr::filter(week > episettings$nweeks) %>%
  dplyr::mutate(period = "endemic")
fitdat <- bind_rows(fit_epidat, fit_enddat)

plt_fullfit <- plot_fit_full(fitdat)
ggsave(paste0(fig_dir, "fitfull_if_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_fullfit, width = 6, height = 4, units = "in")
print(plt_fullfit)

plt_fullfit_zm <- plt_fullfit +
  scale_x_continuous(limits = c(max(fitdat$week)-40, max(fitdat$week))) +
  scale_y_continuous(limits = c(0, 2000))
ggsave(paste0(fig_dir, "fitfullzm_if_", episettings$str_if, "_", endsettings$str_if, ".pdf"), plt_fullfit_zm, width = 4, height = 4, units = "in")
print(plt_fullfit_zm)

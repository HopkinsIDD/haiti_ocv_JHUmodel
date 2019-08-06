## rm(list=ls())
source("haiti_ocv_JHUmodel/Source/utils_clean.R")
reload_source()
source("haiti_ocv_JHUmodel/Source/src_generate_epi_params.R")
source("haiti_ocv_JHUmodel/Source/pomp_fit_epi.R")

## -------------------------------------- ##
## model settings ##
episettings <- fix_epi_settings()
nsamps <- episettings$nsamps
fig_dir <- paste0("haiti_ocv_JHUmodel/figures/", episettings$mcode, "/")
out_dir <- paste0("haiti_ocv_JHUmodel/GeneratedData/", episettings$mcode, "/")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

## -------------------------------------- ##
## data import and setup ##
pop.haiti <- get.haiti.pop()
haiti.dat <- get.mspp.agg.data() 
haiti.dat.epi <- haiti.dat %>% dplyr::filter(week <= episettings$nweeks)
covartab <- make.covartab(0, nrow(haiti.dat.epi)+1, byt=1, degree=6, nbasis=6, per=52.14)
num_betas <- ncol(covartab)-1

starts <- generate_epi_params(out_dir = out_dir, priors = episettings$priors)

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
process_fitstates <- function(fitstate_fns, hti_data){
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
## import epidemic fits from all tms
fn_tm <- list.files(out_dir, "fitstates_tm_")
fn_tm_full <- paste0(out_dir, fn_tm)
fit_tm_summ <- process_fitstates(fn_tm_full, haiti.dat.epi)

plt_tm_summ <- plot_fit(fit_tm_summ)
ggsave(paste0(fig_dir, "fit_tm_", episettings$str_tm, ".pdf"), plt_tm_summ, width = 4, height = 4, units = "in")
print(plt_tm_summ)

plt_tm_summ_zm <- plot_fit(fit_tm_summ) +
  scale_x_continuous(limits = c(episettings$nweeks-80,episettings$nweeks)) +
  scale_y_continuous(limits = c(0, 6000))
ggsave(paste0(fig_dir, "fitzm_tm_", episettings$str_tm, ".pdf"), plt_tm_summ_zm, width = 4, height = 4, units = "in")
print(plt_tm_summ_zm)

## -------------------------------------- ##
## import epidemic fits only from all converged tm
prof_tm <- readRDS(paste0(out_dir, "tmfits_", episettings$str_tm, ".rds")) %>%
  arrange(desc(loglik))
parid_ls <- prof_tm %>% dplyr::filter(conv==0) %>% dplyr::select(parid) %>% unlist %>% unname
fn_tm_conv <- paste0(out_dir, "fitstates_tm_", episettings$str_tm2, "_parid", parid_ls,  "_seed", episettings$seed, ".csv")
fit_tm_lim <- process_fitstates(fn_tm_conv, haiti.dat.epi)

plt_tm_lim <- plot_fit(fit_tm_lim)
ggsave(paste0(fig_dir, "fitconv_tm_", episettings$str_tm, ".pdf"), plt_tm_lim, width = 4, height = 4, units = "in")
print(plt_tm_lim)

## -------------------------------------- ##
## import epidemic fits from all ifs
fn_if <- list.files(out_dir, "fitstates_if_")
fn_if_full <- paste0(out_dir, fn_if)
fit_if_summ <- process_fitstates(fn_if_full, haiti.dat.epi)

plt_if_summ <- plot_fit(fit_if_summ)
ggsave(paste0(fig_dir, "fit_if_", episettings$str_if, ".pdf"), plt_if_summ, width = 4, height = 4, units = "in")
print(plt_if_summ)

plt_if_summ_zm <- plot_fit(fit_if_summ) +
  scale_x_continuous(limits = c(episettings$nweeks-80,episettings$nweeks)) +
  scale_y_continuous(limits = c(0, 6000))
ggsave(paste0(fig_dir, "fitzm_if_", episettings$str_if, ".pdf"), plt_if_summ_zm, width = 4, height = 4, units = "in")
print(plt_if_summ_zm)

## -------------------------------------- ##
## import epidemic fits from all converged ifs
prof_if <- readRDS(paste0(out_dir, "iffits_", episettings$str_if, ".rds"))

parid_conv_ls <- prof_if %>% filter(nfail.max == 0) %>% select(parid) %>% unlist %>% unname
fn_if_conv <- paste0(out_dir, "fitstates_if_", episettings$str_if2, "_parid", parid_conv_ls, "_seed", episettings$seed, ".csv")
fit_if_lim <- process_fitstates(fn_if_conv, haiti.dat.epi)

plt_if_lim <- plot_fit(fit_if_summ)
ggsave(paste0(fig_dir, "fitconv_if_", episettings$str_if, ".pdf"), plt_if_lim, width = 4, height = 4, units = "in")
print(plt_if_lim)


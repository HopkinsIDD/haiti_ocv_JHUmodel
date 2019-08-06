# Achieving coordinated national immunity and cholera elimination in Haiti through vaccination

## Model Supplement - Elizabeth C. Lee, Andrew S. Azman, Justin Lessler

This repository presents the code and data to reproduce the results for the model from the Johns Hopkins Bloomberg School of Public Health (Model-1). Please read the following descriptions about the contents of each folder.

### Supplement_Text

This folder contains documents describing the methods and results for Model-1.

### Source

This folder contains all of the R codes to fit the model, run diagnostics, and project the impact of mass OCV campaigns in Haiti. Model fits and projections were made using the R POMP package ($<$version 2).

The following steps will reproduce the model fit and projection pipeline:
1. Source `fit_epi_tm.R` in R.
1. Source `fit_epi_if.R` in R.
1. Source `fit_end_tm.R` in R.
1. Source `fit_end_if.R` in R.
1. Source `forecast_novac_stoch.R` and/or  Set `src_generate_fc_params.R` to the desired vaccination deployment scenario(s) and Source `forecast_vac_stoch.R` in R.
1. Source `main-outputs.R` and `write_elimdate_sims.R` to produce summaries of model outputs for comparison to other models.

#### Fit the model to reported Haiti incidence data
* `fit_epi_tm.R`: Use POMP trajectory matching to fit data for the epidemic period.
* `fit_epi_if.R`: Use POMP iterated filtering 2 algorithm to fit data for the epidemic period.
* `fit_end_tm.R`: Use POMP trajectory matching to fit data for the endemic period.
* `fit_end_if.R`: Use POMP iterated filtering 2 algorithm to fit data for the endemic period.
* `pomp_fit_epi.R`: POMP model for the epidemic period.
* `pomp_fit_end.R`: POMP model for the endemic period.
* `src_generate_epi_params.R`: Set parameter inputs for model fits to the epidemic period.
* `src_generate_end_params.R`: Set parameter inputs for model fits to the endemic period.

#### Explore model fits
* `main-explore-epi-fits.R`: Examine fits to epidemic period.
* `main-explore-end-fits.R`: Examine fits to endemic period.

#### Project the impact of mass OCV campaigns
* `forecast_novac_stoch.R`:
* `forecast_vac_stoch.R`: Use POMP to project a vaccination campaign.
* `pomp_fc_novac.R`: POMP model for status quo scenario projections.
* `pomp_fc_national.R`: POMP model for national vaccination campaign projections.
* `pomp_fc_2dept.R`: POMP model for 2-department vaccination campaign projections.
* `pomp_fc_3dept.R`: POMP model for 3-department vaccination campaign projections.
* `src_generate_fc_novac_params.R`: Fixed parameter inputs for model projections in the status quo scenario.
* `src_generate_fc_params.R`: Tunable parameter inputs for model projections in any given vaccination scenario.

#### Modeling utilities
* `utils_clean.R`: Utility functions to support model fitting and projection.
* `utils_eval.R`: Utility functions to support model diagnostics.
* `main-outputs.R`: Export model outputs in a format commonly agreed upon across modeling teams. 
* `write_elimdate_sims.R`: Export elimination dates for each simulation in a format commonly agreed upon across modeling teams.

### Data
* `haiti-data-from-2010-10-to-2019-01.csv`: Weekly department-level data of reported cholera cases in Haiti from October 2010 to January 2019.
* `ve_decay_bymonth.csv`: Vaccine effectiveness by month after vaccination for multiple vaccine effectiveness assumptions.

### GeneratedData
* `final/`: Generated data for fits to epidemic period
* `final1/`: Generated data for fits to endemic period
* `final1_id1/`: Generated data for projections of fast national OCV campaign
* `final1_id2/`: Generated data for projections of 2-department OCV campaign
* `final1_id3/`: Generated data for projections of slow national OCV campaign
* `final1_id4/`: Generated data for projections of 3-department OCV campaign
* `final1_id25/`: Generated data for projections of fast national high coverage OCV campaign
* `final1_novac/`: Generated data for projections of status quo scenario

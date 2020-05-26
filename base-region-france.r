library(rstan)
library(data.table)
library(lubridate,warn.conflicts = FALSE)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)

source("utils/arg-parser.r")
source("utils/read-interventions.r")
source("utils/process-covariates.r")
source("utils/process-covariates-region.r")
source("utils/ifr-tools.r")
source("utils/log-and-process.r")

VERSION="v5"

# Commandline options and parsing
parsedargs <- base_arg_parse()
DEBUG <- parsedargs[["DEBUG"]]
FULL <- parsedargs[["FULL"]]
StanModel <- parsedargs[["StanModel"]]
new_sub_folder <- parsedargs[["new_sub_folder"]]
max_date <- parsedargs[["max_date"]]
mobility_source <- parsedargs[["mobility_source"]]
formula_pooling <- parsedargs[["formula_pooling"]]
formula_partialpooling <- parsedargs[["formula_partialpooling"]]
zone_definition_file <- parsedargs[["activezones"]]

run_name <- create_analysis_folder(FULL, DEBUG, StanModel)

region_to_country_map <- read_country_region_map(zone_definition_file)

## Reading data from region file and world data and trimming it to max_date
data_files <- c(
  "data/COVID-19-up-to-date.rds",
  "data/all-france.rds"
)
countries <- names(region_to_country_map)
d <- read_obs_data(countries, data_files, max_date)
# Trim countries and regions that fail the number of death test.
death_thresh_epi_start = 10

trimmed_map = trim_country_map(
  d, region_to_country_map, max_date, death_thresh_epi_start)
region_to_country_map <- trimmed_map$region_to_country_map
countries <- names(region_to_country_map)

# Modelling + Forecasting
min_forecast <- 14
N2 <- trimmed_map$max_epi_data + min_forecast

## get IFR and population from same file
ifr.by.country <- return_ifr()
interventions <- read_interventions('data/interventions.csv', max_date)
mobility <- read_mobility(mobility_source, zone_definition_file)



formula = as.formula(formula_pooling)
formula_partial = as.formula(formula_partialpooling)
processed_data <- process_covariates_regions(
  region_to_country_map = region_to_country_map,
  mobility = mobility,
  interventions = interventions,
  d = d,
  ifr.by.country = ifr.by.country,
  N2 = N2,
  formula = formula,
  formula_partial = formula_partial,
  death_thresh_epi_start=death_thresh_epi_start
)

stan_data <- processed_data$stan_data
dates <- processed_data$dates
reported_deaths <- processed_data$deaths_by_country
reported_cases <- processed_data$reported_cases
infection_to_onset = processed_data$infection_to_onset
onset_to_death = processed_data$onset_to_death

log_simulation_inputs(run_name, region_to_country_map,  ifr.by.country,
    infection_to_onset, onset_to_death, VERSION, parsedargs)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(paste0('stan-models/',StanModel,'.stan'))


if(DEBUG) {
  fit = sampling(m,data=stan_data,iter=40,warmup=20,chains=2)
} else if (FULL) {
  fit = sampling(m,data=stan_data,iter=2000,warmup=1500,chains=4,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 15))
} else { 
  fit = sampling(m,data=stan_data,iter=600,warmup=300,chains=4,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 10))
}  


out = rstan::extract(fit)
estimated_cases_raw = out$prediction
estimated_deaths_raw = out$E_deaths
estimated_deaths_cf = out$E_deaths0

covariate_data = list(interventions, mobility)
save.image(paste0('results/',run_name,'.Rdata'))
save(
  fit, estimated_cases_raw, dates,reported_cases,reported_deaths,countries,
  region_to_country_map, estimated_deaths_raw, estimated_deaths_cf, 
  out,interventions,covariate_data, infection_to_onset, onset_to_death, VERSION,
  formula_pooling, formula_partialpooling,
  file=paste0('results/',run_name,'-stanfit.Rdata'))

postprocess_simulation(run_name, out, countries, dates)

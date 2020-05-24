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

regions <- read_country_file(parsedargs[["activeregions"]])
active_countries <- read_country_file(parsedargs[["activecountries"]])

run_name <- create_analysis_folder(FULL, DEBUG, StanModel)

region_to_country_map = list()
for(Region in regions){
  region_to_country_map[[Region]] <- "France"
}
for(Country in active_countries){
  region_to_country_map[[Country]] <- Country
}

## Reading data from region file and world data and trimming it to max_date
data_files <- c(
  "data/COVID-19-up-to-date.rds",
  "data/all-france.rds"
)
countries <- names(region_to_country_map)
d <- read_obs_data(countries, data_files, max_date)
# Trim countries and regions that fail the number of death test.
death_thresh_epi_start = 10
keep_regions = logical(length = length(region_to_country_map))
for(i in 1:length(region_to_country_map))
{
  Region <- names(region_to_country_map)[i]
  Country = region_to_country_map[[Region]]  
  d1=d[d$Country==Region,c(1,5,6,7)] 
  keep_regions[i] = !is.na(which(cumsum(d1$Deaths)>=death_thresh_epi_start)[1]) # also 5
  if (!keep_regions[i]) {
    message(sprintf(
      "WARNING: Region %s in country %s has not reached 10 deaths on %s, it cannot be processed\nautomatically removed from analysis\n",
      Region, Country, max_date))
  }
}
region_to_country_map <- region_to_country_map[keep_regions]
countries <- names(region_to_country_map)
## get IFR and population from same file
ifr.by.country <- return_ifr()
interventions <- read_interventions('data/interventions.csv', max_date)
mobility <- read_mobility(mobility_source)

# Modelling + Forecasting range needs serious revamp
forecast <- 21
N2 = 120
# N2 <- (max(d$DateRep) - min(d$DateRep) + 1 + forecast)[[1]]


formula = as.formula(formula_pooling)
formula_partial = as.formula(formula_partialpooling)
processed_data <- process_covariates_regions(
  regions = regions,
  mobility = mobility,
  intervention = interventions,
  d = d,
  ifr.by.country = ifr.by.country,
  N2 = N2,
  formula = formula,
  formula_partial = formula_partial
)

stan_data = processed_data$stan_data
dates = processed_data$dates
deaths_by_country = processed_data$deaths_by_country
reported_cases = processed_data$reported_cases
infection_to_onset = processed_data$infection_to_onset
onset_to_death = processed_data$onset_to_death

log_simulation_inputs(run_name, region_to_country_map,  ifr.by.country,
    infection_to_onset, onset_to_death, VERSION)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m = stan_model(paste0('stan-models/',StanModel,'.stan'))


if(DEBUG) {
  fit = sampling(m,data=stan_data,iter=40,warmup=20,chains=2)
} else if (FULL) {
  fit = sampling(m,data=stan_data,iter=1800,warmup=1000,chains=5,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 15))
} else { 
  fit = sampling(m,data=stan_data,iter=1000,warmup=500,chains=4,thin=1,control = list(adapt_delta = 0.95, max_treedepth = 10))
}   


out = rstan::extract(fit)
prediction = out$prediction
estimated.deaths = out$E_deaths
estimated.deaths.cf = out$E_deaths0

save.image(paste0('results/',run_name,'.Rdata'))


save(
  fit, prediction, dates,reported_cases,deaths_by_country,countries,
  region_to_country_map, estimated.deaths, estimated.deaths.cf, 
  out,interventions, infection_to_onset, onset_to_death, VERSION,
  file=paste0('results/',run_name,'-stanfit.Rdata'))

postprocess_simulation(run_name, out, countries, dates)

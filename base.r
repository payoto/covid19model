library(rstan)
library(data.table)
library(lubridate,warn.conflicts = FALSE)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(stringr)
library(abind)
library(optparse)
library(ggplot2)
library(ggrepel)
library(gtable)
library(zoo)

source('Italy/code/utils/read-data-subnational.r')
# provide functions for pre and post processing
source("utils/arg-parser.r")
source('utils/read-data.r')
source('utils/read-interventions.r')
source('utils/process-covariates.r')

VERSION="v5"

# Commandline options and parsing
# Needs cleaning
parsedargs <- base_arg_parse()
DEBUG <- parsedargs[["DEBUG"]]
FULL_RUN <- parsedargs[["FULL"]]
StanModel <- parsedargs[["StanModel"]]
if(length(cmdoptions$args) == 0) {
  StanModel = 'base-italy'
} else {
  StanModel = cmdoptions$args[1]
}

args = cmdoptions$args
# if using rstudio change this, 1: stan-file, 2: mobility data to use, 3:interventions to use
# 4:formula to use for full pool, 5:formula for partial pool
# debug and full are either commandline or sys variable
# these are the intervention names
#"schools_universities" "public_events""lockdown" "social_distancing_encouraged" "self_isolating_if_ill"  
# these are the mobility names
# "grocery"  "parks"  "residential"  "retail" "transit"   "workplace"  "averageMobility"
if(length(args) == 0) {
  args = c('base-italy', 'google', 'interventions',
           '~ -1 + residential + transit + averageMobility',
           '~ -1 + residential + transit + averageMobility'
  ) 
} 

StanModel = args[1]
cat(sprintf("Running:\nStanModel = %s\nMobility = %s\nInterventions = %s\nFixed effects:%s\nRandom effects:%s\nDebug: %s\n",
            StanModel,args[2],args[3], args[4],args[5], DEBUG))

# Read which countires to use
countries <- read.csv('data/regions.csv', stringsAsFactors = FALSE)
# Read deaths data for regions
d <- read_obs_data()
regions <- countries

# Read ifr 
ifr.by.country <- read_ifr_data()

# Read google mobility, apple mobility, interventions, stringency
google_mobility <- read_google_mobility("Italy")
mobility<-google_mobility[which(google_mobility$country!="Italy"),]

# Read interventions
interventions <- read_interventions('data/interventions.csv')


# Table 1 and top 7
regions_sum <- d %>% group_by(country) %>% summarise(Deaths=sum(Deaths)) %>%
  inner_join(ifr.by.country) %>% mutate(deathsPer1000=Deaths/popt) %>% 
  arrange(desc(deathsPer1000))
regions_sum <- regions_sum[,-which(colnames(regions_sum) %in% c("X"))]
regions_sum$ifr<-signif(regions_sum$ifr*100,2)
regions_sum$deathsPer1000 <- signif(regions_sum$deathsPer1000*1000,2)

top_7 <- regions_sum[1:7,]

forecast <- 7 # increase to get correct number of days to simulate
# Maximum number of days to simulate
N2 <- (max(d$DateRep) - min(d$DateRep) + 1 + forecast)[[1]]

formula = as.formula(args[4])
formula_partial = as.formula(args[5])
processed_data <- process_covariates(regions = regions, mobility = mobility, intervention = interventions, 
                                     d = d , ifr.by.country = ifr.by.country, N2 = N2, formula = formula, formula_partial = formula_partial)

stan_data <- processed_data$stan_data
dates <- processed_data$dates
reported_deaths <- processed_data$deaths_by_country
reported_cases <- processed_data$reported_cases

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

out <- rstan::extract(fit)
estimated_cases_raw <- out$prediction
estimated_deaths_raw <- out$E_deaths
estimated_deaths_cf <- out$E_deaths0

JOBID = Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID = as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s",JOBID))

countries <- countries$Regions
covariate_data = list(interventions, mobility)
save.image(paste0('results/',StanModel,'-',JOBID,'.Rdata'))
save(fit, prediction, dates,reported_cases,deaths_by_country,countries, estimated_deaths_raw, estimated_deaths_cf, 
  out,interventions, VERSION,file=paste0('results/',StanModel,'-',JOBID,'-stanfit.Rdata'))

## Ensure that output directories exist
dir.create("results/", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/", showWarnings = FALSE, recursive = TRUE)
dir.create("web/", showWarnings = FALSE, recursive = TRUE)
dir.create("web/data", showWarnings = FALSE, recursive = TRUE)

filename <- paste0(StanModel,'-',JOBID)

regions <- unique(d$country)


# This is a hack to get it to save
states = regions

covariate_data = list(interventions, mobility)

save(fit, dates, reported_cases, reported_deaths, regions, states, JOBID ,
     estimated_cases_raw, estimated_deaths_raw, estimated_deaths_cf, stan_data, covariate_data,
     file=paste0('Italy/results/',StanModel,'-',JOBID,'-stanfit.Rdata'))

source("Italy/code/plotting/make-plots.r")
make_plots_all(paste0('Italy/results/', StanModel, '-', JOBID, '-stanfit.Rdata'), 
               last_date_data = max(dates[[1]]))
source("Italy/code/utils/make-table.r")
# Prints attackrates to console
make_table(paste0('Italy/results/', StanModel, '-', JOBID, '-stanfit.Rdata'), 
           date_till_percentage = max(dates[[1]]))

# code for scenarios runs only in full mode
if (FULL){
  source("Italy/code/utils/simulate-regional.r")
  len_forecast <- 8*7
  
  # Can make plots = TRUE if you want to see 3 panel plots for simulations and rt_plots
  simulate_scenarios(JOBID = JOBID,  StanModel, plots = TRUE, scenario_type = "increase-mob-current", len_forecast = len_forecast,
                     subdir='Italy',
                     simulate_code='Italy/code/stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 40)
  simulate_scenarios(JOBID = JOBID,  StanModel, plots = TRUE, scenario_type = "increase-mob-current", len_forecast = len_forecast,
                     subdir='Italy',
                     simulate_code='Italy/code/stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 20)
  simulate_scenarios(JOBID = JOBID,  StanModel, plots = TRUE, scenario_type = "constant-mob", len_forecast = len_forecast,
                     subdir='Italy',
                     simulate_code='Italy/code/stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 0)
  
  source("Italy/code/utils/make-table.r")
  # Prints attack rates to console
  scenario_type = "constant-mob"
  mobility_increase = 0
  make_table_simulation(paste0('Italy/results/sim-', scenario_type, "-", StanModel, '-', len_forecast, '-', mobility_increase, '-', JOBID, '-stanfit.Rdata'), 
                        date_till_percentage = max(dates[[1]]) + len_forecast)
  
  
  source("Italy/code/plotting/make-scenario-plots-top7.r")
  
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 20,top=7)
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 40,top=7)
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 20,top=8)
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 40,top=8)
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 20,top=9)
  make_scenario_comparison_plots_mobility(JOBID = JOBID, StanModel, len_forecast = len_forecast, 
                                          last_date_data = max(dates[[1]]) + len_forecast, baseline = FALSE, 
                                          mobility_increase = 40,top=9)
print("Generating covariate size effects plot")
covariate_size_effects_error <- system(paste0("Rscript covariate-size-effects.r ", filename,'-stanfit.Rdata'),intern=FALSE)
if(covariate_size_effects_error != 0){
  stop(sprintf("Error while plotting covariate size effects! Code: %d", covariate_size_effects_error))
}

mu = (as.matrix(out$mu))
colnames(mu) = countries
g = bayesplot::mcmc_intervals(mu,prob = .9)
ggplot2::ggsave(sprintf("results/%s-mu.png",filename),g,width=4,height=6)
tmp = lapply(1:length(countries), function(i) (out$Rt_adj[,stan_data$N[i],i]))
Rt_adj = do.call(cbind,tmp)
colnames(Rt_adj) = countries
g = bayesplot::mcmc_intervals(Rt_adj,prob = .9)
ggsave(sprintf("results/%s-final-rt.png",filename),g,width=4,height=6)

print("Generate 3-panel plots")
plot_3_panel_error <- system(paste0("Rscript plot-3-panel.r ", filename,'-stanfit.Rdata'),intern=FALSE)
if(plot_3_panel_error != 0){
  stop(sprintf("Generation of 3-panel plots failed! Code: %d", plot_3_panel_error))
}

print("Generate forecast plot")
plot_forecast_error <- system(paste0("Rscript plot-forecast.r ",filename,'-stanfit.Rdata'),intern=FALSE)
if(plot_forecast_error != 0) {
  stop(sprintf("Generation of forecast plot failed! Code: %d", plot_forecast_error))
}

print("Make forecast table")
make_table_error <- system(paste0("Rscript make-table.r results/",filename,'-stanfit.Rdata'),intern=FALSE)
if(make_table_error != 0){
  stop(sprintf("Generation of alpha covar table failed! Code: %d", make_table_error))
}


verify_result_error <- system(paste0("Rscript web-verify-output.r ", filename,'.Rdata'),intern=FALSE)
if(verify_result_error != 0){
  stop(sprintf("Verification of web output failed! Code: %d", verify_result_error))
}

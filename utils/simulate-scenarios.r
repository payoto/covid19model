source("Italy/code/utils/simulate-regional.r")
source("Italy/code/utils/make-table.r")
source("Italy/code/plotting/make-scenario-plots-top7.r")

mobility_scenarios <- function (run_name){
  len_forecast <- 8*7
  
  # Can make plots = TRUE if you want to see 3 panel plots for simulations and rt_plots
  simulate_scenarios(run_name, plots = TRUE, scenario_type = "increase-mob-current", len_forecast = len_forecast,
                     subdir='.',
                     simulate_code='stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 40)
  simulate_scenarios(run_name, plots = TRUE, scenario_type = "increase-mob-current", len_forecast = len_forecast,
                     subdir='.',
                     simulate_code='stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 20)
  simulate_scenarios(run_name, plots = TRUE, scenario_type = "constant-mob", len_forecast = len_forecast,
                     subdir='.',
                     simulate_code='stan-models/simulate.stan', mobility_vars=c(1,2,3), 
                     mobility_increase = 0)
  
  # Prints attack rates to console
  scenario_type = "constant-mob"
  mobility_increase = 0
  make_table_simulation(paste0('results/', run_name ,'sim-', scenario_type, "-", StanModel, '-', len_forecast, '-', mobility_increase, '-', JOBID, '-stanfit.Rdata'), 
                        date_till_percentage = max(dates[[1]]) + len_forecast)
  
  
  
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
}
library(optparse)


base_arg_parse <- function (){
	# Commandline options and parsing
	parser <- OptionParser()
	parser <- add_option(parser, c("-D", "--debug"), action="store_true",
	                     help="Perform a debug run of the model")
	parser <- add_option(parser, c("-F", "--full"), action="store_true",
	                     help="Perform a full run of the model")
	parser <- add_option(parser, c("--nosubdir"), action="store_true",
	                     help="Do not create subdirectories for generated data.")
	parser <- add_option(parser, c("--maxdate"), default="",
	                     help="Consider only data up to max date 'dd/mm/yy' format.")
	parser <- add_option(parser, c("--activezones"), default="config/covid19model_analysis_zones.csv",
	                     help="Parameter containing the active countries.")
	parser <- add_option(parser, c("--activeregions"), default="config/active-regions.cfg",
	                     help="Parameter containing the active regions.")
	parser <- add_option(parser, c("--activecountries"), default="config/active-countries.cfg",
	                     help="Parameter containing the active countries.")
	
	cmdoptions <- parse_args(parser, args = commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)

	# Default run parameters for the model
	if(is.null(cmdoptions$options$debug)) {
	  DEBUG = Sys.getenv("DEBUG") == "TRUE"
	} else {
	  DEBUG = cmdoptions$options$debug
	}

	if(is.null(cmdoptions$options$full)) {
	  FULL = Sys.getenv("FULL") == "TRUE"
	} else {
	  FULL = cmdoptions$options$full
	}

		
	new_sub_folder = "TRUE"
	if(!is.null(cmdoptions$options$nosubdir)){
		new_sub_folder = !cmdoptions$options$nosubdir
	}

	if(DEBUG && FULL) {
	  stop("Setting both debug and full run modes at once is invalid")
	}

    std_args <- c( # Default ordered arguments
        StanModel = 'base-italy',
        mobility_source = 'google',
        intervention_source = 'interventions',
        formula_pooling = '~ -1 + residential + transit + averageMobility',
        formula_partialpooling = '~ -1 + residential + transit + averageMobility'
    )
    if (length(cmdoptions$args)>0){
	    for (i in 1:length(cmdoptions$args)){ # overwerite defaults
	        std_args[i] = cmdoptions$args[i]
	    }	 
    }

	parsedargs <- c(
			DEBUG=DEBUG,
			FULL=FULL,
			new_sub_folder=new_sub_folder ,
			max_date = cmdoptions$options$maxdate,
			activeregions = cmdoptions$options$activeregions,
			activecountries = cmdoptions$options$activecountries,
			activezones = cmdoptions$options$activezones,
            std_args
		)
    message("Configured values are:")
    for (i in 1:length(parsedargs)){
        message(sprintf("\t%s : %s",names(parsedargs)[i], parsedargs[i]))
    }
    # print(names(parsedargs))
    # print(parsedargs)
	return(parsedargs)
}

read_country_file <- function (filename){
	countries <- scan(filename, what="", sep="\n")
	for (i in 1:length(countries)){
		countries[i] = trimws(countries[i])
	}
	
	return(countries)
}

trim_data_to_date_range <- function (data, max_date, date_field="DateRep", 
	format_field='%d/%m/%Y', format_max='%d/%m/%y'){
    
	if (max_date == "" || is.null(max_date)){
		return(data)
	} else {

		return (data[
			as.Date(data[[date_field]], format=format_field) 
				<= as.Date(max_date, format=format_max),
			])
	}
}

create_analysis_folder <- function(FULL, DEBUG, StanModel){
    JOBID = Sys.getenv("PBS_JOBID")
    if(JOBID == "")
      JOBID = as.character(abs(round(rnorm(1) * 1000000)))
    print(sprintf("Jobid = %s",JOBID))
    fullstr <- "short"
    if (FULL){
      fullstr <- "fullrun"
    } else if (DEBUG) {
      fullstr <- "debug"
    }

    run_name <- paste0(StanModel,'-',fullstr,'-', format(Sys.time(), '%Y%m%dT%H%M%S'),'-',JOBID)
    if (new_sub_folder){
      result_folders <- c(
          "results", "figures"
        )
      for (fold in result_folders){
        dir.create(paste0(fold ,'/', run_name))
      }
      run_name <- paste0(run_name ,'/', run_name)
    }
    return(run_name)
}
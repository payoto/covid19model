library(tidyr)
library(lubridate)
library(stringr)
library(dplyr)

source("utils/arg-parser.r")

read_obs_data <- function(countries, file_list=c('data/COVID-19-up-to-date.rds'), max_date=""){
  # Read the deaths and cases data
  d <- trim_data_to_date_range(
    do.call('rbind', lapply(file_list, readRDS)),
    max_date  # optional arguments allow data customisation
  )
  colnames(d)[colnames(d) == "Countries.and.territories"] <- "Country"
  tryCatch({
    d <-d[d$Country %in% countries$Regions, c(1,5,6,7)]
  },error = function(e) {
    d <-d[d$Country %in% names(countries), c(1,5,6,7)]
  })
  d$DateRep <- as.Date(d$DateRep, format = '%d/%m/%Y')
  return(d)
}

read_ifr_data <- function(){
  ifr.by.country <- read.csv("data/popt_ifr.csv")
  ifr.by.country$country <- as.character(ifr.by.country[,2])
  ifr.by.country$country[ifr.by.country$country == "United Kingdom"] <- "United_Kingdom"
  return(ifr.by.country)
  
}

translate_region_names <- function(
  data,
  translation_file="config/covid19model_analysis_zones.csv",
  data_name_field="country",
  model_name_field = "covid19model_region",  # the column in the translation csv with the target names
  mobility_name_field = "GCMR_region"  # the column in the translation csv with the source names
){
  # Need revamp
  # stop("rad-data.r : translate region names not coded")

  nametrans <- read.csv(translation_file)
  # Replaces the names in the mobility data with that of the model 
  for (i in 1:length(nametrans[[mobility_name_field]])) {
    GCMR_name = nametrans[[mobility_name_field]][i]
    model_name = nametrans[[model_name_field]][i]
    data[[data_name_field]][data[[data_name_field]]==GCMR_name] <- model_name
  }
  return(data)
}

select_and_translate_region_names <- function(
  data,
  selection_file,
  data_name_field="country",
  model_name_field = "covid19model_region",  # the column in the translation csv with the target names
  mobility_name_field = "GCMR_region"  # the column in the translation csv with the source names
){
  # Need revamp
  # stop("rad-data.r : translate region names not coded")

  nametrans <- read.csv(selection_file)
  # Replaces the names in the mobility data with that of the model 
  data_out <- c()
  for (i in 1:length(nametrans[[mobility_name_field]])) {
    GCMR_name = nametrans[[mobility_name_field]][i]
    model_name = nametrans[[model_name_field]][i]
    temp <- data[data[[data_name_field]]==GCMR_name,]
    temp[[data_name_field]] <- model_name
    data_out <- bind_rows(data_out, temp)
  }
  return(data_out)
}

read_google_mobility <- function(Country="Italy"){
  google_mobility <- read.csv('data/Global_Mobility_Report.csv', stringsAsFactors = FALSE)
  google_mobility$date = as.Date(google_mobility$date, format = '%Y-%m-%d')
  google_mobility[, c(6,7,8,9,10,11)] <- google_mobility[, c(6,7,8,9,10,11)]/100
  google_mobility[, c(6,7,8,9,10)] <- google_mobility[, c(6,7,8,9,10)] * -1
  google_mobility<-google_mobility[,c(2,3,5,6,7,8,9,10,11)]
  colnames(google_mobility)[which(colnames(google_mobility)=="country_region")]<-"state"
  colnames(google_mobility)[which(colnames(google_mobility)=="sub_region_1")]<-"country"
  colnames(google_mobility)[which(colnames(google_mobility)=="grocery_and_pharmacy_percent_change_from_baseline")]<-"grocery.pharmacy"
  colnames(google_mobility)[which(colnames(google_mobility)=="parks_percent_change_from_baseline")]<-"parks"
  colnames(google_mobility)[which(colnames(google_mobility)=="transit_stations_percent_change_from_baseline")]<-"transitstations"
  colnames(google_mobility)[which(colnames(google_mobility)=="workplaces_percent_change_from_baseline")]<-"workplace"
  colnames(google_mobility)[which(colnames(google_mobility)=="residential_percent_change_from_baseline")]<-"residential"
  colnames(google_mobility)[which(colnames(google_mobility)=="retail_and_recreation_percent_change_from_baseline")]<-"retail.recreation"
  google_mobility$country[which(google_mobility$country =="")]<-google_mobility$state[which(google_mobility$country =="")]
  google_mobility <- google_mobility %>% select(
    "country",
    "date",
    "grocery.pharmacy",
    "parks",
    "residential",
    "retail.recreation",
    "transitstations",
    "workplace"
  ) 
  
  return (google_mobility)
}

read_mobility <- function(
  mobility_source='google', 
  selection_file="config/covid19model_analysis_zones.csv"
){
  if(mobility_source == 'google'){
    mobility = read_google_mobility()
  } else {
    stop(sprintf("%s is not a recognised mobility_source.", mobility_source))
  }
  mobility = select_and_translate_region_names(mobility, selection_file)
  # mobility<-google_mobility[which(google_mobility$country!="Italy"),]
  mobility <- mobility %>% select(
    "country",
    "date",
    "grocery.pharmacy",
    "parks",
    "residential",
    "retail.recreation",
    "transitstations",
    "workplace"
  ) 
  return(mobility)
}


read_country_region_map <- function (zone_definition_file){
  analysis_zones <- read.csv(zone_definition_file)

  region_to_country_map = list()
  for(i in 1:length(analysis_zones$covid19model_region)){
    region_to_country_map[[analysis_zones$covid19model_region[i]]] <- analysis_zones$covid19model_country[i]
  }

  return(region_to_country_map)
}

trim_country_map <- function(
  d,
  region_to_country_map,
  max_date=Sys.Date(),
  death_thresh_epi_start=10,
  epi_death_offset=30
){
  # Removes contry map with not enough data 
  # (a total number of deaths below `death_thresh_epi_start`)
  # And returns the length of the longest epidemic
  keep_regions = logical(length = length(region_to_country_map))
  max_epi_data = 0
  for(i in 1:length(region_to_country_map))
  {
    Region <- names(region_to_country_map)[i]
    Country = region_to_country_map[[Region]]  
    d1=d[d$Country==Region,c(1,5,6,7)]
    index_10_deaths = which(cumsum(d1$Deaths)>=death_thresh_epi_start)[1]
    keep_regions[i] = !is.na(index_10_deaths) # also 5
    if (!keep_regions[i]) {
      message(sprintf(
        "WARNING: Region %s in country %s has not reached 10 deaths on %s, it cannot be processed\nautomatically removed from analysis\n",
        Region, Country, max_date))
    } else {
      region_epidemic_start = d1$DateRep[index_10_deaths] - days(epi_death_offset)
      region_epidemic_length = max(d1$DateRep) - region_epidemic_start[1] + 1
      max_epi_data = max(max_epi_data, region_epidemic_length[1])
    }
  }

  region_to_country_map <- region_to_country_map[keep_regions]
  return (list(
    region_to_country_map=region_to_country_map,
    max_epi_data=max_epi_data
  ))
}
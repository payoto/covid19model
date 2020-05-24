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
  translation_file="Italy/data/province_name_translation.csv",
  field_name="country",
  trans_field_start="google_county",
  trans_field_end="country"
){
  # Need revamp
  stop("rad-data.r : translate region names not coded")

  nametrans <- read.csv(translation_file)
  # add Italy row to translation layer
  Italy<-data.frame(denominazione_regione=Country,google_county=Country,county=Country)
  nametrans<-bind_rows(nametrans,Italy)
  mobility$country<-as.factor(mobility$country)
  nametrans$google_county<-as.factor(nametrans$google_county)
  #mobility <- mobility %>% filter(country !="")
  colnames(nametrans)[which(colnames(nametrans)=="google_county")]<-"country"
  mobility <- inner_join(mobility,nametrans,by.x="country",by.y="country") # fix names of regions
  mobility$country<-str_replace_all(mobility$denominazione_regione, " ", "_")
  mobility <- mobility %>% select(country,date,grocery.pharmacy,parks,residential,retail.recreation,transitstations,workplace)
  
  # Changing region names
  nametrans <- read.csv("Italy/data/province_name_translation.csv")
  colnames(nametrans)[which(colnames(nametrans)=="denominazione_regione")]<-"country"
  nametrans$country<-as.factor(nametrans$country)
  nametrans$country<-str_replace_all(nametrans$country, " ", "_")
  mobility <- mobility %>% filter(country !=Country)
  mobility <- inner_join(mobility,nametrans,by.x="country",by.y="country") # fix names of regions
  mobility <- mobility[,-which(colnames(mobility) %in% c("country","google_county"))]
  colnames(mobility)[which(colnames(mobility)=="county")] <- "country"
  mobility$country<-str_replace_all(mobility$country, " ", "_")
  mobility <- mobility %>% select("country","date","grocery.pharmacy","parks","residential","retail.recreation","transitstations","workplace") 
  return(mobility)
}


read_google_mobility <- function(Country="Italy"){
  google_mobility <- read.csv('Italy/data/Global_Mobility_Report.csv', stringsAsFactors = FALSE)
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
  google_mobility$country[which(google_mobility$country =="")]<-Country
  colnames <- colnames %>% select(
    "state",
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

read_mobility <- function(mobility_source){
  if(mobility_source == 'google'){
    mobility = read_google_mobility()
  } else {
    stop(sprintf("%s is not a recognised mobility_source.", mobility_source))
  }
  mobility = translate_region_names(mobility)
  # mobility<-google_mobility[which(google_mobility$country!="Italy"),]
  return(mobility)
}

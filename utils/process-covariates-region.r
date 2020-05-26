library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(scales)
library(stringr)
library(abind)
library(zoo)
library(tidyverse)
library(magrittr)


process_covariates_regions <- function(
  region_to_country_map,
  mobility,
  interventions,
  d,
  ifr.by.country,
  N2,
  formula,
  formula_partial,
  death_thresh_epi_start=10
){
  interventions$Country <- factor(interventions$Country)
  
  serial.interval = read.csv("data/serial_interval.csv")
  # Pads serial interval with 0 if N2 is greater than the length of the serial
  # interval array
  if (N2 > length(serial.interval$fit)) {
    pad_serial.interval <- data.frame(
      "X"=(length(serial.interval$fit)+1):N2,
      "fit"=rep(1e-17, max(N2-length(serial.interval$fit), 0 ))
    )
    serial.interval = rbind(serial.interval, pad_serial.interval)
  }
  # various distributions required for modeling
  infection_to_onset <- c("mean"=5.1, "deviation"=0.86)
  onset_to_death <- c("mean"=18.8, "deviation"=0.45)
  # infection-to-onset distribution
  x1 = rgammaAlt(1e6,infection_to_onset["mean"], infection_to_onset["deviation"])
  # onset-to-death distribution
  x2 = rgammaAlt(1e6,onset_to_death["mean"], onset_to_death["deviation"])
  ecdf.saved <- ecdf(x1+x2)
  # stan data definition
  stan_data <- list(M=length(names(region_to_country_map)),N=NULL,deaths=NULL,f=NULL,
                    N0=6,cases=NULL,SI=serial.interval$fit[1:N2],
                    EpidemicStart = NULL, pop = NULL)
  # Other values
  forecast <- 0
  dates <- list()
  reported_cases <- list()
  deaths_by_country <- list()
  intervention_length <- length(colnames(interventions))
  mobility_length <- length(colnames(mobility))
  
  covariate_list <- list()
  covariate_list_partial <- list()
  processed_mobility <-list()
  preprocess_error = FALSE

  k=1
  # going over each region
  for(Region in names(region_to_country_map))
  {
    Country = region_to_country_map[[Region]]
    print(sprintf("Region: %s in country: %s ",Region,Country))
    if(any(ifr.by.country$country == Region)){
        # to add
      IFR <- ifr.by.country$ifr[ifr.by.country$country == Country]
      region_pop <- ifr.by.country[ifr.by.country$country==Country,]
    } else {
        IFR <- ifr.by.country$ifr[ifr.by.country$country == Country]
        region_pop <- ifr.by.country[ifr.by.country$country==Country,]

    }

    if(!any(d$Country==Region)){
      preprocess_error = TRUE
      message(sprintf(
        "ERROR: Region %s in country %s had no data (region length(region)==0)", Region, Country))
      next
    }
    region <- d[d$Country==Region,]
    region$DateRep <-region$DateRep
    region <-region[order(as.Date(region$DateRep)),]  # ensure date ordering

    mobility1 <- mobility[mobility$country==Region,]
    if(!any(mobility$country==Region)){
      preprocess_error = TRUE
      message(sprintf(
        "ERROR: Region %s in country %s had no mobility data (region length(mobility1)==0)", Region, Country))
      next
    }
    mobility1$date<-ymd(mobility1$date)
    mobility1<-na.locf(mobility1)
    mobility1 <- mobility1[order(mobility1$date),]  # ensure date ordering
    # padding in raw data backwards ex. portugal
    date_min <- dmy('31/12/2019') 
    if (region$DateRep[1] > date_min){
      print(paste(Region,'In padding death data'))
      pad_days <-region$DateRep[1] - date_min
      pad_dates <- date_min + days(1:pad_days[[1]]-1)
      padded_data <- data.frame("Country" = rep(Region, pad_days),
                                "DateRep" = as.Date(pad_dates,format="%Y-%m-%d"),
                                "Cases" = as.integer(rep(0, pad_days)),
                                "Deaths" = as.integer(rep(0, pad_days)),
                                stringsAsFactors=F)
      
     region <- bind_rows(padded_data,region)
    }
    
    # replace NA in mobility data
    mobility1$grocery.pharmacy <- na.locf(mobility1$grocery.pharmacy)
    mobility1$residential <- na.locf(mobility1$residential)
    mobility1$parks <- na.locf(mobility1$parks)
    mobility1$workplace <- na.locf(mobility1$workplace)
    mobility1$retail.recreation <- na.locf(mobility1$retail.recreation)
    mobility1$transitstations <- na.locf(mobility1$transitstations)
    
    # Padding in mobility data for dates before first time data exists
    if (mobility1$date[1] > date_min){
      print(paste(Region,'In padding mobility backwards'))
      pad_days <- mobility1$date[1] - date_min
      pad_dates <- date_min + days(1:pad_days[[1]]-1)
      PAD = rep(0,pad_days)
      if (mobility_length==8){
        padded_data <- data.frame("country" = rep(Region, pad_days),
                                  "date" = pad_dates,
                                  "grocery.pharmacy" = PAD, #rep(mobility1$grocery.pharmacy[1],pad_days),
                                  "parks" = PAD, #rep(mobility1$parks[1],pad_days),
                                  "residential" = PAD, #rep(mobility1$residential[1],pad_days),
                                  "retail.recreation" =PAD, #rep(mobility1$retail.recreation[1],pad_days),
                                  "transitstations" = PAD, #rep(mobility1$transitstations[1],pad_days),
                                  "workplace" = PAD, #rep(mobility1$workplace[1],pad_days),
                                  stringsAsFactors=F)
      }
      mobility1 <- bind_rows(padded_data, mobility1)
    }
    mobility1 = mobility1[order(mobility1$date),]  # ensure date ordering
    
    # forecasting mobility data
    date_max = max(region$DateRep)
    date_mobility_max = max(mobility1$date)
    
    if (date_mobility_max < date_max){
      print(paste(Region,'In padding mobility forward'))
      forecast_days <- date_max - date_mobility_max
      forecast_dates <- date_mobility_max + days(1:forecast_days[[1]])
      forecast_data <- data.frame("country"=rep(Region,forecast_days[[1]]),"date" = as.Date(forecast_dates,format="%Y-%m-%d"))
      fore<-list()
      
      for(i in 3:mobility_length){
        mob<- mobility1 %>% select(date,colnames(mobility[i]))
        
        # impute last week
        #fore[[i]] <- c(0:(as.numeric(forecast_days)-1)) %>% map(function(s) ts[length(ts)-6+mod(s,7)])
        f <- c(0:(as.numeric(forecast_days)-1)) %>% map(function(s) mob[nrow(mob)-6+mod(s,7),]) %>% bind_rows()
        fore[[i]]<-f[,2]
      }
      
      padded_data <- data.frame("country" = rep(Region, forecast_days[[1]]),
                                "date" = forecast_dates,
                                "grocery.pharmacy" = fore[[3]],
                                "parks" = fore[[4]],
                                "residential" = fore[[5]],
                                "retail.recreation" = fore[[6]],
                                "transitstations" = fore[[7]],
                                "workplace" = fore[[8]],
                                stringsAsFactors=F)
      mobility1 <- bind_rows(mobility1,padded_data)
    }
    mobility1 = mobility1[order(mobility1$date),]  # ensure date ordering
    mobility1<-na.locf(mobility1)
    
    index = which(region$Cases>0)[1]
    index1 = which(cumsum(region$Deaths)>=death_thresh_epi_start)[1] # also 5
    if (is.na(index1)) {
      preprocess_error = TRUE
      message(sprintf(
        "ERROR: Region %s in country %s has not reached 10 deaths on %s, it cannot be processed\nremove from 'active-countries.cfg' or 'active-regions.cfg'\n",
        Region, Country, max_date))
      next
    }
    index2 = index1-30
    
    print(sprintf("First non-zero cases is on day %d, and 30 days before 10 deaths is day %d",index,index2))
    region=region[index2:nrow(region),]
    stan_data$EpidemicStart = c(stan_data$EpidemicStart,index1+1-index2)
    stan_data$pop = c(stan_data$pop,region_pop$popt)
    mobility1 = mobility1[index2:nrow(mobility1),]
    # NPI interventionss are being used
    interventions_region <- interventions[interventions$Country == Country, c(2,3,4,5,6)] # school, self-isolation, public, lockdown, social-distancing
    for (ii in 1:ncol(interventions_region)) {
      covariate = names(interventions_region)[ii]
      region[covariate] <- (region$DateRep >= as.data.frame(interventions_region)[1,covariate])*1  # should this be > or >=?
    }
    
    dates[[Region]] =region$DateRep
    # hazard estimation
    N = length(region$Cases)
    print(sprintf("%s has %d days of data",Region,N))
    forecast = N2 - N
    # IFR is the overall probability of dying given infection
    convolution = function(u) (IFR * ecdf.saved(u))
    
    f = rep(0,N2) # f is the probability of dying on day i given infection
    f[1] = (convolution(1.5) - convolution(0))
    for(i in 2:N2) {
      f[i] = (convolution(i+.5) - convolution(i-.5)) 
    }
    reported_cases[[Region]] = as.vector(as.numeric(region$Cases))
    deaths=c(as.vector(as.numeric(region$Deaths)),rep(-1,forecast))
    cases=c(as.vector(as.numeric(region$Cases)),rep(-1,forecast))
    deaths_by_country[[Region]] = as.vector(as.numeric(region$Deaths))
    region_intervention <- as.data.frame(region[, colnames(interventions_region)])
    # This line will prevent modelling of a future change of interventions
    region_intervention[N:(N+forecast),] <- region_intervention[N,]
    mobility1[N:(N+forecast),] <- mobility1[N,]
    ## append data
    stan_data$N = c(stan_data$N,N)
    # stan_data$x = cbind(stan_data$x,x)
    stan_data$f = cbind(stan_data$f,f)
    stan_data$deaths = cbind(stan_data$deaths,deaths)
    stan_data$cases = cbind(stan_data$cases,cases)
    
    stan_data$N2=N2
    stan_data$x=1:N2
    if(length(stan_data$N) == 1) {
      stan_data$N = as.array(stan_data$N)
    }
    #parsing features
    df_features <- data.frame('schools_universities' = region_intervention$schools_universities, 
                              'self_isolating_if_ill' = region_intervention$self_isolating_if_ill, 
                              'social_distancing_encouraged' = region_intervention$social_distancing_encouraged, 
                              'public_events' = region_intervention$public_events,
                              'lockdown' = region_intervention$lockdown,
                              'first' = 1*((region_intervention$schools_universities+
                                              region_intervention$self_isolating_if_ill+
                                              region_intervention$social_distancing_encouraged+
                                              region_intervention$public_events+
                                              region_intervention$lockdown) >= 1),
                              'residential' = mobility1$residential, 
                              'transit' = mobility1$transitstations, 
                              'grocery' = mobility1$grocery.pharmacy,
                              'parks' = mobility1$parks,
                              'retail' = mobility1$retail.recreation,
                              'workplace' = mobility1$workplace,
                              'averageMobility' = (mobility1$grocery.pharmacy + mobility1$parks +
                                                     mobility1$retail.recreation + mobility1$workplace)/4)
    features <- model.matrix(formula, df_features)
    features_partial <- model.matrix(formula_partial, df_features)
    processed_mobility <- bind_rows(processed_mobility, bind_cols(mobility1,region_intervention))
    covariate_list[[k]] <- features
    covariate_list_partial[[k]] <- features_partial
    k <- k+1
  }
  
  stan_data$P = dim(features)[2]
  stan_data$X = array(NA, dim = c(stan_data$M , stan_data$N2 ,stan_data$P ))
  stan_data$P_partial = dim(features_partial)[2]
  # If there is no partial pooling of variables set of dim 1 to avoid overflows 
  if(stan_data$P_partial==0){
    stan_data$X_partial = array(0, dim = c(stan_data$M , stan_data$N2, 1))
  }
  # If there partial pooling of variables set X_partial to the correct size:
  # (countries or regions)* N2 * P_partial 
  else{
    stan_data$X_partial = array(NA, dim = c(stan_data$M , stan_data$N2 ,stan_data$P_partial))
  }
  # Unpack covariates and partial covariates to pass to stan model
  for (i in 1:stan_data$M){
    stan_data$X[i,,] = covariate_list[[i]]
    if(stan_data$P_partial != 0)
      stan_data$X_partial[i,,] = covariate_list_partial[[i]]
  }
  if(stan_data$P_partial == 0) # special case when ther are no partial covariates
    stan_data$P_partial = 1

  # Normalise positive covariate effects to be between 0 and 1
  stan_data$X=normalise_covariate_array(stan_data$X)
  # Normalise positive pooled covariate effects to be between 0 and 1
  stan_data$X_partial=normalise_covariate_array(stan_data$X_partial)
  
  return(list("stan_data" = stan_data, "dates" = dates, "reported_cases"=reported_cases,
   "deaths_by_country" = deaths_by_country, "processed_mobility"=processed_mobility))
}

normalise_covariate_array <- function(covariate_array){
  # Normalise positive covariate effects to be between 0 and 1
  dm=dim(covariate_array)
  for(j in 1:dm[3]){ # for covariates
    for(i in 1:dm[1]){ # countries
      raw=covariate_array[i,,j]
      if(all(raw!=0)){
        top = raw[raw>=0]
        bottom = raw[raw<=0]
        adjusted=raw
        # only act on real valued covariates (not the discrete interventions)
        if(sum(top==1)!=length(top)){
          top=rescale(top,to=c(0,1))
          adjusted[raw>=0]=top
        }
        covariate_array[i,,j] = adjusted
      }
    }
  }
  return(covariate_array)
}

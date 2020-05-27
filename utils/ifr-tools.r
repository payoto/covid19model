
source("utils/read-data.r")

return_ifr <- function (){
    return(read_ifr_data())
}

ifr_from_array <- function(ifr.by.country, Country, Region){

	if(any(ifr.by.country$country == Region)){
        # to add
      IFR <- ifr.by.country$ifr[ifr.by.country$country == Region]
    } else if(any(ifr.by.country$country == Country)) {
      IFR <- ifr.by.country$ifr[ifr.by.country$country == Country]
    } else {
    	stop(message(sprintf(
    		"ERROR: Neither %s (region) or %s (country) found in `ifr.by.country`",
    		Region, Country
		)))
    }

    return(IFR)
}

population_from_array <- function(ifr.by.country, region_data, Country, Region){

	if(any(ifr.by.country$country == Region)){
        # to add
      region_pop <- ifr.by.country$popt[ifr.by.country$country==Region]
    } else if(any(ifr.by.country$country == Country))  {
      region_pop <- region_data$popData[1]
    } else {
    	stop(message(sprintf(
    		"ERROR: Neither %s (region) or %s (country) found in `ifr.by.country`",
    		Region, Country
		)))
    }

    return(region_pop)
}
# Other tools coming
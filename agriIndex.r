#script to create agriculture capital map  (e.g. to compare to production data)

#requires maps from soilMap.r, slopeMap.r and climate data from cru_ts4.zip (via https://crudata.uea.ac.uk/cru/data/hrg/)

#To prepare this map of land use capacity we used two key variables that influences large-scale annual crop production in Brazil.
#Variable 1: Slope
#Variable 2: Climate (using method described in Victoria et al. 2007 DOI: 10.1175/EI198.1)
#These are combined into a final Agriculture Capital. 


rm(list = ls())
library(raster)
library(tidyverse)
library(ncdf4)


#read munis.r as latlong
#unzip(zipfile="Data/sim10_BRmunis_latlon_5km_2018-04-27.zip",exdir="Data")  #unzip
munis.r <- raster("SpatialData/sim10_BRmunis_latlon_5km_2018-04-27.asc")
latlong <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs "
crs(munis.r) <- latlong

#Variable 1: Slope
#In Brazil up to 90% of the agricultural areas of annual crops, such as soybean and maize, are in slopes less than 5%. 
#According to Pereira and Lombardi (2004), mechanization limitation due to slope, in an area without rocks or stones issues, are as follows:
#	0 – (0 – 3% No limitation)
#	1 – (3 – 8% Slight limitation )
#	2 – (8 – 13% Moderate limitation)
#	3 – (13 – 20% Strong limitation)
#	4* – (20 – 45 % Very strong limitation)
#	4* – (>45% Unsuitable)
#NB: *We use the code 4 for both slope classes because for the purposes of soybean/maize production, both are unsuitable.

#the slope map with the above classes is created in slopeMAp.r

#read soil texture and set bucket size (see email from Daniel Victoria 2017-11-21)
#unzip(zipfile="Data/vslope_2018-04-30.zip",exdir="Data")  #unzip
vslope.m<-raster("SpatialData/vslope_2018-04-30.asc")  



#Variable 2: Climate
#Climate has an importance in all stages of the plant growth, from the planting season to harvest. 
#To understand climatic limitations associated with plant growth and agricultural production, we used the dryness index, which describe the relation between water deficit and potential evapotranspiration (Pereira and Lombardi, 2004), both obtained from the Thornthwaite – Matter climatic water balance. 
#This index express the water deficit in percentage of potential evapotranspiration and is calculated by the equation:
#Di = 100 DEF / PET
#where Di (%) is the dryness index; DEF represents the water deficit; and PET the potential evapotranspiration. 
#In addition, the dryness index was  combined with the number of months where the  water deficit was greater than 5 mm (Pereira and Lombardi, 2004). 

#To cultivate at least one crop year-round such as soybean or maize, the region of production needs around 4 to 5 months of rainfall and the distribution and amount of available water is key factor. 
#Therefore, the proposed method use the Table 1 to calculate per pixel, the combination of Di with the length of the dry season. 
#For these calculations, the Worldclim monthly climatic dataset (30 year monthly mean) was used.


#soil is used to calculate plant available water
#read soil texture and set bucket size (see email from Daniel Victoria 2017-11-21)
#unzip(zipfile="Data/soilT_2018-05-01.zip",exdir="Data")  #unzip
soil<-raster("SpatialData/soilT_2018-05-01.asc")  

#PAW is Plant Available Water
PAW<-soil
PAW[PAW==1]<-0
PAW[PAW==2]<-0
PAW[PAW==3]<-75
PAW[PAW==4]<-55
PAW[PAW==5]<-35
#plot(PAW)


#unzip cruts climate data files (if needed)
#unzip(zipfile="Data/cru_ts4.zip",exdir="Data/cruts")  # unzip all files 

#this function January to December (NH agricultural year)
nc2raster <- function(ncyear, ncvar)
{
  #this is hacked version of cruts2raster function in library(cruts) see https://rdrr.io/cran/cruts/src/R/import-export.R
  #returns a raster stack of 12 layers for a single year from cru_ts data (e.g. http://data.ceda.ac.uk//badc/cru/data/cru_ts/cru_ts_4.01/data/)
  #needed becase cruts2raster returns same data (year) regardless of timeRange provided 
  
  #ncname is a character string of the ncdf file (e.g."cru_ts4.01.2001.2010.pre.dat.nc")
  #ncyear is an integer indicating which year from the ncdf file is to be returned (so for above file ncyear <- 1 would return data for 2001, ncyear <- 4 would give 2004 etc)
  #ncvar is a character string indicating which variable from the nc file to access (e.g. "pre", "tmp")
  
  ncname <- fn_fromYearVar(ncyear,ncvar)
  nc <- nc_open(ncname)
  pre_array <- ncvar_get(nc,ncvar)
  lon <- ncvar_get(nc,"lon")
  lat <- ncvar_get(nc,"lat")
  
  M <- length(lon)
  N <- length(lat)
  dx <- diff(lon[1:2])
  dy <- diff(lat[1:2])
  
  tr <- y_fromYear(ncyear)
  startmonth <- (tr * 12) - 11
  
  s1 <- raster(t(pre_array[,,startmonth][,N:1]), xmn=lon[1]-dx/2, xmx=lon[M]+dx/2, ymn=lat[1]-dy/2, ymx=lat[N]+dy/2, crs=CRS("+init=epsg:4326"))
  
  startmonth <- startmonth + 1
  endmonth <- startmonth + 10 
  
  for(mon in startmonth:endmonth)
  {
    s1 <- stack(s1, raster(t(pre_array[,,mon][,N:1]), xmn=lon[1]-dx/2, xmx=lon[M]+dx/2, ymn=lat[1]-dy/2, ymx=lat[N]+dy/2, crs=CRS("+init=epsg:4326")))
  }
  
  names(s1) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  
  return(s1)
}



#this function July to June (SH agricultural year)
nc2rasterSH <- function(ncyear, ncvar)
{
  #this is hacked version of cruts2raster function in library(cruts) see https://rdrr.io/cran/cruts/src/R/import-export.R
  #returns a raster stack of 12 layers for a single year from cru_ts data (e.g. http://data.ceda.ac.uk//badc/cru/data/cru_ts/cru_ts_4.01/data/)
  #needed becase cruts2raster returns same data (year) regardless of timeRange provided 
  
  #ncname is a character string of the ncdf file (e.g."cru_ts4.01.2001.2010.pre.dat.nc")
  #ncyear is an integer indicating which year from the ncdf file is to be returned (so for above file ncyear <- 1 would return data for 2001, ncyear <- 4 would give 2004 etc)
  #ncvar is a character string indicating which variable from the nc file to access (e.g. "pre", "tmp")
  

  #for SH neeed to get two years of data! 
  #agricultural year is July ncyear-1 (yrA) to June ncyear (yrB)
  
  #as long as ncyear is not the first year in a ncdf file (e.g. not 2001) we can use just a single ncdf file
  #otherwise we need to get data two ncdf file (this is handled below by the if statement)
  
  #so first get the ncdf file in which yrB is located
  
  ncnameB <- fn_fromYearVar(ncyear,ncvar)
  ncB <- nc_open(ncnameB)  #get data for yrB

  #get these parameters here (if we need to get for another ncdf file we will below)
  pre_arrayB <- ncvar_get(ncB,ncvar)
  lonB <- ncvar_get(ncB,"lon")
  latB <- ncvar_get(ncB,"lat")
  
  MB <- length(lonB)
  NB <- length(latB)
  dxB <- diff(lonB[1:2])
  dyB <- diff(latB[1:2])
  
  #need to split the remaining stacking into two loops, one for yrB, one for yrA
  #startmonth becomes July of yrA
  #endmonth of yrA is December 
  
  #startmonth of yrB is January
  #endmonth of yrB is June
  
  #startmonth <- (ncyear * 12) - 11
  
  tr <- y_fromYear(ncyear)
  
  #works as long as ncyear (yrB) is not the first in the ncdf file
  if(tr != 1) {
    
    startmonth <- (tr * 12) - 17  
    
    s1 <- raster(t(pre_arrayB[,,startmonth][,NB:1]), xmn=lonB[1]-dxB/2, xmx=lonB[MB]+dxB/2, ymn=latB[1]-dyB/2, ymx=latB[NB]+dyB/2, crs=CRS("+init=epsg:4326"))
    
    startmonth <- startmonth + 1
    endmonth <- startmonth + 10 
    
    for(mon in startmonth:endmonth)
    {
      s1 <- stack(s1, raster(t(pre_arrayB[,,mon][,NB:1]), xmn=lonB[1]-dxB/2, xmx=lonB[MB]+dxB/2, ymn=latB[1]-dyB/2, ymx=latB[NB]+dyB/2, crs=CRS("+init=epsg:4326")))
    }
  }
  
  #otherwise data for yrA is in an entirely different ncdf file
  else { 
    
    ncnameA <- fn_fromYearVar(ncyear-1, ncvar)
    print(paste0("reading file ",ncnameA))

    ncA <- nc_open(ncnameA)  #get data for this year

    pre_arrayA <- ncvar_get(ncA,ncvar)
    lonA <- ncvar_get(ncA,"lon")
    latA <- ncvar_get(ncA,"lat")
  
    MA <- length(lonA)
    N.A <- length(latA)
    dxA <- diff(lonA[1:2])
    dyA <- diff(latA[1:2])

    #from earliest nc file get last 6 months
    trA <- y_fromYear(ncyear-1)
    startmonthA <- (trA * 12) - 5 
    
    s1 <- raster(t(pre_arrayA[,,startmonthA][,N.A:1]), xmn=lonA[1]-dxA/2, xmx=lonA[MA]+dxA/2, ymn=latA[1]-dyA/2, ymx=latA[N.A]+dyA/2, crs=CRS("+init=epsg:4326"))
        
    startmonthA <- startmonthA + 1
    endmonthA <- startmonthA + 4 
    
    for(mon in startmonthA:endmonthA)
    {
      s1 <- stack(s1, raster(t(pre_arrayA[,,mon][,N.A:1]), xmn=lonA[1]-dxA/2, xmx=lonA[MA]+dxA/2, ymn=latA[1]-dyA/2, ymx=latA[N.A]+dyA/2, crs=CRS("+init=epsg:4326")))
    }

    #from later nc file get first 6 months
    startmonth <- 1  
    endmonth <- 6 
    
    for(mon in startmonth:endmonth)
    {
      s1 <- stack(s1, raster(t(pre_arrayB[,,mon][,NB:1]), xmn=lonB[1]-dxB/2, xmx=lonB[MB]+dxB/2, ymn=latB[1]-dyB/2, ymx=latB[NB]+dyB/2, crs=CRS("+init=epsg:4326")))
    }
    
  }
  
  names(s1) <- c(paste0("Jul",ncyear-1),paste0("Aug",ncyear-1),paste0("Sep",ncyear-1),paste0("Oct",ncyear-1),paste0("Nov",ncyear-1),paste0("Dec",ncyear-1),paste0("Jan",ncyear),paste0("Feb",ncyear),paste0("Mar",ncyear),paste0("Apr",ncyear),paste0("May",ncyear),paste0("Jun",ncyear))
  
  return(s1)
}



#function to return final digit from a four-digit integer (for use in nc2raster)
y_fromYear <- function(year)
{
  y <- year %% 10
  if(y == 0) { y <- 10 }
  return(y)
}


#function to set climate file name (for use in nc2raster) from a year integer and variable (pre, tmn, tmx)
fn_fromYearVar <- function(year, var)
{
  yr = "SpatialData/cruts/cru_ts4.01.1991.2000."
  if(year > 2000 & year <= 2010) { yr = "SpatialData/cruts/cru_ts4.01.2001.2010."}
  if(year > 2010 & year <= 2016) { yr = "SpatialData/cruts/cru_ts4.01.2011.2016."}
  
  return(paste0(yr,var,".dat.nc"))
}


#extent object to use in pdf map plots
BRA.ext <- extent(-62.39713, -35.43949, -33.89756, -4.06125)


calcAgriMaps <- function(munis.r, PAW, year, BRA.e, hemi, season)
{
  #generate timeRange and filenames for this year
  tr <- y_fromYear(year)
  prefn <- fn_fromYearVar(year,"pre") #precipitation,	millimetres per month  see: https://crudata.uea.ac.uk/cru/data/hrg/#info 
  tmnfn <- fn_fromYearVar(year,"tmn") #monthly average daily minimum temperature,	degrees Celsiusunits are
  tmxfn <- fn_fromYearVar(year,"tmx") #monthly average daily maximum temperature,	degrees Celsius
  
  #month labels for plots
  monlab <- c(paste0("Jan",year),paste0("Feb",year),paste0("Mar",year),paste0("Apr",year),paste0("May",year),paste0("Jun",year),paste0("Jul",year),paste0("Aug",year),paste0("Sep",year),paste0("Oct",year),paste0("Nov",year),paste0("Dec",year))
  if(hemi == "S")  monlab <- c(paste0("Jul",year-1),paste0("Aug",year-1),paste0("Sep",year-1),paste0("Oct",year-1),paste0("Nov",year-1),paste0("Dec",year-1),paste0("Jan",year),paste0("Feb",year),paste0("Mar",year),paste0("Apr",year),paste0("May",year),paste0("Jun",year))
  
  #create season filelabel
  season_label <- paste0(season, collapse="")
  if(hemi == "S" && season_label == "JulAugSepOctNovDecJanFebMarAprMayJun") season_label <- "All"
  if(hemi == "N" && season_label == "JanFebMarAprMayJunJulAugSepOctNovDec") season_label <- "All"
  
  
  #for testing
  #print(tr)
  #print(prefn)
  #print(tmnfn)
  #print(tmxfn)
  
  #read climate files
  
  #northern hemisphere
  if(hemi == "N") {
    pre <- nc2raster(ncyear=year,ncvar="pre")
    tmn <- nc2raster(ncyear=year,ncvar="tmn")
    tmx <- nc2raster(ncyear=year,ncvar="tmx")
  }
  
  #souther hemisphere
  if(hemi == "S") {
    pre <- nc2rasterSH(ncyear=year,ncvar="pre")
    tmn <- nc2rasterSH(ncyear=year,ncvar="tmn")
    tmx <- nc2rasterSH(ncyear=year,ncvar="tmx")
  }
  
  #pdf(paste0("Data/pre",year,".pdf"))
  #plot(pre, ext = BRA.e)
  #dev.off()
  
  #project and crop bricks to extent we want for Brazil
  pre.b <- projectRaster(pre, munis.r)
  pre.b <- mask(x=pre.b, mask=munis.r)
  
  tmn.b <- projectRaster(tmn, munis.r)
  tmn.b <- mask(x=tmn.b, mask=munis.r)
  
  tmx.b <- projectRaster(tmx, munis.r)
  tmx.b <- mask(x=tmx.b, mask=munis.r)
  
  #caclulate average temperature by month (brick)
  avtemp.b <- 0.36*(3*tmx.b-tmn.b)
  
  #annual temp raster layer
  Ta<-mean(avtemp.b)
  
  #total pptn raster layer
  Pa<-sum(pre.b)
  
  #set params needed to calculate PET
  Days <- list(31,28,31,30,31,30,31,31,30,31,30,31) #list of days in the month 
  Idex <- 12*((0.2*Ta)^1.514)#1st regional thermal index
  Adex <- 0.49239+1.7912*10^-2*Idex-7.71*10^-5*Idex^2+6.75*10^-7*Idex^3#2nd regional thermal index
  Ndex <- 12##NEED REAL VALUE
  
  
  #initialize PET with mean temperatures
  PET.b <- avtemp.b 
  
  #function to calculate Potential Evapotranspiration (PET)
  calcPET <- function(PET, D, I, a, N)
  {
    #if mean temperature <= 0
    PET[PET<=0]<-0
    
    #if mean temperature > 26.5
    PET[PET>26.5]<-(-415.85+32.24*PET[PET>26.5]-0.43*(PET[PET>26.5]^2))*(N/12)*(D/30)
    
    #else
    PET[PET>0&PET<=26.5]<-0.0444*((10*(PET[PET>0&PET<=26.5]/I[PET[]>0&PET[]<=26.5]))^a[PET[]>0&PET[]<=26.5])*N*D
    
    return(PET)
  }
  
  
  #map2 to loop over raster brick (layer per month) and Days list (from purrr, see http://r4ds.had.co.nz/iteration.html)
  #brick needs to passed as a list (see https://geocompr.robinlovelace.net/location.html)
  PET.b <- 
    map2(as.list(PET.b), Days, calcPET, I = Idex, a = Adex, N = Ndex) %>% 
    stack()  #remember to re-stack the list after function
  
  names(PET.b) <- monlab
 
  #see Victoria et al. 2007 DOI: 10.1175/EI198.1 Table 2 for equations
  #initialise water storage variables 
  
  Stoi <- PAW  #Stoii is month i-1 storage
  Stoii <- Stoi #Stoi is month i storage (in first month use same values)

  allmeanStoi <- vector("double", 12)  #vector to hold meanStoi for each month
  
  #for creating empty rasters and bricks
  nullRaster <- munis.r
  nullRaster[!is.na(nullRaster)] <- 0  #set anywhere that is not 'no data' in munis.r to 0

  DEF.b <- stack(replicate(12, nullRaster)) #empty brick to save all month's DEF
  ET.b <- stack(replicate(12, nullRaster))  #empty brick to save all month's ET
  
  DEF <- nullRaster #empty layer for temp useage in loop
  ET <- nullRaster #empty layer for temp useage in loop
  
  #par(mfrow=c(1,1))
  
  #see loopProofs (need to use loop, cannot use map)
  for(i in 1:12)
  {
    #hold current values of Stoi to set Stoii for next month below (this is why we can't use map)
    tempStoi <- Stoi
    
    P <- pre.b[[i]]    #get this month's precipitation (for clarity in equations below)
    PET <- PET.b[[i]]  #get this month's PET (for clarity in equations below)
    
    #if pptn < PET set storage
    Stoi[P<PET] <- Stoii[P<PET] * exp(P[P<PET] - PET[P<PET]/PAW[P<PET])
    
    #if pptn >= PET set storage
    Stoi[P>=PET] <- Stoii[P>=PET] + (P[P>=PET] - PET[P>=PET])
    
    #update Stoii ready for next month
    Stoii<-tempStoi
    
    #where Sto > PAW
    Stoi[Stoi[]>PAW[]] <- PAW[Stoi[]>PAW[]]
    
    #save mean Stoi value for this month
    allmeanStoi[i] <- cellStats(Stoi, "mean")
    
    #calculate delta storage
    trSto <- Stoi - Stoii

    #reset ET for this loop
    ET <- nullRaster
    
    #where pptn < PET
    ET[P<PET] <- P[P<PET] - trSto[P<PET]
    
    #where P >= PET
    ET[P>=PET] <- PET[P>=PET]
    
    #reset DEF for this loop 
    DEF <- nullRaster
    
    #where pptn < PET
    DEF[P<PET] <- PET[P<PET] - ET[P<PET]
    
    #where P >= PET
    DEF[P>=PET]<-0
    
    #copy DEF to DEF brick
    DEF.b[[i]] <- DEF
    ET.b[[i]] <- ET
  
  }
  
  names(DEF.b) <- monlab    #apply month-year names to the layes (for plotting and access below)
  names(ET.b) <- monlab

  si_yr <- lapply(season, paste0, year)  #add the year to month names to match monlab format
  season_indices <- ifelse(monlab %in% si_yr,1,2)  #create index of months to use in stackApply below (1 is in season, 2 is not)
  
  #calculate Dryness Index
  #avDEF<-mean(DEF.b)#mean annual DEF
  avDEF <- stackApply(DEF.b, season_indices, mean)  #mean DEF (for specified months), #creates a stack of two layers (season months and non-season months)
  #avPET<-mean(PET.b)#mean annual PET
  avPET <- stackApply(PET.b, season_indices, mean)  #mean PET (for specified months), #creates a stack of two layers (season months and non-season months)

  avDi <- (100*avDEF) / avPET  #creates a stack of two layers (season months and non-season months)
  Di <- (100*DEF.b) / PET.b
  
  
  #pptn and temp by season if needed
  avPptn <- stackApply(pre.b, season_indices, mean)
  avTemp <- stackApply(avtemp.b, season_indices, mean)
  
  
  #Number of months with water deficit - helper function
  countWD <- function(vect, na.rm=T) { return(sum(vect > 5)) }
  
  #Number of months in the season - helper function
  countMonths <- function(vect, na.rm=T) { return(length(vect)) }

  #calculate various ways of defining water deficif months
  DEFmonths <- stackApply(DEF.b,season_indices, countWD)  #creates a stack of two layers (season months and non-season months)
  allmonths <- stackApply(DEF.b,season_indices, countMonths) #creates a stack of two layers (season months and non-season months)
  DEFmonths_prop <- DEFmonths / allmonths   #proportion,  #creates a stack of two layers (season months and non-season months)
  
  #Stoidiffc <- allmeanStoi[12] - allmeanStoi[1]  #this does not seem to be used elsewhere...
  #Stoidiffc
  
  #write data to files
  if(writeClimRast)
  {
    writeRaster(avDEF[["index_1"]], paste0(outputDir,"/",className,"/MeanDEF_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avPET[["index_1"]], paste0(outputDir,"/",className,"/MeanPET_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avTemp[["index_1"]], paste0(outputDir,"/",className,"/MeanTemp_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avPptn[["index_1"]], paste0(outputDir,"/",className,"/MeanPrecip_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avDi[["index_1"]], paste0(outputDir,"/",className,"/MeanDI_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(DEFmonths[["index_1"]], paste0(outputDir,"/",className,"/CountDEFmonths_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(DEFmonths_prop[["index_1"]], paste0(outputDir,"/",className,"/PropDEFmonths_",season_label,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
  }
  
  #write pdfs
  if(writeClimPdf)
  {

    pdf(paste0(outputDir,"/",className,"/DEF_",year,hemi,".pdf"))
    plot(DEF.b, ext = BRA.e)
    dev.off()

    pdf(paste0(outputDir,"/",className,"/PET_",year,hemi,".pdf"))
    plot(PET.b, ext = BRA.e)
    dev.off()

    pdf(paste0(outputDir,"/",className,"/ET_",year,hemi,".pdf"))
    plot(ET.b, ext = BRA.e)
    dev.off()

    pdf(paste0(outputDir,"/",className,"/PPTN_",year,hemi,".pdf"))
    plot(pre.b, ext = BRA.e)
    dev.off()

    pdf(paste0(outputDir,"/",className,"/DI_",year,hemi,".pdf"))
    plot(Di, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/ClimateVariables_",season_label,"_",year,hemi,".pdf"))
    #pdf(paste0(outputDir,"/",className,"/meanDEF_",season_label,"_",year,hemi,".pdf"))
    plot(avDEF[["index_1"]], ext = BRA.e, main=paste("meanDEF",season_label,year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()

    #pdf(paste0(outputDir,"/",className,"/meanPET_",season_label,"_",year,hemi,".pdf"))
    plot(avPET[["index_1"]], ext = BRA.e, main=paste("meanPET",season_label,year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()

    #pdf(paste0(outputDir,"/",className,"/meanDI_",season_label,"_",year,hemi,".pdf"))
    plot(avDi[["index_1"]], ext = BRA.e, main=paste("meanDI",season_label,year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/DEFmonths_",season_label,"_",year,hemi,".pdf"))
    plot(DEFmonths[["index_1"]], ext = BRA.e, main=paste("count DEFmonths",season_label,year,hemi, sep=" "))
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/DEFmonths_prop_",season_label,"_",year,hemi,".pdf"))
    plot(DEFmonths_prop[["index_1"]], ext = BRA.e, main=paste("prop DEFmonths",season_label,year,hemi, sep=" "))
    dev.off()
  
  }
  
  
  #####
  #2019-05-21
  #Switching focus to seasonal analysis means we can ignore ClimLim for now.  

  
  #Classify climate limitation (combination of deficit months and dryness index)
  #Limitation degree:
    #0 – Null, no limitation climatic restrictions
    #1 – Slight 
    #2 – Moderate
    #3 – Strong
    #4 – Very strong
  
  #ClimLim <- nullRaster

  
   
  ###new classification
  # ClimLim[DEFmonths[]<=2]<-0
  # ClimLim[DEFmonths[]>2&DEFmonths[]<=4]<-1
  # ClimLim[DEFmonths[]>4&DEFmonths[]<=6]<-2
  # ClimLim[DEFmonths[]==7]<-3
  # ClimLim[DEFmonths[]>7]<-4
  ###end new classification

  
    
  ###original classification
  # ClimLim[DEFmonths[]<1&avDi[]<20]<-0
  # ClimLim[DEFmonths[]<1&avDi[]>=20&avDi[]<40]<-0
  # ClimLim[DEFmonths[]<1&avDi[]>=40&avDi[]<60]<-1
  # ClimLim[DEFmonths[]<1&avDi[]>=60&avDi[]<80]<-2
  # ClimLim[DEFmonths[]<1&avDi[]>=80]<-3
  # 
  # ClimLim[DEFmonths[]>=1&DEFmonths[]<=3&avDi[]<20]<-0
  # ClimLim[DEFmonths[]>=1&DEFmonths[]<=3&avDi[]>=20&avDi[]<40]<-1
  # ClimLim[DEFmonths[]>=1&DEFmonths[]<=3&avDi[]>=40&avDi[]<60]<-2
  # ClimLim[DEFmonths[]>=1&DEFmonths[]<=3&avDi[]>=60&avDi[]<80]<-3
  # ClimLim[DEFmonths[]>=1&DEFmonths[]<=3&avDi[]>=80]<-4
  #  
  # ClimLim[DEFmonths[]>=4&DEFmonths[]<=6&avDi[]<20]<-1
  # ClimLim[DEFmonths[]>=4&DEFmonths[]<=6&avDi[]>=20&avDi[]<40]<-2
  # ClimLim[DEFmonths[]>=4&DEFmonths[]<=6&avDi[]>=40&avDi[]<60]<-3
  # ClimLim[DEFmonths[]>=4&DEFmonths[]<=6&avDi[]>=60&avDi[]<80]<-4
  # ClimLim[DEFmonths[]>=4&DEFmonths[]<=6&avDi[]>=80]<-4
  #  
  # ClimLim[DEFmonths[]>=7&DEFmonths[]<=9&avDi[]<20]<-2
  # ClimLim[DEFmonths[]>=7&DEFmonths[]<=9&avDi[]>=20&avDi[]<40]<-3
  # ClimLim[DEFmonths[]>=7&DEFmonths[]<=9&avDi[]>=40&avDi[]<60]<-4
  # ClimLim[DEFmonths[]>=7&DEFmonths[]<=9&avDi[]>=60&avDi[]<80]<-4
  # ClimLim[DEFmonths[]>=7&DEFmonths[]<=9&avDi[]>=80]<-4
  #  
  # ClimLim[DEFmonths[]>=10&DEFmonths[]<=12&avDi[]<20]<-3
  # ClimLim[DEFmonths[]>=10&DEFmonths[]<=12&avDi[]>=20&avDi[]<40]<-4
  # ClimLim[DEFmonths[]>=10&DEFmonths[]<=12&avDi[]>=40&avDi[]<60]<-4
  # ClimLim[DEFmonths[]>=10&DEFmonths[]<=12&avDi[]>=60&avDi[]<80]<-4
  # ClimLim[DEFmonths[]>=10&DEFmonths[]<=12&avDi[]>=80]<-4
  ###end original classification
  
  #par(mfrow=c(1,1))
  #plot(ClimLim)
#uncomment  
  # writeRaster(ClimLim, paste0(outputDir,"/",className,"/ClimLim",year,".asc"), format = 'ascii', overwrite=T)
  # 
  # pdf(paste0(outputDir,"/",className,"/ClimLim",year,".pdf"))
  # plot(ClimLim, ext = BRA.e)
  # dev.off()
  # 

  # 
  # #create map of double cropping
  # #discussion with Ramon by email in March 2019 suggested that DC is not possible on limitation 3 and 4
  # dc.f <- ClimLim
  # values(dc.f)[values(dc.f)==0 | values(dc.f)==2 | values(dc.f)==3 ] = 1
  # values(dc.f)[values(dc.f)==4] = 0
  # #writeRaster(dc.f, paste0(outputDir,"/",className,"/DC",year,".asc"), format = 'ascii', overwrite=T)
  # 
  # pdf(paste0(outputDir,"/",className,"/DC",year,".pdf"))
  # plot(dc.f, ext = BRA.e)
  # dev.off()

  
  #Considering the two variables, we join by combining the 6 classes of slope with the 5 classes of climate, resulting in 30 combinations grouped according to the table (below).
  #Each combination is reclassified with a code accordingly to the maximum value of the combination (the value of one of the two variables at least). 
  #The code (suitability classes which equate to Agricultural Capital for CRAFTY) for the combined map was divided by 4 to represent the scale values from 0 to 1.

#uncomment
  # vagri.f <- ClimLim
  # values(vagri.f)[values(ClimLim)==0&values(vslope.m)==0] = 1      #best
  # values(vagri.f)[values(ClimLim)==0&values(vslope.m)==1] = 0.75
  # values(vagri.f)[values(ClimLim)==0&values(vslope.m)==2] = 0.5
  # values(vagri.f)[values(ClimLim)==0&values(vslope.m)==3] = 0.25
  # values(vagri.f)[values(ClimLim)==0&values(vslope.m)==4] = 0.1
  # 
  # values(vagri.f)[values(ClimLim)==1&values(vslope.m)==0] = 0.75
  # values(vagri.f)[values(ClimLim)==1&values(vslope.m)==1] = 0.75
  # values(vagri.f)[values(ClimLim)==1&values(vslope.m)==2] = 0.5
  # values(vagri.f)[values(ClimLim)==1&values(vslope.m)==3] = 0.25
  # values(vagri.f)[values(ClimLim)==1&values(vslope.m)==4] = 0.1
  # 
  # values(vagri.f)[values(ClimLim)==2&values(vslope.m)==0] = 0.5
  # values(vagri.f)[values(ClimLim)==2&values(vslope.m)==1] = 0.5
  # values(vagri.f)[values(ClimLim)==2&values(vslope.m)==2] = 0.5
  # values(vagri.f)[values(ClimLim)==2&values(vslope.m)==3] = 0.25
  # values(vagri.f)[values(ClimLim)==2&values(vslope.m)==4] = 0.1
  # 
  # values(vagri.f)[values(ClimLim)==3&values(vslope.m)==0] = 0.25
  # values(vagri.f)[values(ClimLim)==3&values(vslope.m)==1] = 0.25
  # values(vagri.f)[values(ClimLim)==3&values(vslope.m)==2] = 0.25
  # values(vagri.f)[values(ClimLim)==3&values(vslope.m)==3] = 0.25
  # values(vagri.f)[values(ClimLim)==3&values(vslope.m)==4] = 0.1
  # 
  # values(vagri.f)[values(ClimLim)==4&values(vslope.m)==0] = 0.1
  # values(vagri.f)[values(ClimLim)==4&values(vslope.m)==1] = 0.1
  # values(vagri.f)[values(ClimLim)==4&values(vslope.m)==2] = 0.1
  # values(vagri.f)[values(ClimLim)==4&values(vslope.m)==3] = 0.1
  # values(vagri.f)[values(ClimLim)==4&values(vslope.m)==4] = 0.1
  # values(vagri.f)[values(vagri.f)>1]=0                              #worst

  #!check does this need to be resampled/masked before writing?!
  # writeRaster(vagri.f, paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"), format = 'ascii', overwrite=T)
  # 
  # pdf(paste0(outputDir,"/",className,"/agricultureCapital",year,".pdf"))
  # plot(vagri.f, ext = BRA.e)
  # dev.off()
  
  rm(pre,tmn,tmx,pre.b,tmn.b,tmx.b,PET.b,DEF.b,ET.b
    #ClimLim,vagri.f
    )
  
}



outputDir <- "Output"
className <- "class-A-SH"

#crop_season <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
#soy_season <- c("Oct","Nov","Dec","Jan","Feb","Mar")
mz1_season <- c("Oct","Nov","Dec","Jan","Feb","Mar")
mz2_season <- c("Jan","Feb","Mar","Apr","May","Jun")
crops <- list(mz1_season, mz2_season)


#create the output directory for this classification if it does not exist
if(!dir.exists(paste0(outputDir,"/",className))) { dir.create(paste0(outputDir,"/",className)) }


#yr <-2015
for(crop in crops) {
  for(yr in 2003:2016)
  {
    writeClimRast <- T
    writeClimPdf <- T
    calcAgriMaps(munis.r, PAW, yr, BRA.ext, "S", crop)
    print(paste0(yr," done"))
  }
}


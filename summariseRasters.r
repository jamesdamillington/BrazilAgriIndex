#summarise maps over time (at municipality level)

rm(list = ls())
library(raster)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(viridisLite)

years <- seq(2003,2016,1)

className <- "class-A-SH"
outputDir <- "Output"
munis.r <- raster("SpatialData/sim10_BRmunis_latlon_5km_2018-04-27.asc")
BRmunis <- st_read("SpatialData/BRmunis_sim10_simple2.shp")


meanMuniFromRaster <- function(varName, munis, year) {
  
  #get raster data  
  varRas <- raster(paste0(outputDir,"/",className,"/", varName, year,".asc"))
 
  #summarise raster data by municipality
  meanMtx <- zonal(varRas, munis, fun = "mean")
  
  return(meanMtx)
  
}
  
joinMeans <- function(DataTable, meanMtx, counter) {

  if(counter == 2) DataTable <- DataTable %>% rename_at(c("mean"), funs( paste0(years[counter-1]) ) ) #from https://stackoverflow.com/a/47044471/10219907
    #funs( paste0(., years[counter-1]) )
  DataTable <- DataTable %>% 
    left_join(., as_tibble(meanMtx), by = c("CD_GEOCMUn" ="zone")) %>% 
    rename_at(c("mean"), funs( paste0(years[counter]) ) )  #from https://stackoverflow.com/a/47044471/10219907

  return(DataTable)
}


writeData <- function(data_map, filename) {
  
  st_geometry(data_map) <- NULL
  
  char <- c("NM_MUNICIP","NM_ESTADO","CD_GEOCUF","NM_REGIAO","CD_GEOCMU")
  inte <- c( "CD_GEOCMUn","State")
  
  #following to ensure correct types are written
  #from https://stackoverflow.com/a/27668324/10219907
  data_map <- data_map %>%
    mutate_at(char,as.character) %>% 
    mutate_at(inte, as.integer) %>%
    mutate_at(vars(starts_with("20")),as.double)
  
  write_csv(data_map, paste0(outputDir,"/",className,"/",filename))
  
}

for(i in seq_along(years))
{
  print(years[i])

  #summarise CL
  meanCL <- meanMuniFromRaster("ClimLim", munis.r, years[i])
  if(i == 1) { meanCL_map <- left_join(BRmunis, as_tibble(meanCL), by = c("CD_GEOCMUn" ="zone"))
  } else { meanCL_map <- joinMeans(meanCL_map, meanCL, i) }

  #summarise AC
  meanAC <- meanMuniFromRaster("agricultureCapital", munis.r, years[i])
  if(i == 1) { meanAC_map <- left_join(BRmunis, as_tibble(meanAC), by = c("CD_GEOCMUn" ="zone"))
  } else { meanAC_map <- joinMeans(meanAC_map, meanAC, i) }

  #summarise DEFcount
  DEFc <- meanMuniFromRaster("CountDEFmonths", munis.r, years[i])
  if(i == 1) { DEFc_map <- left_join(BRmunis, as_tibble(DEFc), by = c("CD_GEOCMUn" ="zone"))
  } else { DEFc_map <- joinMeans(DEFc_map, DEFc, i) }

  #summarise DEF mean
  DEFm <- meanMuniFromRaster("MeanAnnualDEF", munis.r, years[i])
  if(i == 1) { DEFm_map <- left_join(BRmunis, as_tibble(DEFm), by = c("CD_GEOCMUn" ="zone"))
  } else { DEFm_map <- joinMeans(DEFm_map, DEFm, i) }

  #summarise DryIdx
  DryIdx <- meanMuniFromRaster("MeanAnnualDI", munis.r, years[i])
  if(i == 1) { DryIdx_map <- left_join(BRmunis, as_tibble(DryIdx), by = c("CD_GEOCMUn" ="zone"))
  } else { DryIdx_map <- joinMeans(DryIdx_map, DryIdx, i) }

  
}
  
writeData(meanCL_map, "/ClimLim_MuniMeans.csv")
writeData(meanAC_map, "/AgriCap_MuniMeans.csv")
writeData(DEFc_map, "/countDEFmonths_MuniMeans.csv")
writeData(DEFm_map, "/meanAnnualDEF_MuniMeans.csv")
writeData(DryIdx_map, "/MeanAnnualDI_MuniMeans.csv")


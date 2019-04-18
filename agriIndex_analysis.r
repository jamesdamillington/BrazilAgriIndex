#script to compare agriculture capital (maps) to municipality production data)

rm(list = ls())
library(raster)
library(tidyverse)



munis.r <- raster("SpatialData/sim10_BRmunis_latlon_5km_2018-04-27.asc")

#soy
soy <- readr::read_csv("ProductionData/soybean_Crop.csv",col_types=cols(.default=col_integer()), na="NA")

for(year in 2011:2015)
{
  #year <- 2011
  className <- "class-A"
  outputDir <- "Output"
  
  ClimLim <- raster(paste0(outputDir,"/",className,"/ClimLim",year,".asc"))
  AgriCap <- raster(paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"))
  
  meanCL <- zonal(ClimLim, munis.r, fun = 'mean')
  meanAC <- zonal(AgriCap, munis.r, fun = 'mean')
  
  soy_yr <- soy %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), soy_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), soy_yr, by = c("zone" = "CODE"))
  
  print(year)
  print("ClimLim Soy")
  print(cor.test(joinCL$mean, joinCL$prod, method="kendall"))
  print("AgriCap Soy")
  print(cor.test(joinAC$mean, joinAC$prod, method="kendall"))
  
  plot(joinCL$mean, joinCL$prod, main=paste0("ClimLim soy",year))
  plot(joinAC$mean, joinAC$prod, main=paste0("AgriCap soy",year))
}


#maize1
maize1 <- readr::read_csv("ProductionData/maize_first.csv",col_types=cols(.default=col_integer()), na="NA")
for(year in 2011:2015)
{
  #year <- 2011
  className <- "class-A"
  outputDir <- "Output"
  
  ClimLim <- raster(paste0(outputDir,"/",className,"/ClimLim",year,".asc"))
  AgriCap <- raster(paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"))
  
  meanCL <- zonal(ClimLim, munis.r, fun = 'mean')
  meanAC <- zonal(AgriCap, munis.r, fun = 'mean')
  
  maize1_yr <- maize1 %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), maize1_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), maize1_yr, by = c("zone" = "CODE"))
  
  print(year)
  print("ClimLim maize1")
  print(cor.test(joinCL$mean, joinCL$prod, method="kendall"))
  print("AgriCap maize1")
  print(cor.test(joinAC$mean, joinAC$prod, method="kendall"))
  
  plot(joinCL$mean, joinCL$prod, main=paste0("ClimLim maize1 ",year))
  plot(joinAC$mean, joinAC$prod, main=paste0("AgriCap maize1 ",year))
}



#maize2
maize2 <- readr::read_csv("ProductionData/maize_second.csv",col_types=cols(.default=col_integer()), na="NA")
for(year in 2011:2015)
{
  #year <- 2011
  className <- "class-A"
  outputDir <- "Output"
  
  ClimLim <- raster(paste0(outputDir,"/",className,"/ClimLim",year,".asc"))
  AgriCap <- raster(paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"))
  
  meanCL <- zonal(ClimLim, munis.r, fun = 'mean')
  meanAC <- zonal(AgriCap, munis.r, fun = 'mean')
  
  maize2_yr <- maize2 %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), maize2_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), maize2_yr, by = c("zone" = "CODE"))
  
  print(year)
  print("ClimLim maize2")
  print(cor.test(joinCL$mean, joinCL$prod, method="kendall"))
  print("AgriCap maize2")
  print(cor.test(joinAC$mean, joinAC$prod, method="kendall"))
  
  plot(joinCL$mean, joinCL$prod, main=paste0("ClimLim maize2 ",year))
  plot(joinAC$mean, joinAC$prod, main=paste0("AgriCap maize2 ",year))
}

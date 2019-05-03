#script to compare agriculture capital (maps) to municipality production data)

rm(list = ls())
library(raster)
library(tidyverse)
library(sf)
library(RColorBrewer)
library(viridisLite)

years <- seq(2003,2016,1)

className <- "Original-NH"
outputDir <- "Output"
munis.r <- raster("SpatialData/sim10_BRmunis_latlon_5km_2018-04-27.asc")
BRmunis <- st_read("SpatialData/BRmunis_sim10_simple2.shp")

soy <- readr::read_csv("ProductionData/soybean_Crop.csv",col_types=cols(.default=col_integer()), na="NA")
maize1 <- readr::read_csv("ProductionData/maize_first.csv",col_types=cols(.default=col_integer()), na="NA")
maize2 <- readr::read_csv("ProductionData/maize_second.csv",col_types=cols(.default=col_integer()), na="NA")

#containers to hold summary stats
soyCL <- c()
soyAC <- c()
soyCLp <- c()
soyACp <- c()
mz1CL <- c()
mz1AC <- c()
mz1CLp <- c()
mz1ACp <- c()
mz2CL <- c()
mz2AC <- c()
mz2CLp <- c()
mz2ACp <- c()

#year <- 2015

for(year in years)
{
  print(year)

  ClimLim <- raster(paste0(outputDir,"/",className,"/ClimLim",year,".asc"))
  AgriCap <- raster(paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"))
  
  meanCL <- zonal(ClimLim, munis.r, fun = 'mean')
  meanAC <- zonal(AgriCap, munis.r, fun = 'mean')
  
  soy_yr <- soy %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), soy_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), soy_yr, by = c("zone" = "CODE"))
  
  cl.test <- cor.test(joinCL$mean, joinCL$prod, method="kendall")
  ac.test <- cor.test(joinAC$mean, joinAC$prod, method="kendall")

  soyCL <- c(soyCL, cl.test$estimate)
  soyAC <- c(soyAC, ac.test$estimate)
  soyCLp <- c(soyCLp, cl.test$p.value)
  soyACp <- c(soyACp, ac.test$p.value)
  
  
  soyCLs <- ggplot(joinCL, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Soy ",year)) +
    xlab("Climate Limitation") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinCL$mean, na.rm=T),y=max(joinCL$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(cl.test$estimate,3),", p=",round(cl.test$p.value,3)))
  
  
  soyACs <- ggplot(joinAC, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Soy ",year)) +
    xlab("Agriculture Capital") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinAC$mean, na.rm=T),y=max(joinAC$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(ac.test$estimate,3),", p=",round(ac.test$p.value,3)))


  #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)
  
  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone")) 
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone")) 
  
  pdf(paste0(outputDir,"/",className,"/",year, "_ClimateAndCapitalMaps.pdf"))
  
  #CL plot
  pbrks = seq(from=0,to=max(CL_map$mean, na.rm=T),length.out=101)
  plot(CL_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Climate Limitation ",year))  

  #AC plot
  pbrks = seq(from=0,to=max(AC_map$mean, na.rm=T),length.out=101)
  plot(AC_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Agriculture Capital ",year))  

  dev.off()
  
  pdf(paste0(outputDir,"/",className,"/",year,"_SoyPlots.pdf"))
  
  print(soyCLs)
  print(soyACs)
  
  #production plot
  pbrks = seq(from=0,to=max(CL_map$prod, na.rm=T),length.out=101)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy Productivity ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy Productivity ",year), add = T)
  

  #CL residual plot
  coolwarm_hcl <- colorspace::diverge_hcl(8,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
  brks <- seq(from=-2000,to=2000,by=500)  #9 values
  
  #CL residual plot
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy vs Climate Limitation Residuals ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy vs Climate Limitation Residuals ",year), add = T)  

  #ACresidual plot
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy vs Agri Capital Residuals ",year), reset = F)  
  plot(st_geometry(AC_map), col="lightgray", lty=0, reset = F, add=T)
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Soy vs Agri Capital Residuals ",year), add = T)  

  dev.off()

  #maize1
  maize1_yr <- maize1 %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), maize1_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), maize1_yr, by = c("zone" = "CODE"))
  
  cl.test <- cor.test(joinCL$mean, joinCL$prod, method="kendall")
  ac.test <- cor.test(joinAC$mean, joinAC$prod, method="kendall")

  mz1CL <- c(mz1CL, cl.test$estimate)
  mz1AC <- c(mz1AC, ac.test$estimate)
  mz1CLp <- c(mz1CLp, cl.test$p.value)
  mz1ACp <- c(mz1ACp, ac.test$p.value)
  
  mz1CLs <- ggplot(joinCL, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Maize First ",year)) +
    xlab("Climate Limitation") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinCL$mean, na.rm=T),y=max(joinCL$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(cl.test$estimate,3),", p=",round(cl.test$p.value,3)))


  mz1ACs <- ggplot(joinAC, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Maize First ",year)) +
    xlab("Agriculture Capital") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinAC$mean, na.rm=T),y=max(joinAC$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(ac.test$estimate,3),", p=",round(ac.test$p.value,3)))

  
    #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)
  
  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone")) 
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone")) 

  
  pdf(paste0(outputDir,"/",className,"/",year,"_Maize1Plots.pdf"))
  
  print(mz1CLs)
  print(mz1ACs)
  
  #production plot
  pbrks = seq(from=0,to=max(CL_map$prod, na.rm=T),length.out=101)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First Productivity ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First Productivity ",year), add = T)

  #palette
  coolwarm_hcl <- colorspace::diverge_hcl(10,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
  brks <- seq(from=-7500,to=7500,by=1500)  #11 values
  
  #CL residual plot
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First vs Climate Limitation Residuals ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First vs Climate Limitation Residuals ",year), add = T)  

  #ACresidual plot
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First vs Agri Capital Residuals ",year), reset = F)  
  plot(st_geometry(AC_map), col="lightgray", lty=0, reset = F, add=T)
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize First vs Agri Capital Residuals ",year), add = T)  

  dev.off()

  

  #maize2
  maize2_yr <- maize2 %>% 
    select(CODE, as.character(year)) %>%
    rename(prod = as.character(year))

  joinCL <- left_join(as_tibble(meanCL), maize2_yr, by = c("zone" = "CODE"))
  joinAC <- left_join(as_tibble(meanAC), maize2_yr, by = c("zone" = "CODE"))
  
  cl.test <- cor.test(joinCL$mean, joinCL$prod, method="kendall")
  ac.test <- cor.test(joinAC$mean, joinAC$prod, method="kendall")

  mz2CL <- c(mz2CL, cl.test$estimate)
  mz2AC <- c(mz2AC, ac.test$estimate)
  mz2CLp <- c(mz2CLp, cl.test$p.value)
  mz2ACp <- c(mz2ACp, ac.test$p.value)
  
  mz2CLs <- ggplot(joinCL, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Maize Second ",year)) +
    xlab("Climate Limitation") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinCL$mean, na.rm=T),y=max(joinCL$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(cl.test$estimate,3),", p=",round(cl.test$p.value,3)))


  mz2ACs <- ggplot(joinAC, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste0("Maize Second ",year)) +
    xlab("Agriculture Capital") +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(joinAC$mean, na.rm=T),y=max(joinAC$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(ac.test$estimate,3),", p=",round(ac.test$p.value,3)))

  
      #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)
  
  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone")) 
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone")) 

  
  pdf(paste0(outputDir,"/",className,"/",year,"_Maize2Plots.pdf"))
  
  print(mz2CLs)
  print(mz2ACs)
  
  #production plot
  pbrks = seq(from=0,to=max(CL_map$prod, na.rm=T),length.out=101)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second Productivity ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second Productivity ",year), add = T)

  #palette
  coolwarm_hcl <- colorspace::diverge_hcl(10,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
  brks <- seq(from=-7500,to=7500,by=1500)  #11 values
  
  #CL residual plot
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second vs Climate Limitation Residuals ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second vs Climate Limitation Residuals ",year), add = T)  

  #ACresidual plot
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second vs Agri Capital Residuals ",year), reset = F)  
  plot(st_geometry(AC_map), col="lightgray", lty=0, reset = F, add=T)
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Maize Second vs Agri Capital Residuals ",year), add = T)  


  dev.off()

}

dat <- cbind(years,soyCL,soyCLp,soyAC,soyACp,mz1CL,mz1CLp,mz1AC,mz1ACp,mz2CL,mz2CLp,mz2AC,mz2ACp)
write_csv(as.data.frame(dat), paste0(outputDir,"/",className,"/Correlations.csv"))

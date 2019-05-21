#########################
#  This script calculates annual correlations for between agricultural productivity 
#  and various variables:
#       - 'climate limitation'
#       - 'agricultural capital'
#       - count of months with 'water deficit'
#       - mean annual 'water deficit' (mean of monthly values)
#       - mean annual dryness index (mean of monthly values)
#  These correlations are for each of three crops (soy, maize as first  crop, maize second crop)
#  Correlations are calculated for municipalities in all 10 states, by region and by state
#########################


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

soy <- readr::read_csv("ProductionData/soybean_Crop.csv",col_types=cols(.default=col_integer()), na="NA")
maize1 <- readr::read_csv("ProductionData/maize_first.csv",col_types=cols(.default=col_integer()), na="NA")
maize2 <- readr::read_csv("ProductionData/maize_second.csv",col_types=cols(.default=col_integer()), na="NA")

#containers to hold summary stats
soyData <- data.frame()
mz1Data <- data.frame()
mz2Data <- data.frame()
soyStates <- data.frame()
mz1States <- data.frame()
mz2States <- data.frame()
soyRegions <- data.frame()
mz1Regions <- data.frame()
mz2Regions <- data.frame()


datSummary <- function(data) {
  
  #method from https://stackoverflow.com/a/41811350/10219907
  #(unclear why safely does not work as here: https://stackoverflow.com/a/42282866/10219907)
  
  analysis <- data %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(prod ~ mean, data = .)),
    cor = map(data, possibly(
      ~broom::tidy(cor.test(.$prod, .$mean, method="kendall"), 3), otherwise = data.frame())
    )
  )
  return(analysis)
}

stateSummary <- function(data) {
  
  #method from https://stackoverflow.com/a/41811350/10219907
  #(unclear why safely does not work as here: https://stackoverflow.com/a/42282866/10219907)
  
  analysis <- data %>%
  group_by(state) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(prod ~ mean, data = .)),
    cor = map(data, possibly(
      ~broom::tidy(cor.test(.$prod, .$mean, method="kendall"), 3), otherwise = data.frame())
    )
  )
  return(analysis)
}

regionSummary <- function(data) {
  
  #method from https://stackoverflow.com/a/41811350/10219907
  #(unclear why safely does not work as here: https://stackoverflow.com/a/42282866/10219907)

  analysis <- data %>%
  group_by(regionNM) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(prod ~ mean, data = .)),
    cor = map(data, possibly(
      ~broom::tidy(cor.test(.$prod, .$mean, method="kendall"), 3), otherwise = data.frame())
    )
  )
  return(analysis)
}

collateCor <- function(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, taus) {
  
  taus <- data.frame(year, 
    analysisCL$cor[[1]]$estimate, 
    analysisCL$cor[[1]]$p.value,
    analysisAC$cor[[1]]$estimate, 
    analysisAC$cor[[1]]$p.value,
    analysisDEFc$cor[[1]]$estimate, 
    analysisDEFc$cor[[1]]$p.value,
    analysisDEFm$cor[[1]]$estimate, 
    analysisDEFm$cor[[1]]$p.value,
    analysisDI$cor[[1]]$estimate, 
    analysisDI$cor[[1]]$p.value
  )

  colnames(taus) <- c("year","CLt","CLp","ACt","ACp","DEFct","DEFcp","DEFmt","DEFmp","DIt","DIp")

  return(taus)
  
}

collateStateCor <- function(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, taus) {
  
  #to access output data for a single state
  #s51 <- filter(analysis, state == 51)
  
  #correlation estimate
  #s51$cor[[1]]$estimate

  #model estimates
  #summary(s51$model[[1]])$coefficients
  
  if(length(analysisCL[[1]]) != length(analysisAC[[1]])) { print("warning") } else {
    for(i in 1:length(analysisCL$cor))
    {
      taus <- rbind(taus, 
        c(analysisCL$cor[[i]]$estimate, 
          analysisCL$cor[[i]]$p.value,
          analysisAC$cor[[i]]$estimate, 
          analysisAC$cor[[i]]$p.value,
          analysisDEFc$cor[[i]]$estimate, 
          analysisDEFc$cor[[i]]$p.value,
          analysisDEFm$cor[[i]]$estimate, 
          analysisDEFm$cor[[i]]$p.value,
          analysisDI$cor[[i]]$estimate, 
          analysisDI$cor[[i]]$p.value
          )
        )
    }
  
    taus <- cbind(analysisCL$state, rep(year,length(analysisCL$state)),taus)
    colnames(taus) <- c("state","year","CLt","CLp","ACt","ACp","DEFct","DEFcp","DEFmt","DEFmp","DIt","DIp")
  }
  
  return(taus)
}

collateRegionCor <- function(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, taus) {
  
  if(length(analysisCL[[1]]) != length(analysisAC[[1]])) { print("warning") } else {
    for(i in 1:length(analysisCL$cor))
    {
      taus <- rbind(taus, 
        c(analysisCL$cor[[i]]$estimate, 
          analysisCL$cor[[i]]$p.value,
          analysisAC$cor[[i]]$estimate, 
          analysisAC$cor[[i]]$p.value,
          analysisDEFc$cor[[i]]$estimate, 
          analysisDEFc$cor[[i]]$p.value,
          analysisDEFm$cor[[i]]$estimate, 
          analysisDEFm$cor[[i]]$p.value,
          analysisDI$cor[[i]]$estimate, 
          analysisDI$cor[[i]]$p.value
          )
        )
    }
  
    taus <- cbind(analysisCL$regionNM, rep(year,length(analysisCL$regionNM)),taus)
    colnames(taus) <- c("region","year","CLt","CLp","ACt","ACp","DEFct","DEFcp","DEFmt","DEFmp","DIt","DIp")
  }
  
  return(taus)
}

cleanData <- function(data) {
  
    cdata <- data %>% 
    dplyr::select(CODE, as.character(year)) %>%
    rename(prod = as.character(year)) %>%
    filter(!is.na(prod)) %>%
    filter(prod < quantile(prod, 0.995)) %>%
    filter(prod > quantile(prod, 0.005)) %>%
    mutate(state = (CODE %/% 100000)) %>%
    mutate(stateNM = if_else(state == 17, "TO", 
       if_else(state == 29, "BA",
       if_else(state == 31, "MG",
       if_else(state == 35, "SP",
       if_else(state == 41, "PR",
       if_else(state == 42, "SC",
       if_else(state == 43, "RS", 
       if_else(state == 50, "MS",
       if_else(state == 51, "MT",
       if_else(state == 52, "GO", "NA"
       ))))))))))
     ) %>%
    
    #regions as discussed with Ramon 2019-05-07
    mutate(regionNM = if_else(state == 17, "N", 
       if_else(state == 29, "N",
       if_else(state == 31, "E",
       if_else(state == 35, "E",
       if_else(state == 41, "S",
       if_else(state == 42, "S",
       if_else(state == 43, "S", 
       if_else(state == 50, "W",
       if_else(state == 51, "W",
       if_else(state == 52, "W", "NA"
       ))))))))))
     ) 
    
    return(cdata)
  
}

scatterPlot <- function(dat, tau, p, commodity, year, variable) {
  
  myplot <- ggplot(dat, aes(x=mean, y=prod)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    ggtitle(paste(commodity,year)) +
    xlab(variable) +
    ylab("Productivity (kg/ha)") +
    annotate("text",x=max(dat$mean, na.rm=T),y=max(dat$prod, na.rm=T),hjust=1,vjust=0,label=paste0("t=",round(tau,3),", p=",round(p,3)))
    
    return(myplot)
}

pdfPlots <- function(pdfpath, variable, scatters, CL_map, AC_map){
    
  pdf(pdfpath)

  print(scatters)
  
  #production plot
  pbrks = seq(from=0,to=max(CL_map$prod, na.rm=T),length.out=101)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," Productivity ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map["prod"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," Productivity ",year), add = T)

  #palette
  coolwarm_hcl <- colorspace::diverge_hcl(10,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
  brks <- seq(from=-7500,to=7500,by=1500)  #11 values
  
  #CL residual plot
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Climate Limitation Residuals ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Climate Limitation Residuals ",year), add = T)  

  #ACresidual plot
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Agri Capital Residuals ",year), reset = F)  
  plot(st_geometry(AC_map), col="lightgray", lty=0, reset = F, add=T)
  plot(AC_map[".resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Agri Capital Residuals ",year), add = T)  

  #standardised residual plots (different breaks)
  coolwarm_hcl <- colorspace::diverge_hcl(10,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
  brks <- seq(from=-5,to=5,by=1)  #11 values

  #CL standardised residual plot
  plot(CL_map[".std.resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Climate Limitation, Standardised Residuals ",year), reset = F)  
  plot(st_geometry(CL_map), col="lightgray", lty=0, reset = F, add=T)
  plot(CL_map[".std.resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(CL_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Climate Limitation, Standardised Residuals ",year), add = T)  

  #ACresidual plot
  plot(AC_map[".std.resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Agri Capital, Standardised Residuals ",year), reset = F)  
  plot(st_geometry(AC_map), col="lightgray", lty=0, reset = F, add=T)
  plot(AC_map[".std.resid"], pal = coolwarm_hcl, breaks = brks, graticule = st_crs(AC_map), key.pos = 4, axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0(variable," vs Agri Capital, Standardised Residuals ",year), add = T)  

  dev.off()
    
}

stateScatter <- function(data, commodity, year, variable) {

  #scatter plot with line colouring by state
  myplot <- ggplot(data, aes(x=mean, y=prod)) +
    geom_point(aes(colour = stateNM), shape=1) +    # Use hollow circles
    geom_smooth(aes(group=stateNM, colour = stateNM),method=lm) +
    ggtitle(paste(commodity,year)) +
    xlab(variable) +
    ylab("Productivity (kg/ha)")
  
  return(myplot)
}

regionScatter <- function(data, commodity, year, variable) {
    
  #scatter plot with line colouring by region 
  myplot <- ggplot(data, aes(x=mean, y=prod)) +
    geom_point(aes(colour = regionNM), shape=1) +    # Use hollow circles
    geom_smooth(aes(group=regionNM, colour = regionNM),method=lm) +
    ggtitle(paste(commodity,year)) +
    xlab(variable) +
    ylab("Productivity (kg/ha)")
  
  return(myplot)
    
}


year <- 2016

for(year in years)
{
  print(year)

  
  #######CLIMATE AND CAPITAL DATA ###################
  #get raster data  
  ClimLim <- raster(paste0(outputDir,"/",className,"/ClimLim",year,".asc"))
  AgriCap <- raster(paste0(outputDir,"/",className,"/agricultureCapital",year,".asc"))
  DEFcount <- raster(paste0(outputDir,"/",className,"/CountDEFmonths",year,".asc"))
  DEFmean <- raster(paste0(outputDir,"/",className,"/MeanAnnualDEF",year,".asc"))
  DryIdx <- raster(paste0(outputDir,"/",className,"/MeanAnnualDI",year,".asc"))
  
  #summarise raster data by municipality
  meanCL <- zonal(ClimLim, munis.r, fun = 'mean')
  meanAC <- zonal(AgriCap, munis.r, fun = 'mean')
  meanDEFc <- zonal(DEFcount, munis.r, fun = 'mean')
  meanDEFm <- zonal(DEFmean, munis.r, fun = 'mean')
  meanDI <- zonal(DryIdx, munis.r, fun = 'mean')
  
  #join summarised raster data to vector data
  AC_map <- left_join(BRmunis, as_tibble(meanAC), by = c("CD_GEOCMUn" ="zone")) 
  CL_map <- left_join(BRmunis, as_tibble(meanCL), by = c("CD_GEOCMUn" ="zone"))
  DEFc_map <- left_join(BRmunis, as_tibble(meanDEFc), by = c("CD_GEOCMUn" ="zone")) 
  DEFm_map <- left_join(BRmunis, as_tibble(meanDEFm), by = c("CD_GEOCMUn" ="zone")) 
  DI_map <- left_join(BRmunis, as_tibble(meanDI), by = c("CD_GEOCMUn" ="zone")) 

  #print summarised CL and AC data using vector
  pdf(paste0(outputDir,"/",className,"/",year, "_ClimateAndCapitalMaps.pdf"))

  #CL plot
  pbrks = seq(from=0,to=max(CL_map$mean, na.rm=T),length.out=101)
  plot(CL_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(CL_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Climate Limitation ",year))

  #AC plot
  pbrks = seq(from=0,to=max(AC_map$mean, na.rm=T),length.out=101)
  plot(AC_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(AC_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Agriculture Capital ",year))

  #DEFc plot
  pbrks = seq(from=0,to=max(DEFc_map$mean, na.rm=T),length.out=101)
  plot(DEFc_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(DEFc_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Count DEF months ",year))

  #DEFm plot
  pbrks = seq(from=0,to=max(DEFm_map$mean, na.rm=T),length.out=101)
  plot(DEFm_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(DEFm_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Mean Annual DEF ",year))

  #DI plot
  pbrks = seq(from=0,to=max(DI_map$mean, na.rm=T),length.out=101)
  plot(DI_map["mean"], pal = viridis(100), breaks = pbrks, graticule = st_crs(DI_map), axes = T, lty = 0, lwd=0.2, border="lightgray", main=paste0("Mean Annual DI ",year))

  dev.off()
  #######################################################
  
  
  
  ####### SOY  ###################
  #clean soy data (including remove top and bottom 0.5% of points)
  soy_yr <- cleanData(soy)
  
  #join productivity data to summarised raster data
  joinCL <- left_join(as_tibble(meanCL), soy_yr, by = c("zone" = "CODE")) 
  joinCL <- joinCL %>% filter(!is.na(prod))
  
  joinAC <- left_join(as_tibble(meanAC), soy_yr, by = c("zone" = "CODE")) 
  joinAC <- joinAC %>% filter(!is.na(prod))
  
  joinDEFc <- left_join(as_tibble(meanDEFc), soy_yr, by = c("zone" = "CODE")) 
  joinDEFc <- joinDEFc %>% filter(!is.na(prod))
  
  joinDEFm <- left_join(as_tibble(meanDEFm), soy_yr, by = c("zone" = "CODE")) 
  joinDEFm <- joinDEFm %>% filter(!is.na(prod))
  
  joinDI <- left_join(as_tibble(meanDI), soy_yr, by = c("zone" = "CODE")) 
  joinDI <- joinDI %>% filter(!is.na(prod))
  
  #fit correlations
  analysisCL <- datSummary(joinCL)
  analysisAC <- datSummary(joinAC)
  analysisDEFc <- datSummary(joinDEFc)
  analysisDEFm <- datSummary(joinDEFm)
  analysisDI <- datSummary(joinDI)
  
  CorSoy <- data.frame()
  CorSoy <- collateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, CorSoy) 
  
  #attach table for all years
  soyData <- rbind(soyData, CorSoy)
  

  ##States
  #create tables of taus
  analysisCL <- stateSummary(joinCL)
  analysisAC <- stateSummary(joinAC)
  analysisDEFc <- stateSummary(joinDEFc)
  analysisDEFm <- stateSummary(joinDEFm)
  analysisDI <- stateSummary(joinDI)

  #create state taus table for this year
  stateCorSoy <- data.frame()
  stateCorSoy <- collateStateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, stateCorSoy) 

  #attach table for all years
  soyStates <- rbind(soyStates, stateCorSoy)
  
  
  #hack to create tables of taus
  analysisCL <- regionSummary(joinCL)
  analysisAC <- regionSummary(joinAC)
  analysisDEFc <- regionSummary(joinDEFc)
  analysisDEFm <- regionSummary(joinDEFm)
  analysisDI <- regionSummary(joinDI)

  #create state taus table for this year
  regionCorSoy <- data.frame()
  regionCorSoy <- collateRegionCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, regionCorSoy) 

  #attach table for all years
  soyRegions <- rbind(soyRegions, regionCorSoy)
  

  #create plots
  #all data
  soyCLs <- scatterPlot(joinCL, CorSoy$CLt, CorSoy$CLp, "Soy", year, "Climate Limitation")
  soyACs <- scatterPlot(joinAC, CorSoy$ACt, CorSoy$ACp, "Soy", year, "Agriculture Captial")
  soyDEFcs <- scatterPlot(joinDEFc, CorSoy$DEFct, CorSoy$DEFcp, "Soy", year, "Count DEF months")
  soyDEFms <- scatterPlot(joinDEFm, CorSoy$DEFmt, CorSoy$DEFmp, "Soy", year, "Mean Annual DEF")
  soyDIs <- scatterPlot(joinDI, CorSoy$DIt, CorSoy$DIp, "Soy", year, "Mean Annual DI")

  #region
  soyCLrs <- regionScatter(joinCL, "Soy", year, "Climate Limitation")
  soyCLrsf <- soyCLrs + facet_grid(regionNM~.)

  soyACrs <- regionScatter(joinAC, "Soy", year, "Agricultural Capital")
  soyACrsf <- soyACrs + facet_grid(regionNM~.)

  soyDEFcrs <- regionScatter(joinDEFc, "Soy", year, "Count DEF months")
  soyDEFcrsf <- soyDEFcrs + facet_grid(regionNM~.)

  soyDEFmrs <- regionScatter(joinDEFm, "Soy", year, "Mean Annual DEF")
  soyDEFmrsf <- soyDEFmrs + facet_grid(regionNM~.)

  soyDIrs <- regionScatter(joinDI, "Soy", year, "Mean Annual DI")
  soyDIrsf <- soyDIrs + facet_grid(regionNM~.)

  #state
  soyCLss <- stateScatter(joinCL, "Soy", year, "Climate Limitation")
  soyCLssf <- soyCLss + facet_grid(stateNM~.)

  soyACss <- stateScatter(joinAC, "Soy", year, "Agricultural Capital")
  soyACssf <- soyACss + facet_grid(stateNM~.)

  soyDEFcss <- stateScatter(joinDEFc, "Soy", year, "Count DEF months")
  soyDEFcssf <- soyDEFcss + facet_grid(stateNM~.)

  soyDEFmss <- stateScatter(joinDEFm, "Soy", year, "Mean Annual DEF")
  soyDEFmssf <- soyDEFmss + facet_grid(stateNM~.)

  soyDIss <- stateScatter(joinDI, "Soy", year, "Mean Annual DI")
  soyDIssf <- soyDIss + facet_grid(stateNM~.)

  soyScatters <- list(soyCLs, soyACs, soyDEFcs, soyDEFms, soyDIs,
    soyCLrs, soyCLrsf, soyACrs, soyACrsf, soyDEFcrs, soyDEFcrsf, soyDEFmrs, soyDEFmrsf,soyDIrs,soyDIrsf,
    soyCLss, soyCLssf, soyACss, soyACssf, soyDEFcss, soyDEFcssf, soyDEFmss, soyDEFmssf,soyDIss,soyDIssf)

  #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)

  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone"))
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone"))

  #output plots and maps to pdf
  soypath <- paste0(outputDir,"/",className,"/",year,"_SoyPlots.pdf")
  pdfPlots(soypath, "Soy", soyScatters, CL_map, AC_map)
  #######################################################

  
  ####### MAIZE FIRST ###################
  #maize1
  maize1_yr <- cleanData(maize1)
  
  #join productivity data to summarised raster data
  joinCL <- left_join(as_tibble(meanCL), maize1_yr, by = c("zone" = "CODE")) 
  joinCL <- joinCL %>% filter(!is.na(prod))
  
  joinAC <- left_join(as_tibble(meanAC), maize1_yr, by = c("zone" = "CODE")) 
  joinAC <- joinAC %>% filter(!is.na(prod))
  
  joinDEFc <- left_join(as_tibble(meanDEFc), maize1_yr, by = c("zone" = "CODE")) 
  joinDEFc <- joinDEFc %>% filter(!is.na(prod))
  
  joinDEFm <- left_join(as_tibble(meanDEFm), maize1_yr, by = c("zone" = "CODE")) 
  joinDEFm <- joinDEFm %>% filter(!is.na(prod))
  
  joinDI <- left_join(as_tibble(meanDI), maize1_yr, by = c("zone" = "CODE")) 
  joinDI <- joinDI %>% filter(!is.na(prod))
  
  #fit correlations
  analysisCL <- datSummary(joinCL)
  analysisAC <- datSummary(joinAC)
  analysisDEFc <- datSummary(joinDEFc)
  analysisDEFm <- datSummary(joinDEFm)
  analysisDI <- datSummary(joinDI)
  
  CorMz1 <- data.frame()
  CorMz1 <- collateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, CorMz1) 
  
  #attach table for all years
  mz1Data <- rbind(mz1Data, CorMz1)


  ##States
  #create tables of taus
  analysisCL <- stateSummary(joinCL)
  analysisAC <- stateSummary(joinAC)
  analysisDEFc <- stateSummary(joinDEFc)
  analysisDEFm <- stateSummary(joinDEFm)
  analysisDI <- stateSummary(joinDI)

  #create state taus table for this year
  stateCorMz1 <- data.frame()
  stateCorMz1 <- collateStateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, stateCorMz1) 

  #attach table for all years
  mz1States <- rbind(mz1States, stateCorMz1)
  
  
  #hack to create tables of taus
  analysisCL <- regionSummary(joinCL)
  analysisAC <- regionSummary(joinAC)
  analysisDEFc <- regionSummary(joinDEFc)
  analysisDEFm <- regionSummary(joinDEFm)
  analysisDI <- regionSummary(joinDI)

  #create state taus table for this year
  regionCorMz1 <- data.frame()
  regionCorMz1 <- collateRegionCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, regionCorMz1) 

  #attach table for all years
  mz1Regions <- rbind(mz1Regions, regionCorMz1)
  
  
  
  #create plots
  #all data
  mz1CLs <- scatterPlot(joinCL, CorMz1$CLt, CorMz1$CLp, "mz1", year, "Climate Limitation")
  mz1ACs <- scatterPlot(joinAC, CorMz1$ACt, CorMz1$ACp, "mz1", year, "Agriculture Captial")
  mz1DEFcs <- scatterPlot(joinDEFc, CorMz1$DEFct, CorMz1$DEFcp, "mz1", year, "Count DEF months")
  mz1DEFms <- scatterPlot(joinDEFm, CorMz1$DEFmt, CorMz1$DEFmp, "mz1", year, "Mean Annual DEF")
  mz1DIs <- scatterPlot(joinDI, CorMz1$DIt, CorMz1$DIp, "mz1", year, "Mean Annual DI")

  #region
  mz1CLrs <- regionScatter(joinCL, "mz1", year, "Climate Limitation")
  mz1CLrsf <- mz1CLrs + facet_grid(regionNM~.)

  mz1ACrs <- regionScatter(joinAC, "mz1", year, "Agricultural Capital")
  mz1ACrsf <- mz1ACrs + facet_grid(regionNM~.)

  mz1DEFcrs <- regionScatter(joinDEFc, "mz1", year, "Count DEF months")
  mz1DEFcrsf <- mz1DEFcrs + facet_grid(regionNM~.)

  mz1DEFmrs <- regionScatter(joinDEFm, "mz1", year, "Mean Annual DEF")
  mz1DEFmrsf <- mz1DEFmrs + facet_grid(regionNM~.)

  mz1DIrs <- regionScatter(joinDI, "mz1", year, "Mean Annual DI")
  mz1DIrsf <- mz1DIrs + facet_grid(regionNM~.)

  #state
  mz1CLss <- stateScatter(joinCL, "mz1", year, "Climate Limitation")
  mz1CLssf <- mz1CLss + facet_grid(stateNM~.)

  mz1ACss <- stateScatter(joinAC, "mz1", year, "Agricultural Capital")
  mz1ACssf <- mz1ACss + facet_grid(stateNM~.)

  mz1DEFcss <- stateScatter(joinDEFc, "mz1", year, "Count DEF months")
  mz1DEFcssf <- mz1DEFcss + facet_grid(stateNM~.)

  mz1DEFmss <- stateScatter(joinDEFm, "mz1", year, "Mean Annual DEF")
  mz1DEFmssf <- mz1DEFmss + facet_grid(stateNM~.)

  mz1DIss <- stateScatter(joinDI, "mz1", year, "Mean Annual DI")
  mz1DIssf <- mz1DIss + facet_grid(stateNM~.)

  mz1Scatters <- list(mz1CLs, mz1ACs, mz1DEFcs, mz1DEFms, mz1DIs,
    mz1CLrs, mz1CLrsf, mz1ACrs, mz1ACrsf, mz1DEFcrs, mz1DEFcrsf, mz1DEFmrs, mz1DEFmrsf,mz1DIrs,mz1DIrsf,
    mz1CLss, mz1CLssf, mz1ACss, mz1ACssf, mz1DEFcss, mz1DEFcssf, mz1DEFmss, mz1DEFmssf,mz1DIss,mz1DIssf)


  #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)

  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone"))
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone"))

  #output plots and maps to pdf
  mz1path <- paste0(outputDir,"/",className,"/",year,"_Maize1Plots.pdf")
  pdfPlots(mz1path, "Maize First",mz1Scatters, CL_map, AC_map)
  #######################################################
  
  
  
  ####### MAIZE SECOND ###################
  #maize2
  maize2_yr <- cleanData(maize2)
  
  #join productivity data to summarised raster data
  joinCL <- left_join(as_tibble(meanCL), maize2_yr, by = c("zone" = "CODE")) 
  joinCL <- joinCL %>% filter(!is.na(prod))
  
  joinAC <- left_join(as_tibble(meanAC), maize2_yr, by = c("zone" = "CODE")) 
  joinAC <- joinAC %>% filter(!is.na(prod))
  
  joinDEFc <- left_join(as_tibble(meanDEFc), maize2_yr, by = c("zone" = "CODE")) 
  joinDEFc <- joinDEFc %>% filter(!is.na(prod))
  
  joinDEFm <- left_join(as_tibble(meanDEFm), maize2_yr, by = c("zone" = "CODE")) 
  joinDEFm <- joinDEFm %>% filter(!is.na(prod))
  
  joinDI <- left_join(as_tibble(meanDI), maize2_yr, by = c("zone" = "CODE")) 
  joinDI <- joinDI %>% filter(!is.na(prod))
  
  #fit correlations
  analysisCL <- datSummary(joinCL)
  analysisAC <- datSummary(joinAC)
  analysisDEFc <- datSummary(joinDEFc)
  analysisDEFm <- datSummary(joinDEFm)
  analysisDI <- datSummary(joinDI)
  
  CorMz2 <- data.frame()
  CorMz2 <- collateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, CorMz2) 
  
  #attach table for all years
  mz2Data <- rbind(mz2Data, CorMz2)


  ##States
  #create tables of taus
  analysisCL <- stateSummary(joinCL)
  analysisAC <- stateSummary(joinAC)
  analysisDEFc <- stateSummary(joinDEFc)
  analysisDEFm <- stateSummary(joinDEFm)
  analysisDI <- stateSummary(joinDI)

  #create state taus table for this year
  stateCorMz2 <- data.frame()
  stateCorMz2 <- collateStateCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, stateCorMz2) 

  #attach table for all years
  mz2States <- rbind(mz2States, stateCorMz2)
  
  
  #hack to create tables of taus
  analysisCL <- regionSummary(joinCL)
  analysisAC <- regionSummary(joinAC)
  analysisDEFc <- regionSummary(joinDEFc)
  analysisDEFm <- regionSummary(joinDEFm)
  analysisDI <- regionSummary(joinDI)

  #create state taus table for this year
  regionCorMz2 <- data.frame()
  regionCorMz2 <- collateRegionCor(analysisCL, analysisAC, analysisDEFc, analysisDEFm, analysisDI, regionCorMz2) 

  #attach table for all years
  mz2Regions <- rbind(mz2Regions, regionCorMz2)
  
    
  
  # #create plots
  # #all data
  mz2CLs <- scatterPlot(joinCL, CorMz2$CLt, CorMz2$CLp, "mz2", year, "Climate Limitation")
  mz2ACs <- scatterPlot(joinAC, CorMz2$ACt, CorMz2$ACp, "mz2", year, "Agriculture Captial")
  mz2DEFcs <- scatterPlot(joinDEFc, CorMz2$DEFct, CorMz2$DEFcp, "mz2", year, "Count DEF months")
  mz2DEFms <- scatterPlot(joinDEFm, CorMz2$DEFmt, CorMz2$DEFmp, "mz2", year, "Mean Annual DEF")
  mz2DIs <- scatterPlot(joinDI, CorMz2$DIt, CorMz2$DIp, "mz2", year, "Mean Annual DI")

  #region
  mz2CLrs <- regionScatter(joinCL, "mz2", year, "Climate Limitation")
  mz2CLrsf <- mz2CLrs + facet_grid(regionNM~.)

  mz2ACrs <- regionScatter(joinAC, "mz2", year, "Agricultural Capital")
  mz2ACrsf <- mz2ACrs + facet_grid(regionNM~.)

  mz2DEFcrs <- regionScatter(joinDEFc, "mz2", year, "Count DEF months")
  mz2DEFcrsf <- mz2DEFcrs + facet_grid(regionNM~.)

  mz2DEFmrs <- regionScatter(joinDEFm, "mz2", year, "Mean Annual DEF")
  mz2DEFmrsf <- mz2DEFmrs + facet_grid(regionNM~.)

  mz2DIrs <- regionScatter(joinDI, "mz2", year, "Mean Annual DI")
  mz2DIrsf <- mz2DIrs + facet_grid(regionNM~.)

  #state
  mz2CLss <- stateScatter(joinCL, "mz2", year, "Climate Limitation")
  mz2CLssf <- mz2CLss + facet_grid(stateNM~.)

  mz2ACss <- stateScatter(joinAC, "mz2", year, "Agricultural Capital")
  mz2ACssf <- mz2ACss + facet_grid(stateNM~.)

  mz2DEFcss <- stateScatter(joinDEFc, "mz2", year, "Count DEF months")
  mz2DEFcssf <- mz2DEFcss + facet_grid(stateNM~.)

  mz2DEFmss <- stateScatter(joinDEFm, "mz2", year, "Mean Annual DEF")
  mz2DEFmssf <- mz2DEFmss + facet_grid(stateNM~.)

  mz2DIss <- stateScatter(joinDI, "mz2", year, "Mean Annual DI")
  mz2DIssf <- mz2DIss + facet_grid(stateNM~.)

  mz2Scatters <- list(mz2CLs, mz2ACs, mz2DEFcs, mz2DEFms, mz2DIs,
    mz2CLrs, mz2CLrsf, mz2ACrs, mz2ACrsf, mz2DEFcrs, mz2DEFcrsf, mz2DEFmrs, mz2DEFmrsf,mz2DIrs,mz2DIrsf,
    mz2CLss, mz2CLssf, mz2ACss, mz2ACssf, mz2DEFcss, mz2DEFcssf, mz2DEFmss, mz2DEFmssf,mz2DIss,mz2DIssf)

  #fit regressions
  CLfit <- lm(prod ~ mean, joinCL)
  CLfit <- broom::augment(CLfit, joinCL)
  ACfit <- lm(prod ~mean, joinAC)
  ACfit <- broom::augment(ACfit, joinAC)

  #join model data to spatial data
  AC_map <- left_join(BRmunis, ACfit, by = c("CD_GEOCMUn" ="zone"))
  CL_map <- left_join(BRmunis, CLfit, by = c("CD_GEOCMUn" ="zone"))

  #output plots and maps to pdf
  mz2path <- paste0(outputDir,"/",className,"/",year,"_Maize2Plots.pdf")
  pdfPlots(mz2path, "Maize Second",mz2Scatters, CL_map, AC_map)
  #######################################################

}

#dat <- cbind(years,soyCL,soyCLp,soyAC,soyACp,mz1CL,mz1CLp,mz1AC,mz1ACp,mz2CL,mz2CLp,mz2AC,mz2ACp)
write_csv(soyData, paste0(outputDir,"/",className,"/CorrelationsSoy.csv"))
write_csv(mz1Data, paste0(outputDir,"/",className,"/CorrelationsMz1.csv"))
write_csv(mz2Data, paste0(outputDir,"/",className,"/CorrelationsMz2.csv"))

write_csv(soyStates, paste0(outputDir,"/",className,"/StateCorrelationsSoy.csv"))
write_csv(mz1States, paste0(outputDir,"/",className,"/StateCorrelationsMz1.csv"))
write_csv(mz2States, paste0(outputDir,"/",className,"/StateCorrelationsMz2.csv"))

write_csv(soyRegions, paste0(outputDir,"/",className,"/RegionCorrelationsSoy.csv"))
write_csv(mz1Regions, paste0(outputDir,"/",className,"/RegionCorrelationsMz1.csv"))
write_csv(mz2Regions, paste0(outputDir,"/",className,"/RegionCorrelationsMz2.csv"))


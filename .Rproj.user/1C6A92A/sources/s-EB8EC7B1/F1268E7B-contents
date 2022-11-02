rm(list = ls())

library(ncdf4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lubridate)
library(Yoko.regrowth)
library(zoo)
library(tidyr)
library(PEcAn.ED2)

start_year <- 1850 ; end_year <- 2021

start_date <- paste0(start_year,"-01-01")
end_date <- paste0(end_year,"-12-31")
slat = 0.3
slon = 25.3
overwrite = TRUE

destination <- "/data/gent/vo/000/gvo00074/felicien/R/data/"

site.XTS <- readRDS(file = file.path(destination,"ERA5_reanalysis.RDS"))

for (year in seq(start_year,end_year)){

  print(year)

  out <- Yoko.regrowth::met2CF.ERA5(
    lat = slat,
    long = slon,
    start_date =  paste0(year,"-01-01"),
    end_date = paste0(year,"-12-31"),
    sitename = "Yoko",
    outfolder = destination,
    out.xts = list(site.XTS %>% subset(year(index(.)) == year)),
    overwrite = overwrite,
    verbose = FALSE,
    type = "reanalysis")
}

for (year in seq(start_year,end_year)){

  print(year)
  Yoko.regrowth::mod.met2model.ED2(
    in.path =  file.path(destination,"ERA5_Yoko_1"),
    in.prefix = "ERA5.1",
    outfolder = file.path(destination,"ERA5_Yoko_1","ED2"),
    start_date =  paste0(year,"-01-01"),
    end_date = paste0(year,"-12-31"),
    lat = slat,
    lon = slon,
    lst = 0,
    overwrite = overwrite)
}

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/process.ERA5.Yoko.HPC.R hpc:/data/gent/vo/000/gvo00074/felicien/R/

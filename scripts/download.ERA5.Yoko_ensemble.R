rm(list = ls())

library(reticulate)
library(future)
library(purrr)
library(furrr)

library(ncdf4)

prefix <- "ERA5_Yoko_h_"

setwd("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko")
# setwd("/data/gent/vo/000/gvo00074/ED_common_data/met/Yoko/ERA5/")

plan(multiprocess)

files.already.dow <- tools::file_path_sans_ext(list.files(getwd(),
                                                          pattern = prefix))
yrs <- as.numeric(sub(".*\\_","",files.already.dow))
yrs2download <- 1960:2021
years <- yrs2download[!(yrs2download %in% yrs)]

c(years) %>%
  future_map(function(year) {

    # you need to have an account for downloaing the files
    # Read the documantion for how to setup your account and settings before trying this
    # https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-3-DownloadERA5datathroughtheCDSAPI
    cdsapi <-import("cdsapi")
    c <- cdsapi$Client()

    c$retrieve(
      'reanalysis-era5-single-levels',
      list(
        'product_type' = 'reanalysis',
        'format' = 'netcdf',
        'day' = list('01','02','03',
                     '04','05','06',
                     '07','08','09',
                     '10','11','12',
                     '13','14','15',
                     '16','17','18',
                     '19','20','21',
                     '22','23','24',
                     '25','26','27',
                     '28','29','30',
                     '31'),
        # 'time' = list('00:00','03:00','06:00',
        #               '09:00','12:00','15:00',
        #               '18:00','21:00'),
        'time' = list('00:00','01:00','02:00',
                      '03:00','04:00','05:00',
                      '06:00','07:00','08:00',
                      '09:00','10:00','11:00',
                      '12:00','13:00','14:00',
                      '15:00','16:00','17:00',
                      '18:00','19:00','20:00',
                      '21:00','22:00','23:00'),
        'month' = list('01','02','03',
                       '04','05','06',
                       '07','08','09',
                       '10','11','12'),
        'year' = as.character(year),
        'area' = "0.35/25.25/0.25/25.35",
        'grid' = "0.1/0.1",
        'variable' = list( "2m_temperature",
                           "surface_pressure",
                           "2m_dewpoint_temperature",
                           "total_precipitation",
                           "10m_u_component_of_wind",
                           "10m_v_component_of_wind",
                           "surface_solar_radiation_downwards",
                           "surface_thermal_radiation_downwards")
      ),
      paste0(prefix,"_",year,'.nc')
    )
  })


ncfile <- paste0("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/",
                 prefix,"_1959.nc")
nc <- nc_open(ncfile)
lon <- ncvar_get(nc,"longitude")
lat <- ncvar_get(nc,"latitude")
nc_close(ncfile)

# reanalysis, ensemble_members
# "10/-0/-10/40"

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/download.ERA5.Yoko_ensemble.R hpc:/data/gent/vo/000/gvo00074/felicien/R

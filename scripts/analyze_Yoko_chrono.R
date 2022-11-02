rm(list = ls())

library(dplyr)
library(tidyr)
library(PEcAn.ED2)
library(stringr)
library(raster)
library(rhdf5)
library(pracma)

source("/data/gent/vo/000/gvo00074/felicien/R/h5read_opt.r")

DIR <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/histo/"

cdf <- data.frame()

years <- seq(1700,2300,5)
for (iyear in seq(1,length(years))){

  print(iyear/length(years))

  for (imonth in seq(1,12)){

    h5file <- file.path(DIR,
                        paste("history","S",as.character(years[iyear]),sprintf("%02d",imonth),"01-000000-g01.h5",sep = "-"))

    if (file.exists(h5file)){
      mymont    = lapply(h5read_opt(h5file),FUN=aperm)
      names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".")

      cdf <- bind_rows(list(cdf,
                            data.frame(yr = years[iyear],
                                       month = imonth,
                                       AGB = sum(mymont$AGB.PY),
                                       LAI = sum(mymont$LAI.PY))))

    }
  }
}


saveRDS(cdf,"df_Yoko_chronosequence.RDS")

# scp /home/femeunier/Documents/projects/YGB/scripts/analyze_Yoko_chrono.R hpc:/data/gent/vo/000/gvo00074/felicien/R

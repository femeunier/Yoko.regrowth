# Step 1
# ml purge ; source load.sh ; R

rm(list = ls())

library(reticulate)
library(ncdf4)

use_python("/kyukon/home/apps/RHEL8/skylake-ib/software/Python/3.9.6-GCCcore-11.2.0/bin/python",
           required = T)
py_config()

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/mod_netcdf_var.py hpc:/data/gent/vo/000/gvo00074/felicien/R/scripts/
source_python("./scripts/mod_netcdf_var.py")

maindir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble/ERA5_Yoko_ensemble_"
Nensemble <- 10

# maindir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.reanalysis/ERA5_Yoko_"
# Nensemble <- 1

years2add <- 1700:1709

for (i in seq(1,Nensemble)){
  dirpath <- paste0(maindir,i)
  years2replace <- seq(1850,1859)
  for (iyear in seq(1,length(years2add))){

    if (lubridate::leap_year(years2add[iyear])){
      yearselect <- sample(rep(years2replace[lubridate::leap_year(years2replace)],2),1)
    } else {
      yearselect <- sample(rep(years2replace[!lubridate::leap_year(years2replace)],2),1)
    }
    print(yearselect)
    years2replace <- years2replace[years2replace != yearselect]

    fname <- paste0("ERA5.",i,".",yearselect,".nc")
    fname.mod <- paste0("ERA5.",i,".",years2add[iyear],".nc")

    met_driver <- file.path(dirpath,fname)
    met_driver.mod <- file.path(dirpath,fname.mod)

    if (file.exists(met_driver.mod)) system2("rm",paste(met_driver.mod))
    system2("cp",paste(met_driver,met_driver.mod))

    y = mod_netcdf_var(met_driver.mod,
                       "mole_fraction_of_carbon_dioxide_in_air",
                       280./1e6)
  }
}

# Step 2

# ml purge ; ml R/4.1.0-foss-2021a ; ml NCO/5.0.1-foss-2021a ; R

rm(list = ls())

maindir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble/ERA5_Yoko_ensemble_"
Nensemble <- 10

years2add <- 1700:1709

for (i in seq(1,Nensemble)){
  dirpath <- paste0(maindir,i)

  for (iyear in seq(1,length(years2add))){

    fname.mod <- paste0("ERA5.",i,".",years2add[iyear],".nc")
    met_driver.mod <- file.path(dirpath,fname.mod)

    system2("ncatted",paste0("-h -a units,time,o,c,'hours since ",years2add[iyear],
                             "-01-01T00:00' ",
                             maindir,i,"/ERA5.",i,".",years2add[iyear],".nc"))


  }
}


# Step 3
# ml purge ; source load.sh ; R

rm(list = ls())

maindir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble/ERA5_Yoko_ensemble_"
Nensemble <- 10

slat = 0.3
slon = 25.3
overwrite = TRUE
years2add <- 1700:1709

for (i in seq(1,Nensemble)){

  for (year in seq(min(years2add),max(years2add))){

    print(year)
    Yoko.regrowth::mod.met2model.ED2(
      in.path =  paste0(maindir,i),
      in.prefix = paste0("ERA5.",i),
      outfolder = paste0(maindir,i,"/ED2"),
      start_date =  paste0(year,"-01-01"),
      end_date = paste0(year,"-12-31"),
      lat = slat,
      lon = slon,
      lst = 0,
      overwrite = overwrite)
  }
}

# nc <- nc_open(met_driver.mod)

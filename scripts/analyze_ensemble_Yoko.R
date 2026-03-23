rm(list = ls())

library(dplyr)
library(LidarED)
library(purrr)
library(stringr)
library(reshape2)
library(rhdf5)

rundir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run/ERA5.ensemble"   # Directory for the run folders
outdir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble"                      # Directory for the outfolders folders

Nensemble <- 10
df <- data.frame()

source("/data/gent/vo/000/gvo00074/felicien/R/h5read_opt.r")

for (iens in seq(1,Nensemble)){

  run_name <- paste0("Ensemble_",iens)

  run_ref <- file.path(rundir,run_name)
  out_ref <- file.path(outdir,run_name)

  details.file <- file.info(list.files(path = file.path(out_ref,"histo"), full.names = TRUE,pattern = ".h5"))

  if (nrow(details.file)>0){
    files.OP.ordered <- details.file[with(details.file, order(as.POSIXct(mtime),decreasing = TRUE)), ]
    h5file <- file.path(rownames(files.OP.ordered)[1])
    h5file.name <- basename(h5file)

    final.year <- as.numeric(stringr::str_split(h5file.name,"-")[[1]][3])

    for (year in seq(1550,final.year,5)){

      h5file <- file.path(out_ref,"histo",paste0("history-S-",year,"-01-01-000000-g01.h5"))

      mymont    = lapply(h5read_opt(h5file),FUN=aperm)
      names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".")

      AGB <- sum(mymont$AGB.PY)
      AGB.tree <- sum(mymont$AGB.PY[1,,c(2,3,4)])
      LAI <- sum(mymont$LAI.PY)
      LAI.tree <- sum(mymont$LAI.PY[1,,c(2,3,4)])

      df <- bind_rows(list(df,
                           data.frame(ens = iens,
                                      year,
                                      AGB,
                                      AGB.tree,
                                      LAI,
                                      LAI.tree)))
    }
  }
}

out_ref <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/"

details.file <- file.info(list.files(path = file.path(out_ref,"histo"), full.names = TRUE,pattern = ".h5"))

if (nrow(details.file)>0){
  files.OP.ordered <- details.file[with(details.file, order(as.POSIXct(mtime),decreasing = TRUE)), ]
  h5file <- file.path(rownames(files.OP.ordered)[1])
  h5file.name <- basename(h5file)

  final.year <- as.numeric(stringr::str_split(h5file.name,"-")[[1]][3])

  for (year in seq(1550,final.year,5)){

    h5file <- file.path(out_ref,"histo",paste0("Yoko_default-S-",year,"-01-01-000000-g01.h5"))

    mymont    = lapply(h5read_opt(h5file),FUN=aperm)
    names(mymont) <- gsub(x = names(mymont), pattern = "\\_", replacement = ".")

    AGB <- sum(mymont$AGB.PY)
    AGB.tree <- sum(mymont$AGB.PY[1,,c(2,3,4)])

    LAI <- sum(mymont$LAI.PY)
    LAI.tree <- sum(mymont$LAI.PY[1,,c(2,3,4)])

    df <- bind_rows(list(df,
                         data.frame(ens = 0,
                                    year,
                                    AGB,
                                    AGB.tree,
                                    LAI,
                                    LAI.tree)))
  }
}

saveRDS(object = df,file = file.path('.',"df_OP_ensemble_Yoko.RDS"))

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/analyze_ensemble_Yoko.R hpc:/data/gent/vo/000/gvo00074/felicien/R


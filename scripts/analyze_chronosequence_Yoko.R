rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(ED2scenarios)
library(PEcAn.ED2)
library(purrr)
library(cowplot)
library(pracma)
library(lubridate)
library(rhdf5)

source("/data/gent/vo/000/gvo00074/felicien/R/h5read_opt.r")

ref_dir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run/"
rundir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/run/"
outdir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/out/"

t.since.disturbance <- c(5,12,20,30,40,50,60,Inf)
all.OP <- data.frame()


for (i in seq(1,length(t.since.disturbance))){

  print(i/length(t.since.disturbance))

  if (t.since.disturbance[i] < Inf){
    run_name <- paste0("plot.",t.since.disturbance[i],".yr.old")

    run_ref <- file.path(rundir,run_name)
    out_ref <- file.path(outdir,run_name)

    analy.dir <- file.path(out_ref,"analy")
    histo.dir <- file.path(out_ref,"histo")

    years2check <- seq(2020-t.since.disturbance[i],2020,1)

    basename <- "history"

  } else {

    out_ref <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/"

    basename <- "Yoko_default"

    analy.dir <- file.path(out_ref,"analy")
    histo.dir <- file.path(out_ref,"histo")

    years2check <- 1550:2020
  }

  for (cyear in years2check){

    h5file <- file.path(histo.dir,paste0(basename,"-S-",cyear,"-01-01-000000-g01.h5"))

    if (!file.exists(h5file)) next()
    mymont <- lapply(h5read_opt(h5file),FUN=aperm)

    cAGB.census <- sum(mymont$AGB_PY[,2:11,])
    cAGB <- sum(mymont$AGB_PY)
    all.OP <- bind_rows(all.OP,
                        data.frame(year = cyear,
                                   timing = t.since.disturbance[i],
                                   AGB = cAGB,
                                   AGB.census = cAGB.census))
  }
}

saveRDS(all.OP,"./outputs/AGB_chronoseq_Yoko.RDS")
# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/analyze_chronosequence_Yoko.R hpc:/data/gent/vo/000/gvo00074/felicien/R


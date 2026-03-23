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

# http://openknowledge.nau.edu/id/eprint/2916/6/Kauffman_J_etal_1993_Biomass_and_Nutrient_Dynamics_Associated_with_Slash_Fires_in_Neotropical_Dry_Forests(1).pdf
# https://link.springer.com/content/pdf/10.1007/BF00341336.pdf


# If we want to modify the initial file
# source("/home/femeunier/Documents/ED2/R-utils/h5read_opt.r")
#
# use_python("/usr/bin/python3.8", required = T)
# py_config()
#
# source_python("./scripts/mod_netcdf_var.py")


ref_dir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run/"
rundir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/run/"
outdir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/out/"

ed2in <- read_ed2in(file.path(ref_dir,"ED2IN_Yoko_default_history"))
ed2in$IMONTHA = 1
ed2in$IDATEA = 1
ed2in$IYEARA = 1700

ed2in$IMONTHZ = 2
ed2in$IDATEZ = 1
ed2in$IYEARZ = 2020

ed2in$RUNTYPE <- "HISTORY"
ed2in$IMOUTPUT <- 3
ed2in$IED_INIT_MODE <- 6

# ed2in$SFILIN <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/histo/Yoko_historical"
ed2in$SFILIN <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/histo/Yoko_disturbed"

sfilin2copy <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/histo/Yoko_default-S-1550-01-01-000000-g01.h5"

list_dir <- list()

Nsimuperjob = 1
isimu = 0

t.since.disturbance <- c(5,12,20,60,30,40,50)

for (i in seq(1,length(t.since.disturbance))){

  run_name <- paste0("plot.",t.since.disturbance[i],".yr.old")

  run_ref <- file.path(rundir,run_name)
  out_ref <- file.path(outdir,run_name)

  isimu = isimu + 1

  if(!dir.exists(run_ref)) dir.create(run_ref)
  if(!dir.exists(out_ref)) dir.create(out_ref)
  if(!dir.exists(file.path(out_ref,"analy"))) dir.create(file.path(out_ref,"analy"))
  if(!dir.exists(file.path(out_ref,"histo"))) dir.create(file.path(out_ref,"histo"))

  # ED2IN
  ed2in_scenar <- ed2in

  ed2in_scenar$IEDCNFGF <- file.path(run_ref,"config.xml")
  ed2in_scenar$FFILOUT = file.path(out_ref,"analy","analysis")
  ed2in_scenar$SFILOUT = file.path(out_ref,"histo","history")

  ed2in_scenar$IYEARH <- 2020 - t.since.disturbance[i]
  ed2in_scenar$IMONTHH <- 1


  # We need to modify the history file here, for now we just copy an empty one
  system2("cp",c(sfilin2copy,
                 paste0(ed2in$SFILIN,"-S-", ed2in_scenar$IYEARH,"-01-01-000000-g01.h5")))

  write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN"))

  if (isimu == 1){
    isfirstjob = TRUE
    dir_joblauncher = run_ref
    list_dir[[run_name]] = run_ref
  } else{
    isfirstjob = FALSE
  }

  # job.sh

  write_joblauncher(file =  file.path(dir_joblauncher,"job.sh"),
                    nodes = 1,ppn = 18,mem = 16,walltime = 12,
                    prerun = "ml purge ; ml intel-compilers/2021.4.0 HDF5/1.12.1-iimpi-2021b UDUNITS/2.2.28-GCCcore-11.2.0; ulimit -s unlimited",
                    CD = run_ref,
                    ed_exec = "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2.2/ED2/ED/build/ed_2.2-opt-master-fa80dab6",
                    Rplot_function = NULL,
                    ED2IN = "ED2IN",
                    firstjob = isfirstjob,
                    clean = TRUE,
                    in.line = 'ml purge; ml R/4.1.2-foss-2021b',
                    reload = TRUE)

  if (isimu == Nsimuperjob){
    isimu = 0
  }
}

dumb <- write_bash_submission(file = file.path(rundir,"all_jobs_chronosequence.sh"),
                              list_files = list_dir,
                              job_name = "job.sh")

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/generate_chronosequence_Yoko.R hpc:/data/gent/vo/000/gvo00074/felicien/R


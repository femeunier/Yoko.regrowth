rm(list = ls())

# Libraries
library(dplyr)
library(tidyr)
library(purrr)
library(ED2scenarios)
library(PEcAn.ED2)
library(BayesianTools)

# Directories
ref_dir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run"
ed2in <- read_ed2in(file.path(ref_dir,"ED2IN_Yoko_default"))                                # reference ED2IN file

# Global config
ed2in$RUNTYPE  <- "INITIAL"
ed2in$IED_INIT_MODE <- 0
ed2in$ITOUTPUT <- 0
ed2in$IYEARZ <- 1850

rundir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run/ERA5.ensemble"   # Directory for the run folders
outdir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble"                       # Directory for the outfolders folders

dir.create(rundir,showWarnings = FALSE)
dir.create(outdir,showWarnings = FALSE)

# Multi jobs --> to group the simulations in a single job file
Nsimuperjob = 1
isimu = 0

# Config file
PREFIX_XML <- "<?xml version=\"1.0\"?>\n<!DOCTYPE config SYSTEM \"ed.dtd\">\n"
defaults <- list_dir <- list()

# Default settings
settings <- list(model = list(revision = "git",
                              config.header = NULL),
                 pfts = list(pft = list(num = 2,
                                        ed2_pft_number = 2,
                                        name = "Early"),
                             pft = list(num = 3,
                                        ed2_pft_number = 3,
                                        name = "Mid"),
                             pft = list(num = 4,
                                        ed2_pft_number = 4,
                                        name = "Late")))

# Default config
config <- list()
config[["Early"]] <- unlist(list(num = 2))
config[["Mid"]] <- unlist(list(num = 3))
config[["Late"]] <- unlist(list(num = 4))

#################################################################
# Main loops for the ensemble runs

Nensemble <- 10

for (iens in seq(1,Nensemble)){

  # Directories
  run_name <- paste0("Ensemble_",iens)
  isimu = isimu + 1

  run_ref <- file.path(rundir,run_name)
  out_ref <- file.path(outdir,run_name)

  if(!dir.exists(run_ref)) dir.create(run_ref)
  if(!dir.exists(out_ref)) dir.create(out_ref)
  if(!dir.exists(file.path(out_ref,"analy"))) dir.create(file.path(out_ref,"analy"))
  if(!dir.exists(file.path(out_ref,"histo"))) dir.create(file.path(out_ref,"histo"))

  # ED2IN
  ed2in_scenar <- ed2in
  ed2in_scenar$IEDCNFGF <- file.path(run_ref,"config.xml")
  ed2in_scenar$FFILOUT = file.path(out_ref,"analy","analysis")
  ed2in_scenar$SFILOUT = file.path(out_ref,"histo","history")
  ed2in_scenar$ED_MET_DRIVER_DB <- paste0("/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/ERA5.ensemble/ERA5_Yoko_ensemble_",iens,
                                          "/ED2/ED_MET_DRIVER_HEADER")

  write_ed2in(ed2in_scenar,filename = file.path(run_ref,"ED2IN"))

  # Config
  config_simu <- config

  xml <- write.config.xml.ED2(defaults = defaults,
                              settings = settings,
                              trait.values = config_simu)

  XML::saveXML(xml, file = file.path(run_ref,"config.xml"), indent = TRUE,
               prefix = PREFIX_XML)


  # job.sh

  if (isimu == 1){
    isfirstjob = TRUE
    dir_joblauncher = run_ref
    list_dir[[run_name]] = run_ref
  } else{
    isfirstjob = FALSE
  }

  write_joblauncher(file =  file.path(dir_joblauncher,"job.sh"),
                    nodes = 1,ppn = 16,mem = 16,walltime = 72,
                    prerun = "ml purge ; ml intel-compilers/2021.4.0 HDF5/1.12.1-iimpi-2021b UDUNITS/2.2.28-GCCcore-11.2.0; ulimit -s unlimited",
                    CD = run_ref,
                    ed_exec = "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2.2/ED2/ED/build/ed_2.2-opt-master-fa80dab",
                    Rplot_function = '/data/gent/vo/000/gvo00074/felicien/R/read_and_plot_ED2.2_all_tspft_yearly.r',
                    ED2IN = "ED2IN",
                    firstjob = isfirstjob,
                    clean = TRUE,
                    in.line = 'ml purge; ml R/4.1.2-foss-2021b',
                    reload = TRUE)


  if (isimu >= Nsimuperjob){
    isimu = 0
  }
}

dumb <- write_bash_submission(file = file.path(rundir,"all_jobs_ensemble_spinup.sh"),
                              list_files = list_dir,
                              job_name = "job.sh")


# To transfer the files
# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/generate_ensemble_Yoko.R hpc:/data/gent/vo/000/gvo00074/felicien/R

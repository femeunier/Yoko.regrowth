rm(list = ls())

# Instructions go here

library(R.utils)

# 3 days license!!
# https://help.ceda.ac.uk/article/4442-ceda-opendap-scripted-interactions

perm <- "/user/gent/425/vsc42558/ceda_pydap_cert_code/online_ca_client/contrail/security/onlineca/client/sh/creds.pem"
dir <- "https://dap.ceda.ac.uk/badc/cru/data/cru_jra/cru_jra_2.5/data/"

vars <- c("ugrd","vgrd")
years <- seq(1901,2023)

dir.OP <- "/data/gent/vo/000/gvo00074/felicien/TrENDY/inputs"

for (ivar in seq(1,length(vars))){
  for (iyear in seq(1,length(years))){

    var = vars[ivar]
    year = years[iyear]

    print(paste(var,"-",year))

    filename <- paste0("crujra.v2.5.5d.",var,".",year,".365d.noc.nc")
    filename.compressed <- paste(filename,"gz",sep = ".")
    #filename.cropped <- paste0("crujra.v2.5.5d.",var,".",year,".365d.noc.cropped.nc")

    source_file <- file.path(dir,var,filename.compressed)
    dest.file <- file.path(dir.OP,filename.compressed)
    dest.uncompressed <- file.path(dir.OP,var,filename)

    if (file.exists(dest.file)){ next()}

    system2("curl",paste("--cert",
                         perm,"-L",
                         "-c /dev/null",
                         source_file,
                         "--output",
                         dest.file))

    #gunzip(dest.file)
  }
}

# dest.unziped.cropped <- file.path("./outputs",filename.cropped)
# system2("cdo",paste0("sellonlatbox",",",-10,",",45,",",-15,",",10," ",dest.uncompressed," ",dest.unziped.cropped))
# system2("rm",dest.uncompressed)

# scp /home/femeunier/Documents/projects/Yoko.regrowth/scripts/download.CRUJRA.R hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/


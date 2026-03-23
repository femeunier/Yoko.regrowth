rm(list = ls())

library(rGEDI)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(leaflet)
library(leafsync)

b_box <- c(25,25.5,0,0.5)
e <- raster::extent(b_box)
ul_lat <- ymax <- e@ymax
lr_lat <- ymin <- e@ymin
ul_lon <- xmin <- e@xmin
lr_lon <- xmax <- e@xmax
daterange=c("2010-01-01","2030-06-01")
outdir = tempdir()

gLevel2A<-rGEDI::gedifinder(product="GEDI02_A",
                            ul_lat, ul_lon,
                            lr_lat, lr_lon,
                            version="002",
                            daterange=daterange)


outdir="./outputs/GEDI.data/"

netrc = file.path(outdir, ".netrc")
netrc_conn <- file(netrc)

writeLines(c("machine urs.earthdata.nasa.gov",
             sprintf("login %s", "femeunier"),
             sprintf("password %s", "Jleconnaispas0")),
           netrc_conn)


dir.create(outdir,showWarnings = FALSE)
system2("rm",
        c(paste0(outdir,"*.curltmp")))

# Downloading GEDI data
rGEDI::gediDownload(filepath=gLevel2A,outdir=outdir)

system2("rm",file.path(outdir,"level2a_clip_bb.h5"))
gedilevel2a<-readLevel2A(level2Apath = file.path(outdir,basename(gLevel2A[1])))
level2a_clip_bb <- clipLevel2A(gedilevel2a, xmin, xmax, ymin, ymax,
                               output=file.path(outdir,"level2a_clip_bb.h5"))

gedilevel2a<-readLevel2A(level2Apath = file.path(outdir,"level2a_clip_bb.h5"))
level2AM<-getLevel2AM(level2a_clip_bb)
head(level2AM[,c("beam","shot_number","elev_highestreturn","elev_lowestmode","rh100")])

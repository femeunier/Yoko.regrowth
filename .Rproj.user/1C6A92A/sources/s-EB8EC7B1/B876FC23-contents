rm(list = ls())

library(ggplot2)
library(dplyr)
library(rhdf5)

source("/home/femeunier/Documents/ED2/R-utils/h5read_opt.r")

# system2("scp",
#         paste("hpc:/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/analy/Yoko_default-Q-1990-01-00-000000-g01.h5",
#               "./outputs/"))
# system2("scp",
#         paste("hpc:/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/out/plot.5.yr.old/histo/history-S-2015-01-01-000000-g01.h5",
#               "./outputs/"))

csspssfile_name <- "Yoko.init.lat0.500lon25.500"

file <- "Yoko_default-Q-1990-01-00-000000-g01.h5"
h5file <- paste0("./outputs/",file)

mymont <- lapply(h5read_opt(h5file),FUN=aperm)

# mymont2 <- lapply(h5read_opt("./outputs/Yoko_default-S-1550-01-01-000000-g01.h5"),
#                   FUN=aperm)
# plot(mymont$DBH,mymont$NPLANT,log = "xy")
# lines(mymont2$DBH,mymont2$NPLANT, type = "p", col = "red")
#
# hist(mymont$PFT)
# hist(mymont2$PFT)
#
# mymont$AREA
# mymont2$AREA
#
# mymont$AGE
# mymont2$AGE
#
# mymont$NPATCHES_GLOBAL
# mymont2$NPATCHES_GLOBAL
#
# mymont$NCOHORTS_GLOBAL
# mymont2$NCOHORTS_GLOBAL
#
# plot(mymont$DBH,mymont$MMEAN_LAI_CO,log = "xy")
# lines(mymont2$DBH,mymont2$LAI_CO, type = "p", col = "red")
#
# plot(mymont2$DBH,mymont2$NPLANT,log = "xy")
# lines(mymont$DBH,mymont$NPLANT, type = "p", col = "red")
#
# plot(mymont$DBH,mymont$AGB_CO,log = "xy")
# lines(mymont2$DBH,mymont2$AGB_CO, type = "p", col = "red")


###################################################################################################################################
# Patch file

Npatches <- mymont$NPATCHES_GLOBAL
year <- 1990
paco <- rep(1:Npatches,mymont$PACO_N)
pa.area <- (mymont$AREA)/sum((mymont$AREA))

patches <- data.frame(time = year,
                      patch = unique(paco),
                      trk = 2,
                      age = mymont$AGE,
                      area = pa.area,
                      water = 0,
                      fsc = mymont$MMEAN_FAST_SOIL_C_PA + mymont$MMEAN_FAST_GRND_C_PA,
                      stsc = mymont$MMEAN_STRUCT_SOIL_C_PA + mymont$MMEAN_STRUCT_GRND_C_PA,
                      stsl =  mymont$MMEAN_STRUCT_SOIL_L_PA + mymont$MMEAN_STRUCT_GRND_L_PA,
                      ssc = mymont$MMEAN_SLOW_SOIL_C_PA,
                      lai = data.frame(pa = paco,
                                       lai = mymont$MMEAN_LAI_CO) %>% group_by(pa) %>% summarise(lai = sum(lai)) %>% pull(lai),
                      msn = mymont$MMEAN_MINERAL_SOIL_N_PA,
                      fsn = mymont$MMEAN_FAST_GRND_N_PA + mymont$MMEAN_FAST_SOIL_N_PA,
                      nep = 0,
                      gpp = 0,
                      rh = 0)

# weighted.mean(mymont$MMEAN_FAST_SOIL_C_PA + mymont$MMEAN_FAST_GRND_C_PA,
#               pa.area)
# weighted.mean(mymont$MMEAN_STRUCT_SOIL_C_PA + mymont$MMEAN_STRUCT_GRND_C_PA,
#               pa.area)
# weighted.mean(mymont$MMEAN_SLOW_SOIL_C_PA,
#               pa.area)

dir.create("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
           showWarnings = FALSE)


write.table(x = patches,
            file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                             paste0(csspssfile_name,".pss")),
            row.names = FALSE,
            col.names = TRUE)

write.table(x = patches,
            file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                             paste0("Yoko.init.empty.lat0.500lon25.500",".pss")),
            row.names = FALSE,
            col.names = TRUE)


###################################################################################################################################
# Cohort file
Ncohorts <- length(mymont$AGB_CO)

cohorts <- data.frame(time = year,
                      patch = paco,
                      cohort = 1:Ncohorts,
                      dbh = mymont$DBH,
                      hite = mymont$HITE,
                      pft = mymont$PFT,
                      n = mymont$NPLANT, #*rep(pa.area,mymont$PACO_N),
                      bdead = mymont$BDEADA + mymont$BDEADB,
                      balive = mymont$MMEAN_BALIVE_CO,
                      lai = mymont$MMEAN_LAI_CO) %>% filter(dbh >= 0.2)

write.table(x = cohorts,
            file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                             paste0(csspssfile_name,".css")),
            row.names = FALSE,col.names = TRUE)

empty.cohorts <- data.frame(time = year,
                            patch = sort(rep(unique(paco),3)),
                            cohort = 1:(3*Npatches),
                            dbh = rep(c(0.17,0.16,0.15),Npatches),
                            hite = 0,
                            pft = rep(c(2,3,4),Npatches),
                            n = rep(0.1,3*Npatches),
                            bdead = 0,
                            balive = 0,
                            lai = 0)

write.table(x = empty.cohorts,
            file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                             paste0("Yoko.init.empty.lat0.500lon25.500",".css")),
            row.names = FALSE,col.names = TRUE)


system2("scp",paste(paste0(file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                                      paste0(csspssfile_name)),"*"),
                     "hpc:/data/gent/vo/000/gvo00074/ED_common_data/inits/Yoko/"))

system2("scp",paste(paste0(file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/IC/",
                                     paste0("Yoko.init.empty.lat0.500lon25.500")),"*"),
                    "hpc:/data/gent/vo/000/gvo00074/ED_common_data/inits/Yoko/"))


##################################################################################################

# system2("scp",
#         paste("hpc:/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/chronosequence/out/plot.5.yr.old/analy/analysis.RData",
#               "./outputs/"))
#
# load("./outputs/analysis.RData")
# restart <- datum
#
# load("./outputs/Yoko_default.RData")
# original <- datum
#
# matplot(original$szpft$gpp[,12,c(2,3,4,18)], type = "l", ylim = c(0,5))
# matlines(restart$szpft$gpp[,12,c(2,3,4,18)], lty = 2)
#
# matplot(original$szpft$lai[,12,c(2,3,4,18)], type = "l", ylim = c(0,7))
# matlines(restart$szpft$lai[,12,c(2,3,4,18)], lty = 2)
#
# matplot(original$szpft$agb[,12,c(2,3,4,18)], type = "l", ylim = c(0,20))
# matlines(restart$szpft$agb[,12,c(2,3,4,18)], lty = 2)
#
# matplot(original$szpft$nplant[,12,c(2,3,4,18)], type = "l", ylim = c(0,30))
# matlines(restart$szpft$nplant[,12,c(2,3,4,18)], lty = 2)
#
# plot(original$emean$fast.grnd.c + original$emean$fast.soil.c, type = "l", ylim = c(0,1))
# lines(restart$emean$fast.grnd.c + restart$emean$fast.soil.c, type = "l", col = "red")
#
# plot(original$emean$struct.grnd.c + original$emean$struct.soil.c, type = "l")
# lines(restart$emean$struct.grnd.c + restart$emean$struct.soil.c, type = "l", col = "red")
#
# plot(original$emean$slow.soil.c, type = "l")
# lines(restart$emean$slow.soil.c, type = "l", col = "red")


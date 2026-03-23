rm(list = ls())

library(ncdf4)
library(ggplot2)
library(dplyr)

nc <- nc_open(filename = file <- "/home/femeunier/Documents/projects/Yoko.regrowth/data/TRENDY/ORCHIDEE-CNP_S3_cVeg.nc")

all.lats <- ncvar_get(nc,"lat")
all.lons <- ncvar_get(nc,"lon")
all.times <- ncvar_get(nc,"time")
cVeg <- ncvar_get(nc,"cVeg")

nc_close(nc)

nc <- nc_open(filename = file <- "/home/femeunier/Documents/projects/Yoko.regrowth/data/TRENDY/ORCHIDEE-CNP_S2_cVeg.nc")

all.timesS2 <- ncvar_get(nc,"time")
cVegS2 <- ncvar_get(nc,"cVeg")

nc_close(nc)

# example Yoko
clat = 0.25 ; clon = 25.25
clon = 49.625 ; clat =  -14.875
clon = 22.875  ; clat = 4.875

dist.df <- expand.grid(all.lats,all.lons) %>% rename(lat = Var1,
                                                     lon = Var2) %>%
  mutate(dist = sqrt((lat - clat)**2 + (lon - clon)**2)) %>% arrange(dist) %>% slice_head(n = 1)

df <- data.frame(all.times,
                 cVeg = as.vector(cVeg[which.min(abs(as.vector(all.lons) - (dist.df %>% pull(lon)))),
                                       which.min(abs(as.vector(all.lats) - (dist.df %>% pull(lat)))),
                                       ]),
                 cVegS2 = as.vector(cVegS2[which.min(abs(as.vector(all.lons) - (dist.df %>% pull(lon)))),
                                       which.min(abs(as.vector(all.lats) - (dist.df %>% pull(lat)))),
                                       ]))

plot(df$all.times,df$cVegS2,ylim = c(0,1.5*max(df$cVegS2)))
lines(df$all.times,df$cVeg,col = "red")

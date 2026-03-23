rm(list = ls())

library(GEDI4R)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(leaflet)
library(leafsync)

b_box <- c(25,25.5,0,0.5)
e <- raster::extent(b_box)
ul_lat <- e@ymax
lr_lat <- e@ymin
ul_lon <- e@xmin
lr_lon <- e@xmax
daterange=c("2010-01-01","2030-06-01")
outdir = tempdir()

gLevel4A<-rGEDI::gedifinder(product="GEDI04_A",
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
rGEDI::gediDownload(filepath=gLevel4A,outdir=outdir)

gediL4 <- l4_getmulti(file.path(outdir,basename(gLevel4A)),
                           merge=TRUE,catch = FALSE,
                           source = TRUE)
clipped <- l4_clip(gediL4,
                   clip=b_box)

saveRDS(clipped,paste0(outdir,
                       "GEDI.RDS"))

# system2("rsync",
#         c("-avz",
#           "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/GEDI.data/GEDI.RDS",
#           "./outputs/"))

clipped <- readRDS("./outputs/GEDI.RDS")
clipped.filt <- clipped %>%
  filter(l4_quality_flag == 1,
         pft_class == 2,
         degrade_flag == 0,
         sensitivity > 0.9)


# hist(as.numeric(as.Date(clipped.filt$date)))
# unique(clipped.filt$l4_quality_flag)
# unique(clipped.filt$degrade_flag)
# hist(clipped.filt$sensitivity)
# hist(clipped.filt$lat_lowestmode)
# hist(clipped.filt$lon_lowestmode)
# hist(clipped.filt$tree_cover)
# hist(clipped.filt$elev_lowestmode)
# unique(clipped.filt$pft_class)
# hist(clipped.filt$agbd)

ggplot(data = clipped.filt) +
  geom_point(aes(x = lon_lowestmode, y = lat_lowestmode,alpha = agbd,
                 color = as.factor(pft_class))) +
  labs(x = "", y = "") +
  theme_bw() +
  guides(alpha = "none")

ggplot(data = clipped.filt) +
  geom_point(aes(x = lon_lowestmode, y = lat_lowestmode,
                 color = tree_cover)) +
  scale_color_gradient(limits = c(0,100),low = "white",high = "darkgreen",
                       oob = scales::squish) +
  labs(x = "", y = "") +
  theme_bw() +
  guides(alpha = "none")


ggplot(data = clipped.filt) +
  geom_point(aes(x = lon_lowestmode, y = lat_lowestmode,
                 color = agbd)) +
  scale_color_gradient(limits = c(0,500),low = "white",high = "darkgreen",
                       oob = scales::squish) +
  labs(x = "", y = "") +
  theme_bw() +
  guides(alpha = "none")

Delta_m <- 30
Delta = Delta_m*360/(2*pi*6378*1000)  # Conversion 30m into °
lat <- seq(0,0.5,Delta) ; lon <- seq(25,25.5,Delta) ; grid <- expand.grid(lon,lat)
r <- rasterFromXYZ(grid %>%
                     mutate(value = 1))

X <- as.data.frame(clipped.filt %>%
                     dplyr::select(lon_lowestmode,lat_lowestmode,agbd))
Y <- raster(SpatialPixelsDataFrame(points = X[c("lon_lowestmode","lat_lowestmode")],
                       data = X["agbd"],
                       tolerance = 0.716644))

r.rspld <- resample(Y,r)

r.rspld.df <- as.data.frame(r.rspld,xy = TRUE) %>%
  filter(!is.na(agbd))


r.rspld.agg <- raster::aggregate(r.rspld,
                                 fact = 30)
plot(r.rspld.agg)

ggplot(data = r.rspld.df %>%
         filter(!is.na(agbd))) +
  geom_tile(aes(x = x, y = y,
                 color = agbd)) +
  scale_color_gradient(limits = c(0,500),low = "white",high = "darkgreen",
                       oob = scales::squish) +
  labs(x = "", y = "") +
  theme_bw() +
  guides(alpha = "none")

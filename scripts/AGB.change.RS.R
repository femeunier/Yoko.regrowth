rm(list = ls())

library(raster)
library(dplyr)
library(tidyr)
library(ggplot2)

# files <- c("/home/femeunier/Downloads/AGB2007c_1km.tif",
#            "/home/femeunier/Downloads/AGB2008c_1km.tif",
#            "/home/femeunier/Downloads/AGB2010c_1km.tif",
#            "/home/femeunier/Downloads/AGB2017c_1km.tif")
# years <- c(2007,2008,2010,2017)

files <- c("/home/femeunier/Downloads/AGB2007 (1).tif")
years <- c(2007)

e <- as(extent(0, 2.5, 0, 2.5), 'SpatialPolygons')
crs(e) <- "+proj=longlat +datum=WGS84 +no_defs"

df.all <- data.frame()

for (ifile in seq(1,length(files))){
  rst <- raster(files[ifile])
  rst.rspld <- crop(rst, e)

  cdf <- as.data.frame(rst.rspld,
                       xy = TRUE) %>%
    rename(lon = x,
           lat = y) %>%
    mutate(year = years[ifile])

  df.all <- bind_rows(list(df.all,
                           cdf))
}

df.wide <- df.all %>%
  pivot_wider(names_from = "year",
              values_from = "agb") %>%
  mutate(diff = `2017` - `2007`)


coords2keep <- df.wide %>%
  filter(diff < -50) %>%
  mutate(lon.lat = paste(lon,lat,sep = "_"))


df2plot <- df.all %>%
  mutate(lon.lat = paste(lon,lat,sep = "_")) %>%
  dplyr::filter(lon.lat %in% c(coords2keep %>% pull(lon.lat)))

ggplot(data = df2plot) +
  geom_line(aes(x = year,
                y = agb,
                group = interaction(lat,lon))) +
  theme_bw()


df2plot.wide <- df2plot %>%
  pivot_wider(names_from = "year",
              values_from = "agb") %>%
  mutate(diff = `2017` - `2008`)

ggplot(data = df2plot.wide) +
  geom_density(aes(x = diff/9/20)) +  # kgC/m²/yr
  theme_bw()

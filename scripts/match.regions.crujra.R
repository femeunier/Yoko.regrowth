rm(list = ls())

library(raster)

sites <- read.csv("./data/sites.locations.csv") %>%
  rename(lon = lon.m,
         lat = lat.m) %>%
  na.omit() %>%
  dplyr::select(-X) %>%
  mutate(site = 1:n())

map <- readRDS("/home/femeunier/Documents/projects/Congo.ED2/outputs/biome.CRUJRA.1901.2019.RDS") %>%
  filter(model == "ORCHIDEE")

crujra.coord <- data.frame()
for (isite in seq(1,nrow(sites))){
  clat <- sites[isite,"lat"]; clon <- sites[isite,"lon"]
  cdist <- map %>%
    mutate(dist = sqrt((lat - clat)**2 + (lon - clon)**2)) %>%
    arrange(dist) %>%
    slice_head(n = 1)

  crujra.coord <- bind_rows(crujra.coord,
                            cdist %>%
                              dplyr::select(lon,lat,dist) %>%
                              rename(crujra.lon = lon,
                                     crujra.lat = lat) %>%
                              mutate(site = isite))
}

sites.crujra <- sites %>%
  left_join(crujra.coord,
            by = "site")

selected <- c("DRC_Yoko","DRC_Baego","IC_Yaya")

write.csv(sites.crujra %>%
            # dplyr::filter(region %in% selected) %>%
            dplyr::select(region,crujra.lon,crujra.lat) %>%
            rename(site = region,
                   lon = crujra.lon,
                   lat = crujra.lat),
          "./data/data2share/sites.locations.csv")

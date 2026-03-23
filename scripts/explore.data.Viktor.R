rm(list = ls())

library(dplyr)
library(ggplot2)
library(ggthemes)
library(raster)
library(sf)

file.data <- "./data/Meta-analysis_plots_FINAL_10_cm_Felicien.csv"
data <- read.csv(file.data,sep = ";") %>%
  mutate(African.region = case_when(country %in% c("DRC","Cameroon","Gabon") ~ "Central Africa",
                                    country %in% c("Kenya","Tanzania","Mozambique") ~ "Eastern Africa",
                                    country %in% c("Ghana","Ivory Coast","Liberia",
                                                   "Republic of Guinea","Sierra Leone") ~ "Western Africa",
                                    country %in% c("Angola","Zambia") ~ "Southern Africa",
                                    TRUE ~ NA_character_))

data %>%
  filter(region == "DRC_Babagulu")

data %>%
  filter(country %in% c("Gabon","Cameroon")) %>%
  pull(lon)

data.loc <- data %>%
  group_by(region) %>%
  summarise(lon.m = mean(lon,na.rm = TRUE),
            lat.m = mean(lat,na.rm = TRUE),
            .groups = "keep") %>%
  mutate(lat.m = case_when(region == "DRC_Babagulu" ~ 0.484091238,
                           TRUE ~ lat.m),
         lon.m = case_when(region == "DRC_Babagulu" ~ 25.56989943,
                           TRUE ~ lon.m))

write.csv(data.loc,
          "./data/sites.locations.csv")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

r <- raster("~/Documents/projects/LianaRemovalRevisited/data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1_pantropical_aggr.tif")
df.r <- as.data.frame(r,xy = TRUE) %>%
  rename(lon =  x,
         lat = y,
         LU = C3S.LC.L4.LCCS.Map.300m.P1Y.2020.v2.1.1_pantropical_aggr) %>%
  filter(!is.na(LU))


# shape <- readOGR(dsn = "~/Downloads/", layer = "plot_corners_polygons_32635")
# shape <- st_read("~/Downloads/plot_corners_polygons_32635.shp")

ggplot() +
  geom_raster(data = df.r,
              aes(x = lon, y = lat, fill = as.factor(LU)),
              alpha = 0.4,show.legend = FALSE) +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = data ,
             aes(x = lon, y = lat,
                 col = region), shape = 16, size = 0.1) +
  geom_point(data = data.loc,
             aes(x = lon.m, y = lat.m,
                 col = region), size = 1, shape = 0) +
  theme_map() +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        text = element_text(size = 24)) +
  scale_fill_manual(values = c("white",c("#72a83d"),"darkgreen")) +
  scale_x_continuous(limits = c(-20,60),expand = c(0,0)) +
  scale_y_continuous(limits = c(-25,25),expand = c(0,0))


ggplot() +
  geom_raster(data = df.r,
              aes(x = lon, y = lat, fill = as.factor(LU)),
              alpha = 0.4,show.legend = FALSE) +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = data.loc,
             aes(x = lon.m, y = lat.m), color = "grey17",
             fill = "grey17",alpha = 0.8,
             size = 2, shape = 16) +
  theme_map() +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        text = element_text(size = 24)) +
  scale_fill_manual(values = c("white",c("#72a83d"),"darkgreen")) +
  scale_x_continuous(limits = c(-20,60),expand = c(0,0)) +
  scale_y_continuous(limits = c(-25,25),expand = c(0,0))


other.data <- read.csv("~/Downloads/PlotCoordinates_Baego_Yoko_Ituri-meta-1.csv")


ggplot() +
  geom_raster(data = df.r,
              aes(x = lon, y = lat, fill = as.factor(LU)),
              alpha = 0.4,show.legend = FALSE) +
  geom_sf(data = world,fill = NA,color = "grey") +
  geom_point(data = data %>%
               filter(grepl("DRC",region)),
             aes(x = lon, y = lat,
                 col = region), shape = 16, size = 0.1) +

  geom_point(data = other.data,
             aes(x = longitude_degrees, y = latitude_degrees), shape = 16, size = 2) +
  geom_point(data = data.loc %>%
               filter(grepl("DRC",region)),
             aes(x = lon.m, y = lat.m,
                 col = region), size = 1, shape = 0) +
  # theme_map() +
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 24)) +
  geom_sf(data = shape %>%
            filter(grepl("IITA",region)),
          colour='red', fill='red') +
  scale_fill_manual(values = c("white",c("#72a83d"),"darkgreen"))+
  scale_x_continuous(limits = c(25,30),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2),expand = c(0,0))

plot(shape)

data.sum <- data %>%
  group_by(region,year_inventory,age) %>%
  summarise(AGC_ha.m = mean(AGC_ha,na.rm = TRUE),
            AGC_ha.sd = sd(AGC_ha,na.rm = TRUE),
            .groups = "keep")

saveRDS(data.sum,
        "./data/Data.sum.RDS")

ggplot(data = data.sum,
       aes(x = age, y = AGC_ha.m/10, color = region)) +
  geom_point() + # kgC/m²
  scale_x_log10() +
  scale_y_log10() +
  # stat_smooth(method = "nls",se = FALSE,formula = y ~ a*(1 -exp(-b*x))) +
  theme_bw() +
  theme(legend.position = "none")

selected <- c("DRC_Yoko","DRC_Baego","IC_Yaya")
unique(data.sum$region)


ggplot(data = data.sum %>%
         filter(region %in% selected),
       aes(x = age, y = AGC_ha.m/10, color = region)) +
  geom_point() + # kgC/m²
  # scale_x_log10() +
  # scale_y_log10() +
  # stat_smooth(method = "nls",se = FALSE,formula = y ~ a*(1 -exp(-b*x))) +
  theme_bw()

################################################################################
data.numeric <- data %>%
  mutate(year_inventory =
           case_when(year_inventory == "Oct 2016 - May 2017" ~ "2017",
                     year_inventory == "Oct 2016 - Oct 2020" ~ "2020",
                     year_inventory == "" ~ "2024",
                     TRUE ~ year_inventory)) %>%
  mutate(year_inventory_num = as.numeric(year_inventory))

data.oldgrowth <- data.numeric %>%
  filter(structure != "P") %>%
  dplyr::select(African.region,country,region,plotID,age,year_inventory_num) %>%
  left_join(data.numeric %>%
              dplyr::filter(structure == "P") %>%
              mutate(old.growth = TRUE) %>%
              dplyr::select(region,old.growth) %>%
              distinct(),
            by = "region")

ggplot(data = data.numeric %>%
         filter(structure != "P") %>%
         group_by(region) %>%
         arrange(desc(age))) +
  geom_segment(aes(x = year_inventory_num - age, xend = year_inventory_num,
                   y = paste(region,age,sep = "_"),
                   color = region)) +
  geom_point(data = data.oldgrowth %>%
               filter(old.growth),
             aes(x = 2030, y = paste(region,age,sep = "_"),
                 color = region)) +
  theme_bw() +
  facet_wrap(~ African.region, scales = "free_y") +
  theme(legend.position = "none")



ggplot(data = data.numeric %>%
         filter(structure != "P",
                region == "DRC_Yoko") %>%
         group_by(region) %>%
         arrange(desc(age)) %>%
         mutate(region_age = factor(paste0(region,"_",age),
                                    levels = rev(c("DRC_Yoko_5","DRC_Yoko_12","DRC_Yoko_20","DRC_Yoko_60"))))) +
  geom_segment(aes(x = year_inventory_num - age, xend = year_inventory_num,
                   y = region_age,
                   color = region)) +
  geom_vline(aes(xintercept = 2019,
                 color = region),
             linetype = 2) +
  theme_bw() +
  # facet_wrap(~ African.region, scales = "free_y") +
  theme(legend.position = "none") +
  labs(x = "",y = "")


data.numeric %>%
  filter(region == "Gh_Abofour") %>%
  pull(age) %>%
  sort()


ggplot(data = data.numeric %>%
         filter(structure != "P",
                region %in% selected) %>%
         group_by(region) %>%
         arrange(desc(age))) +
  geom_segment(aes(x = year_inventory_num - age, xend = year_inventory_num,
                   y = paste(region,age,sep = "_"),
                   color = region)) +
  geom_point(data = data.oldgrowth %>%
               filter(old.growth,
                      region %in% selected),
             aes(x = 2030, y = paste(region,age,sep = "_"),
                 color = region)) +
  theme_bw() +
  facet_wrap(~ African.region, scales = "free_y") +
  theme(legend.position = "none")


plot.ages <- data.numeric %>%
  filter(structure != "P") %>%
  group_by(region, age) %>%
  summarise(time_ini = unique(year_inventory_num - age),
            time_end = unique(year_inventory_num),
            N = length(unique(year_inventory_num)),
            .groups = "keep") %>%
  rename(site = region)

write.csv(plot.ages,
          "./data/data2share/plot_ages.csv")

unique(plot.ages$site)
nrow(plot.ages)

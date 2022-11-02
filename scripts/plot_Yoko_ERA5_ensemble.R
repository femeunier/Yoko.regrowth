# rm(list = ls())

library(xts)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

met.Yoko <- readRDS(file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko",
                                     paste0("ERA5_Yoko_ensemble_processed"),
                                     "ERA5_ensemble.RDS"))

for (i in seq(1,10)){
  met.Yoko[[i]] <- as.data.frame(cbind(met.Yoko[[i]],ensemble.member = i)) %>%
    tibble::rownames_to_column(var = "t")
}


df.met.Yoko <- do.call(rbind,met.Yoko) %>%
  mutate(year = year(t),
         month = month(t),
         day = day(t),
         h = hour(t),
         min = minute(t),
         sec = second(t))

# Conversion
df.met.Yoko.conv <- df.met.Yoko %>%
  mutate(ssrd = ssrd/(3*3600),
         strd = strd/(3*3600),
         tp = tp*1000/(3*3600),
         temp = t2m - 273.15,
         dewpoint = d2m - 273.15,
         beta = (112 - (0.1 * temp) + dewpoint) / (112 + (0.9 * temp)),
         rh = beta ^ 8) %>%
  mutate(sh = PEcAn.data.atmosphere::rh2qair(rh,
                                             as.numeric(t2m),
                                             as.numeric(sp))) %>%
  dplyr::select(-c(t2m,d2m,beta,rh,dewpoint))

saveRDS(df.met.Yoko.conv,
        "/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/ERA5_Yoko_ensemble_processed/ERA5_ensemble_conv.RDS")

variables <- tibble::tribble(
  ~description, ~units, ~api_name, ~xts_name,
  "air_temperature", "Kelvin", "2m_temperature", "temp",
  "air_pressure", "Pa", "surface_pressure", "sp",
  "specific_humidity", "g/g", "2m_specific_humidity", "sh",
  "precipitation_flux", "kg/m2/s", "total_precipitation", "tp",
  "eastward_wind", "m/s", "10m_u_component_of_wind", "u10",
  "northward_wind", "m/s", "10m_v_component_of_wind", "v10",
  "surface_downwelling_shortwave_flux_in_air", "W/m2", "surface_solar_radiation_downwards", "ssrd",
  "surface_downwelling_longwave_flux_in_air", "W/m2", "surface_thermal_radiation_downwards", "strd",
  "air_mole_fraction","mol/mol","mole_fraction_of_carbon_dioxide_in_air","mole_fraction_of_carbon_dioxide_in_air")


df.all.long <- df.met.Yoko.conv %>%
  pivot_longer(cols = -c(t,ensemble.member,year,month,day,h,min,sec),
               names_to = "var",
               values_to = "value") %>% arrange(t)

l.cycle <- df.all.long %>%
  dplyr::filter(var %in% c("mole_fraction_of_carbon_dioxide_in_air","sh",
                           "ssrd","strd","temp","tp")) %>%
  group_by(year,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(var,ensemble.member) %>%
  mutate(value.m.smooth = rollapply(value.m, width = 10, FUN = mean, align = "center", partial = TRUE, na.rm = TRUE)) %>%
  group_by(var,year) %>%
  mutate(value.m.smooth.m = mean(value.m.smooth)) %>%
  ungroup()

ggplot(l.cycle) +
  geom_rect(xmin = 1960, xmax = 1969 + 11/12, ymin = -Inf, ymax = Inf,
            fill = "grey",
            alpha = 0.4) +
  geom_line(aes(x = year,y = value.m, group = ensemble.member)) +
  geom_line(aes(x = year,value.m.smooth.m),
            color = "red") +
  facet_wrap(~var,scales = "free_y") +
  scale_x_continuous(breaks = seq(1850,2020,50)) +
  labs(x = "",y = "") +
  theme_bw()

s.cycle <- df.all.long %>%
  dplyr::filter(var %in% c("sh","ssrd","strd","temp",
                           "tp","v10")) %>%
  group_by(month,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

s.cycle.sum <- s.cycle %>%
  group_by(month,var) %>%
  summarise(value.m = mean(value.m),
            .groups = "keep")

ggplot(s.cycle) +
  geom_line(aes(x = month,y = value.m, group = ensemble.member)) +
  geom_line(data = s.cycle.sum,
            aes(x = month,y = value.m), color = "red") +

  facet_wrap(~var,scales = "free_y") +
  labs(x = "",y = "") +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J","F","M","A","M","J",
                                "J","A","S","O","N","D")) +
  theme_bw()

d.cycle <- df.all.long %>%
  dplyr::filter(var %in% c("ssrd","strd","sh","temp")) %>%
  group_by(h,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

d.cycle.sum <- d.cycle %>%
  group_by(h,var) %>%
  summarise(value.m = mean(value.m),
            .groups = "keep")

ggplot(d.cycle) +
  geom_line(aes(x = h,y = value.m, group = ensemble.member)) +
  geom_line(data = d.cycle.sum,
            aes(x = h,y = value.m), color = "red") +
  facet_wrap(~var,scales = "free_y",nrow = 1) +
  labs(x = "",y = "") +
  theme_bw()


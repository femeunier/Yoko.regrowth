rm(list = ls())

library(ncdf4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lubridate)


years <- 2001:2010

main.dir <- "/home/femeunier/Documents/projects/YGB/data/CRUNCEP/Yoko"

vars <- c("air_temperature",
          "surface_downwelling_longwave_flux_in_air",
          "air_pressure",
          "surface_downwelling_shortwave_flux_in_air",
          "eastward_wind",
          "northward_wind",
          "specific_humidity",
          "precipitation_flux")
var.names <- c("temperature","longwave","press","shortwave","ewind","nwind","humidity","precip")

df.precip.all <- data.frame()

df.all <- data.frame()
for (iyear in seq(1,length(years))){
  print(years[iyear])

  cruncep.file <- file.path(main.dir,paste0("CRUNCEP.",years[iyear],".nc"))
  ncfile <- nc_open(cruncep.file)

  for (ivar in seq(1,length(vars))){

    if(ivar == 1) {

      temp.time <- ncvar_get(ncfile,"time")
      ctime <- temp.time*86400

      cdf <- data.frame(time = ctime)
    }

    cvar <- ncvar_get(ncfile,vars[ivar])
    cdf <- cdf %>% mutate(value = cvar)

    colnames(cdf)[ivar + 1] <- var.names[ivar]
  }

  nc_close(ncfile)

  df.all <- bind_rows(list(df.all,
                           cdf %>% mutate(origin.yr = years[iyear])))


}

df.all.time <- df.all %>% mutate(t = as.POSIXct(time,
                                                tz = "GMT",
                                                origin = paste0(origin.yr,"-01-01"))) %>%
  mutate(year = year(t),
         month = month(t),
         day = day(t),
         h = hour(t)) %>%
  dplyr::select(-c(time,origin.yr))

df.all.long <- df.all.time %>%
  pivot_longer(cols = -c(t,year,month,day,h),
               names_to = "var",
               values_to = "value")

ggplot(df.all.long) +
  geom_line(aes(x = t,y = value)) +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  theme_bw()

s.cycle <- df.all.long %>%
  group_by(month,var) %>%
  summarise(value.m = mean(value),
            .groups = "keep")

ggplot(s.cycle) +
  geom_line(aes(x = month,y = value.m)) +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  theme_bw()

d.cycle <- df.all.long %>%
  group_by(h,var) %>%
  summarise(value.m = mean(value),
            .groups = "keep")

ggplot(d.cycle) +
  geom_line(aes(x = h,y = value.m)) +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  theme_bw()



# saveRDS(df.precip.all,"CRUNCEP.RDS")
# scp /home/femeunier/Documents/projects/SoilSensitivity/scripts/analyze_CRUNCEP.R hpc:/data/gent/vo/000/gvo00074/felicien/R

# world <- ne_countries(scale = "medium", returnclass = "sf")

# ggplot(data = world) +
#   geom_sf() +
#   geom_tile(data = df.precip %>% filter(!is.na(precip)),
#             aes(x = lon, y = lat,fill = precip),na.rm = TRUE) +
#   coord_sf(xlim = c(-90, -30), ylim = c(-20, 15), expand = FALSE)

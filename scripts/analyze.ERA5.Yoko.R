rm(list = ls())

library(ncdf4)
library(reshape2)
library(dplyr)
library(ggplot2)
library(lubridate)
library(YGB)
library(zoo)
library(tidyr)
library(PEcAn.ED2)

years <- c(1960:2021)                                                          # years available
main.dir <- "/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko"  # location of the file

prefix <- "ERA5_Yoko_h__"

vars <- c("t2m","sp","tp","u10","v10","ssrd","strd","d2m")

var.names <- vars

df.all <- data.frame()
for (iyear in seq(1,length(years))){
  print(years[iyear])
  ERA5.file <- file.path(main.dir,paste0(prefix,years[iyear],".nc"))

  if (!file.exists(ERA5.file)) {
    stop(paste("Missing file:",ERA5.file))
    next()
  }

  ncfile <- nc_open(ERA5.file)

  for (ivar in seq(1,length(vars))){
    if(ivar == 1) {
      temp.time <- ncvar_get(ncfile,"time")
      ctime <- temp.time/24
      cdf <- data.frame(time = ctime)
    }
    cvar <- apply(ncvar_get(ncfile,vars[ivar]),3,mean)   # 4 sites (2lon, 2lat) --> we take the mean
    cdf <- cdf %>% mutate(value = cvar)
    colnames(cdf)[ivar + 1] <- var.names[ivar]
  }
  nc_close(ncfile)
  df.all <- bind_rows(list(df.all,
                           cdf))
}

df.all.time <- df.all %>% mutate(t = as.POSIXct(time*86400,
                                                tz = "UTC",
                                                origin = "1900-01-01")) %>%
  mutate(year = year(t),
         month = month(t),
         day = day(t),
         h = hour(t),
         min = minute(t),
         sec = second(t))



years2add <- 1849:1959                   # we complete the historical period ... (+ 1 year to fix the timezone issue)
years2replace <- seq(1960,1969)          # ... with the first decade of observations (no de-trend at the moment)

replacement = TRUE

df.all.time.all <- df.all.time
for (iyear in seq(1,length(years2add))){

  if (lubridate::leap_year(years2add[iyear])){
    yearselect <- sample(rep(years2replace[lubridate::leap_year(years2replace)],2),1)
  } else {
    yearselect <- sample(rep(years2replace[!lubridate::leap_year(years2replace)],2),1)
  }

  if (!replacement) years2replace <- years2replace[years2replace != yearselect]

  year2replace <- df.all.time %>% filter(year == yearselect)

  df.all.time.all <- bind_rows(list(df.all.time.all,
                                    year2replace %>% mutate(year = years2add[iyear],
                                                            t = as.POSIXct(paste0(years2add[iyear],"-",sprintf("%02d",year2replace$month),"-",sprintf("%02d",year2replace$day)," ",
                                                                                  sprintf("%02d",year2replace$h),":",sprintf("%02d",year2replace$min),":",sprintf("%02d",year2replace$sec)),
                                                                           tz = "GMT",
                                                                           format = "%Y-%m-%d %H:%M:%S"))
    ))

}

df.all.time.all.local.time <- df.all.time.all %>% mutate(t = t + 3600) %>%
  mutate(year = year(t),
         month = month(t),
         day = day(t),
         h = hour(t),
         min = minute(t),
         sec = second(t)) %>%
  dplyr::filter(year > min(year))


# Source = https://keelingcurve.ucsd.edu/permissions-and-data-sources/
CO2.data <- read.csv(file = "./data/CO2.csv",header = TRUE) %>%
  mutate(year = floor(time)) %>%
  mutate(CO2 = case_when(CO2 == -99.99 ~ NA_real_,
                         TRUE ~ CO2)) %>%
  group_by(year) %>%
  summarise(CO2 = mean(CO2,na.rm = TRUE)/1e6,
            .groups = "keep")
CO2.data.all <- data.frame(year = 1700:2022) %>%
  left_join(CO2.data,
            by = "year") %>%
  mutate(CO2.interp = na.approx(CO2, maxgap = 50, rule = 2)) %>%
  dplyr::select(-CO2) %>%
  rename(CO2 = CO2.interp)

df.all.time.CO2 <- df.all.time.all.local.time %>% left_join(CO2.data.all,
                                                            by = "year") %>%
  filter(!is.na(t))


df.all.long <- df.all.time.CO2 %>%
  pivot_longer(cols = -c(time,t,year,month,day,h,min,sec),
               names_to = "var",
               values_to = "value") %>% arrange(t)


l.cycle <- df.all.long %>%
  group_by(year,var) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(var) %>%
  mutate(value.m.smooth = rollapply(value.m, width = 10, FUN = mean, align = "center", partial = TRUE)) %>%
  ungroup()

ggplot(l.cycle) +
  geom_line(aes(x = year,y = value.m)) +
  geom_line(aes(x = year,value.m.smooth),
            color = "red") +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  theme_bw()

s.cycle <- df.all.long %>%
  group_by(month,var) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

ggplot(s.cycle) +
  geom_line(aes(x = month,y = value.m)) +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  scale_x_continuous(breaks = 1:12) +
  theme_bw()

d.cycle <- df.all.long %>%
  group_by(h,var) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

ggplot(d.cycle) +
  geom_line(aes(x = h,y = value.m)) +
  facet_wrap(~var,scales = "free") +
  labs(x = "",y = "") +
  theme_bw()


###################################################################################

start_year <- 1850 ; end_year <- 2021

start_date <- paste0(start_year,"-01-01")
end_date <- paste0(end_year,"-12-31")
slat = 0.3
slon = 25.3
overwrite = TRUE

colnames(df.all.time.CO2)[colnames(df.all.time.CO2) == "CO2"] <- "mole_fraction_of_carbon_dioxide_in_air"
vars <- unique(c(vars,"mole_fraction_of_carbon_dioxide_in_air"))

local.data.point <- vars %>%
  purrr::map_dfc(function(vname) {

    nn <- df.all.time.CO2[,vname]

    if (!is.numeric(nn)) {
      PEcAn.logger::logger.severe(paste0(
        "Expected raster object to be numeric, but it has type `",
        paste0(typeof(nn), collapse = " "),
        "`"
      ))
    }

    nn

  }) %>%
  `colnames<-`(vars)

site.XTS <- xts::xts(local.data.point, order.by = (df.all.time.CO2$t))
site.XTS[site.XTS<0] <- 0

destination <- file.path(main.dir,paste0("ERA5_Yoko_processed"))
dir.create(destination,showWarnings = FALSE)
saveRDS(object = site.XTS,
        file = file.path(destination,"ERA5_reanalysis.RDS"))

system2("scp",paste(file.path(destination,"ERA5_reanalysis.RDS"),
                    "hpc:/data/gent/vo/000/gvo00074/felicien/R/data/"))

stop()

for (year in seq(start_year,end_year)){

  print(year)

  out <- Yoko.regrowth::met2CF.ERA5(
    lat = slat,
    long = slon,
    start_date =  paste0(year,"-01-01"),
    end_date = paste0(year,"-12-31"),
    sitename = "Yoko",
    outfolder = destination,
    out.xts = list(site.XTS %>% subset(year(index(.)) == year)),
    overwrite = overwrite,
    verbose = FALSE,
    type = "reanalysis")
}

for (year in seq(start_year,end_year)){

  print(year)
  Yoko.regrowth::mod.met2model.ED2(
    in.path =  "/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/ERA5_Yoko_processed/ERA5_Yoko_1/",
    in.prefix = "ERA5.1",
    outfolder = "/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/ERA5_Yoko_processed/ERA5_Yoko_1/ED2",
    start_date =  paste0(year,"-01-01"),
    end_date = paste0(year,"-12-31"),
    lat = slat,
    lon = slon,
    lst = 0,
    overwrite = overwrite)
}

# nc <- nc_open("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/ERA5_Yoko_processed/ERA5_Yoko_1/ED2/1700APR.h5")
# # plot(ncvar_get(nc,"surface_downwelling_longwave_flux_in_air"),type = "l")
#
# plot(ncvar_get(nc,"prate"),type = "l")

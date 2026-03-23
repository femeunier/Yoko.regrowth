rm(list = ls())

library(lubridate)
library(tidyr)
library(dplyr)
library(CFtime)
library(ggplot2)

system2("rsync",
        c("hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/Timeseries.regrowth.RDS",
          "./outputs/"))

sites <- read.csv("/home/femeunier/Documents/projects/Yoko.regrowth//data/data2share/sites.locations.csv") %>%
  mutate(lon_lat = paste0(lon,"_",lat))

selected <- c("DRC_Yoko","DRC_Baego","IC_Yaya")
# selected <- c("DRC_Yoko")
selected <- sites$site
# selected <- c("DRC_Babagulu","DRC_Baego","DRC_Djolu","DRC_Sakania",
#               "DRC_Salonga","DRC_Yoko","GU_Ziama","IC_Agbo1","IC_Badenou",
#               "IC_Foumbou","IC_Irobo","IC_Tene","IC_Yaya","Ken_Nyatike",
#               "SL_Outamba","ZAM_Chintumukulu","Gh_Abofour")


Timeseries.regrowth.all <- readRDS("/Users/felicien/Documents/projects/Yoko.regrowth//outputs/Timeseries.regrowth.RDS") %>%
  ungroup() %>%
  mutate(year = year(time),
         month = month(time),
         hour = hour(time)) %>%
  mutate(lon_lat = paste0(lon,"_",lat))

Timeseries.regrowth <- Timeseries.regrowth.all %>%
  dplyr::filter(lon_lat %in%
                  c(sites %>%
                      filter(site %in% selected) %>%
                      pull(lon_lat))) %>%
  left_join(sites %>%
              dplyr::select(lon_lat,site),
            by = "lon_lat") %>%
  dplyr::select(-lon_lat)

Timeseries.regrowth.long <- Timeseries.regrowth %>%
  pivot_longer(cols = -c(lon,lat,time,site,year,month,hour),
               names_to = "variable") %>%
  group_by(variable,site) %>%
  arrange(time)

Timeseries.regrowth.long.fixed <- Timeseries.regrowth.long %>%
  mutate(value = case_when(variable == "tmin" & value <= 200 ~ NA,
                           TRUE ~ value)) %>%
  ungroup() %>%
  mutate(ID = 1:n())

# Timeseries.regrowth.long.fixed <- Timeseries.regrowth.long.fixed %>%
#   group_by(site,variable,month,hour) %>%
#   mutate(value = case_when(is.na(value) ~ mean(value[year %in% c(year + (-15:14))],
#                                                na.rm = TRUE),
#                            TRUE ~ value))

N.na <- Timeseries.regrowth.long.fixed %>%
  filter(is.na(value))

for (i in seq_len(nrow(N.na))){

  print(i/nrow(N.na))

  crow <- N.na[i,]
  cID <- crow$ID

  Timeseries.regrowth.long.fixed[cID,"value"] <-
    Timeseries.regrowth.long.fixed %>%
    filter(variable == crow$variable,
           month == crow$month,
           year == c(crow$year) + (-15:14),
           hour == crow$hour,
           site == crow$site) %>%
    filter(!is.na(value)) %>%
    pull(value) %>%
    mean()

}


DC <- Timeseries.regrowth.long.fixed %>%
  group_by(lon,lat,site,hour,variable) %>%
  summarise(value.m = mean(value),
            .groups = "keep")

ggplot(data = DC) +
  geom_line(aes(x = hour, y = value.m,
                color = site)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "none")

SC <- Timeseries.regrowth.long.fixed %>%
  group_by(lon,lat,site,month,variable) %>%
  summarise(value.m = mean(value),
            .groups = "keep")

ggplot(data = SC %>%
         filter(variable %in% c("dswrf","pre","tmp"))) +
  geom_line(aes(x = month, y = value.m,color = site)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  scale_x_continuous(breaks = c(1:12),
                     labels = c('J',"F","M","A","M","J",
                                "J","A","S","O","N","D")) +
  theme(text = element_text(size = 16)) +
  labs(x = "", y = "value", color = "")

 TS <- Timeseries.regrowth.long.fixed %>%
  group_by(lon,lat,site,year,month,variable) %>%
  summarise(value.m = mean(value),
            .groups = "keep")

ggplot(data = TS) +
  geom_line(aes(x = year + (month-1/2)/12, y = value.m,
                color = site)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "none")

Timeseries.regrowth.wide.fixed <-
  Timeseries.regrowth.long.fixed %>%
  dplyr::select(-c(ID,year,month,hour)) %>%
  pivot_wider(names_from = variable,
              values_from = value) %>%
  arrange(site,time)

ggplot(data = Timeseries.regrowth.long.fixed %>%
         filter(variable %in% c("ugrd","vgrd"))) +
  geom_density(aes(x = value, fill = site), alpha = 0.4, color = NA) +
  facet_wrap(~ variable) +
  theme_bw()

# summary(Timeseries.regrowth.long.fixed %>%
#           filter(variable %in% c("vgrd")) %>% pull(value))


for (csite in selected){

  print(csite)
  cdf <- Timeseries.regrowth.wide.fixed %>%
    filter(site == csite)

  saveRDS(cdf,
          paste0("/home/femeunier/Documents/projects/Yoko.regrowth/data/data2share/drivers/",
                 "climate.",csite,".RDS"))

}


unique(Timeseries.regrowth.wide.fixed$site)

A <- Timeseries.regrowth.wide.fixed %>%
  mutate(year = year(time),
         month = month(time)) %>%
  filter(year %in% c(1971:2000)) %>%
  rename(region = site) %>%
  group_by(region,year,month) %>%
  summarise(pre = mean(pre),
            tmp = mean(tmp),
            tmin = mean(tmin),
            tmax = mean(tmax),
            dswrf = mean(dswrf),
            spfh = mean(spfh),
            .groups = "keep") %>%
  group_by(region,month) %>%
  summarise(pre = mean(pre),
            tmp = mean(tmp),
            tmin = mean(tmin),
            tmax = mean(tmax),
            dswrf = mean(dswrf),
            spfh = mean(spfh),
            .groups = "keep")

saveRDS(A,"./outputs/Climate.sites.RDS")

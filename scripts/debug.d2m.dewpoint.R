rm(list = ls())

library(xts)
library(dplyr)
library(lubridate)

# reanalysis
met.Yoko <- readRDS(file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko",
                                     paste0("ERA5_Yoko_processed"),
                                     "ERA5_reanalysis.RDS"))

df.met.Yoko <- as.data.frame(met.Yoko) %>%
  tibble::rownames_to_column(var = "t") %>%
  mutate(year = year(t),
         month = month(t),
         day = day(t),
         h = hour(t),
         min = minute(t),
         sec = second(t))

df.met.Yoko.conv <- df.met.Yoko %>%
  mutate(ssrd = ssrd/(1*3600),
         strd = strd/(1*3600),
         tp = tp*1000/(1*3600),
         temp = t2m - 273.15,
         dewpoint = d2m - 273.15,
         beta = (112 - (0.1 * temp) + dewpoint) / (112 + (0.9 * temp)),
         rh = beta ^ 8) %>%
  mutate(sh = PEcAn.data.atmosphere::rh2qair(rh,
                                             as.numeric(t2m),
                                             as.numeric(sp)))


# Ensemble
met.Yoko <- readRDS(file = file.path("/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko",
                                     paste0("ERA5_Yoko_ensemble_processed"),
                                     "ERA5_ensemble.RDS"))

# Renalaysis
df.met.Yoko.conv.reanalysis <- readRDS(file = "/home/femeunier/Documents/projects/Yoko.regrowth/data/ERA5/Yoko/ERA5_Yoko_processed/ERA5_reanalysis_conv.RDS") %>%
  filter(year <= 2021)


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
df.met.Yoko.conv.ensemble <- df.met.Yoko %>%
  mutate(ssrd = ssrd/(3*3600),
         strd = strd/(3*3600),
         tp = tp*1000/(3*3600),
         temp = t2m - 273.15,
         dewpoint = d2m - 273.15,
         beta = (112 - (0.1 * temp) + dewpoint) / (112 + (0.9 * temp)),

         rh = beta ^ 8) %>%
  mutate(sh = PEcAn.data.atmosphere::rh2qair(rh,
                                             as.numeric(t2m),
                                             as.numeric(sp)))

################################################################################
# Comparison
df.all.long <- df.met.Yoko.conv %>%
  dplyr::select(t,year,month,h,temp,dewpoint,beta,rh,sh) %>%
  pivot_longer(cols = -c(t,year,month,h),
               names_to = "var",
               values_to = "value") %>% arrange(t)

df.all.long.ensemble <- df.met.Yoko.conv.ensemble %>%
  dplyr::select(t,ensemble.member,year,month,h,temp,dewpoint,beta,rh,sh) %>%
  pivot_longer(cols = -c(t,ensemble.member,year,month,h),
               names_to = "var",
               values_to = "value") %>% arrange(t)

#######################################################
# Diel cycle

d.cycle <- df.all.long.ensemble %>%
  group_by(h,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

d.cycle.sum <- d.cycle %>%
  group_by(h,var) %>%
  summarise(value.m = mean(value.m),
            .groups = "keep")

d.cycle.renanalysis <- df.all.long %>%
  group_by(h,var) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

ggplot(d.cycle) +
  geom_line(aes(x = h,y = value.m, group = ensemble.member)) +
  geom_line(data = d.cycle.sum,
            aes(x = h,y = value.m), color = "red") +
  geom_line(data = d.cycle.renanalysis,
            aes(x = h,y = value.m), color = "red",linetype = 2) +
  facet_wrap(~var,scales = "free_y",nrow = 1) +
  labs(x = "",y = "") +
  theme_bw()

##############################################################
# timeseries


l.cycle <- df.all.long.ensemble %>%
  dplyr::filter(var %in% c("dewpoint")) %>%
  group_by(year,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  group_by(var,ensemble.member) %>%
  mutate(value.m.smooth = rollapply(value.m, width = 10, FUN = mean, align = "center", partial = TRUE, na.rm = TRUE)) %>%
  group_by(var,year) %>%
  mutate(value.m.smooth.m = mean(value.m.smooth)) %>%
  ungroup()

l.cycle.sum <- l.cycle %>%
  group_by(year,var) %>%
  summarise(value.m = mean(value.m),
            .groups = "keep")

l.cycle.reanalysis <- df.all.long %>%
  dplyr::filter(var %in% c("dewpoint")) %>%
  group_by(year,var) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep") %>%
  ungroup()


ggplot(l.cycle) +
  geom_rect(xmin = 1960, xmax = 1969 + 11/12, ymin = -Inf, ymax = Inf,
            fill = "grey",
            alpha = 0.4) +
  geom_line(aes(x = year,y = value.m, group = ensemble.member)) +
  geom_line(data = l.cycle.sum,
            aes(x = year,value.m),
            color = "red") +
  geom_line(data = l.cycle.reanalysis,
            aes(x = year,value.m),
            color = "red",linetype = 2) +
  facet_wrap(~var,scales = "free_y") +
  scale_x_continuous(breaks = seq(1850,2020,50)) +
  labs(x = "",y = "") +
  theme_bw()

#########################################################
# Seasonal cycle


s.cycle <- df.all.long.ensemble %>%
  dplyr::filter(var %in% c("dewpoint")) %>%
  group_by(month,var,ensemble.member) %>%
  summarise(value.m = mean(value,na.rm = TRUE),
            .groups = "keep")

s.cycle.reanalysis <- df.all.long %>%
  dplyr::filter(var %in% c("dewpoint")) %>%
  group_by(month,var) %>%
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
  geom_line(data = s.cycle.reanalysis,
            aes(x = month,y = value.m), color = "red", linetype = 2) +
  facet_wrap(~var,scales = "free_y") +
  labs(x = "",y = "") +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J","F","M","A","M","J",
                                "J","A","S","O","N","D")) +
  theme_bw()

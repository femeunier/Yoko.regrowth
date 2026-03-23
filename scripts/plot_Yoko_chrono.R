rm(list = ls())

library(ncdf4)
library(ggplot2)
library(dplyr)
library(reshape2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(tidyr)
library(plotbiomes)
library(minpack.lm)

####################################################################################################

tmax = 150
threshold = TRUE # To compute C stocks (max growth per year)?

# Data

# file <- "/home/femeunier/Desktop/FWO/site.csv"
# data.site <- read.csv(file)
# data.wide <- data.site %>%
#   mutate(agb = agb) %>%
#   pivot_wider(names_from = type,
#               values_from = agb) %>%
#   mutate(t = case_when (t > 100 ~ as.integer(tmax),
#                         TRUE ~ t)) %>%
#   mutate(timing = case_when(t == tmax ~ "Inf",
#                             TRUE ~ as.character(t)))

if (threshold){
  data.wide <- readRDS("./data/Yoko_AGB.tot_dyn.RDS") %>%
    rename(t = age.num,
           mean = agb.threshold.av,
           sd = agb.threshold.sd) %>%
    mutate(timing = as.character(t),
           up = mean+1.96*(sd)/sqrt(3),
           low = mean-1.96*(sd)/sqrt(3)) %>%
    mutate(timing = case_when(t == max(t) ~ "Inf",
                              TRUE ~ timing))
} else {
  data.wide <- readRDS("./data/Yoko_AGB.tot_dyn.RDS") %>%
    rename(t = age.num,
           mean = agb.av,
           sd = agb.sd) %>%
    mutate(timing = as.character(t),
           up = mean+1.96*(sd)/sqrt(3),
           low = mean-1.96*(sd)/sqrt(3)) %>%
    mutate(timing = case_when(t == max(t) ~ "Inf",
                              TRUE ~ timing))
}




m0 <- nlsLM(data = data.wide,
            mean ~ a*(1 - exp(-b*t))**c,
            start=list(a = 300,
                       b = 0.001,
                       c = 1),
            lower = c(0,0,0),
            upper = c(Inf,Inf,Inf),
            control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                  printEval = TRUE, warnOnly = TRUE))


df.fit <- data.frame(time = seq(0,tmax),
                     agb = coef(m0)[1]*(1 - exp(-coef(m0)[2]*seq(0,tmax)))**coef(m0)[3])

# df.fit <- data.frame(time = seq(0,tmax),
#                      agb = approx(data.wide$t,
#                                   data.wide$mean,
#                                   seq(0,tmax))[["y"]])

########################################################################################

system2("rsync",paste("-avz",
                      "hpc:/data/gent/vo/000/gvo00074/felicien/R/outputs/AGB_chronoseq_Yoko.RDS",
                      "./outputs/"))

df_Yoko <- readRDS(file.path("./outputs/","AGB_chronoseq_Yoko.RDS")) %>%
  filter(phen == 2) %>%
  mutate(AGB = AGB.census)

ggplot(data = df_Yoko %>% filter(year >= 1550)) +
  geom_line(aes(x = year, y = 10*AGB,
                group = interaction(as.factor(phen),timing),
                color = as.factor(timing))) +
  geom_point(data = data.wide,
             aes(x = 2018,
                 y = mean,
                 color = as.factor(timing)),alpha = 0.7,
             size = 2) +
  scale_x_continuous(limits = c(2000,2020)) +
  facet_wrap(~as.factor(phen)) +
  theme_bw()

ggplot(data = df_Yoko %>%
         group_by(timing) %>%
         mutate(year0 = year - min(year)) %>%
         filter(year0 <=  tmax)) +
  geom_line(aes(x = year0, y = 10*AGB,color = as.factor(timing))) +
  geom_point(data = data.wide %>% filter(t <= 60),
             aes(x = t,
                 y = mean,
                 color = as.factor(t)),alpha = 0.7,
             size = 2) +
  facet_wrap(~as.factor(phen)) +
  theme_bw()


df.chrono <- bind_rows(list(data.frame(year = 2020,
                                       timing = 0,
                                       AGB = 0,
                                       AGB.census = 0),
                            df_Yoko %>% filter(year == 2018))) %>%
  mutate(timing = case_when(timing == Inf ~ tmax,
                            TRUE ~ timing))



m1 <- nlsLM(data = df.chrono,
            AGB ~ a*(1 - exp(-b*timing))**c,
            start=list(a = 300,
                       b = 0.001,
                       c = 1),
            lower = c(0,0,0),
            upper = c(Inf,Inf,Inf),
            control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                  printEval = TRUE, warnOnly = TRUE))


df.fit.mod <- data.frame(time = seq(0,tmax),
                         agb = 10*coef(m1)[1]*(1 - exp(-coef(m1)[2]*seq(0,tmax)))**coef(m1)[3])

# df.fit.mod <- data.frame(time = seq(0,tmax),
#                          agb = 10*approx(df.chrono$timing,
#                                          df.chrono$AGB,
#                                          seq(0,tmax))[["y"]])


###################################################################################

ggplot(data = df.chrono,
       aes(x = timing, y = AGB*10)) +
  geom_point(color = "red",shape = 2) +
  geom_line(data = df.fit.mod,
            aes(x = time, y = agb),
            color= "red") +
  geom_line(data = df.fit,
            aes(x = time, y = agb),
            color= "grey") +
  geom_point(data = data.wide,
             aes(x = t,
                 y = mean),
             color = "black",alpha = 0.7,
             size = 2) +
  geom_errorbar(data = data.wide,aes(x = t,y = mean,ymin = low,ymax = up),
                width = 1,
                color = "black",alpha = 0.7) +
  theme_bw()

######################################################################################

df.diff <- df.fit %>%
  rename(obs = agb) %>%
  left_join(df.fit.mod %>%
              rename(mod = agb),
            by = "time") %>%
  mutate(diff = mod - obs,
         diff.rel = 100*(mod - obs)/obs)


ggplot(data = df.diff) +
  geom_line(aes(x = time, y = diff)) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()


ggplot(data = df.diff) +
  geom_line(aes(x = time, y = diff.rel)) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw()



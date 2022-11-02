rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(rhdf5)
library(lubridate)
library(wesanderson)
library(minpack.lm)

t.of.dist <- c(LPJ = 931,
               ED2 = 2200,
               ORCHIDEE = 17)

tmax = 1000


#######################################################################################################
# Data

data.wide <- bind_rows(list(data.frame(age.num = 0, agb.av = 0, agb.sd = 0),
                            readRDS("./data/Yoko_AGB.tot_dyn.RDS"))) %>%
  mutate(age.num = case_when(age.num < 100 ~ age.num,
                             TRUE ~ tmax)) %>%
  mutate(agb.av= agb.av/10,
         agb.sd = agb.sd/10) %>%
  rename(mean = agb.av,
         t = age.num,
         sd = agb.sd) %>%
  mutate(low = mean - sd,
         up = mean + sd) %>%
  mutate(t.since = case_when(t > 100 ~ -5,
                             TRUE ~ t))

#######################################################################################################
# LPJ
raw.LPJ <- read.table("./data/LPJ-GUESS/new_runs_14July/output_default/cmass.out",header = TRUE)
cmass <- raw.LPJ %>%
  dplyr::select(-Lon,-Lat) %>%
  pivot_longer(cols = -Year,
               names_to = "pft",
               values_to = "cmass") %>%
  filter(pft == "Total") %>%
  mutate(t.since = Year - t.of.dist["LPJ"])

raw.fit.LPJ <- raw.LPJ %>%
  filter(Year <= (t.of.dist["LPJ"] + 1000)) %>%
  mutate(Year = case_when((Year < t.of.dist["LPJ"]) ~ tmax,
                          TRUE ~ (as.integer(Year) - t.of.dist["LPJ"])))

m.LPJ <- nlsLM(data = raw.fit.LPJ,
               # Total ~ a*(1 - exp(-b*Year)),
               Total ~ a*(1 - exp(-b*Year))**c,
               start=list(a = 300,
                          b = 0.001,
                          c  = 1),
               lower = c(0,0,0),
               upper = c(Inf,Inf,Inf),
               control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                     printEval = TRUE, warnOnly = TRUE))

df.fit.LPJ <- data.frame(time = seq(0,tmax),
                         agb = coef(m.LPJ)[1]*(1 - exp(-coef(m.LPJ)[2]*seq(0,tmax)))**coef(m.LPJ)[3])

# df.fit.LPJ <- data.frame(time = seq(0,tmax),
#                          agb = approx(raw.fit.LPJ %>% arrange(Year) %>% pull(Year),
#                                       raw.fit.LPJ %>% arrange(Year) %>% pull(Total),
#                                       seq(0,tmax))[["y"]])

plot(raw.fit.LPJ$Year,raw.fit.LPJ$Total)
lines(df.fit.LPJ$time,df.fit.LPJ$agb, col = "red")


#######################################################################################################
# ED2

# system2("rsync",paste("-avz",
#                       "hpc:/data/gent/vo/000/gvo00074/felicien/R/df_Yoko_chronosequence.RDS",
#                       "./outputs/"))

df_OP_SA_Yoko.ref <- readRDS(file.path("./data/","df_Yoko_chronosequence.RDS")) %>%
  mutate(t = yr + (month - 1)/12/2) %>%
  rename(agb = AGB) %>%
  mutate(t.since = t - t.of.dist["ED2"]) %>%
  mutate(agb = case_when(t.since <= 0 ~ agb + 5,
                         TRUE ~ agb))


raw.fit.ED2 <- df_OP_SA_Yoko.ref %>%
  filter(yr <= 2300,yr > 1800) %>%
  mutate(yr = case_when((yr <= t.of.dist["ED2"]) ~ tmax,
                        TRUE ~ (as.integer(yr) - t.of.dist["ED2"] - 1)))

m.ED2 <- nlsLM(data = raw.fit.ED2,
               # agb ~ a*(1 - exp(-b*yr)),
               agb ~ a*(1 - exp(-b*yr))**c,
               start=list(a = 300,
                          b = 0.001,
                          c = 1),
               lower = c(0,0,0),
               upper = c(Inf,Inf,Inf),
               control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                     printEval = TRUE, warnOnly = TRUE))

df.fit.ED2 <- data.frame(time = seq(0,tmax),
                         agb = coef(m.ED2)[1]*(1 - exp(-coef(m.ED2)[2]*seq(0,tmax)))**coef(m.ED2)[3])

# df.fit.ED2 <- data.frame(time = seq(0,tmax),
#                          agb = approx(raw.fit.ED2 %>% arrange(yr) %>% pull(yr),
#                                       raw.fit.ED2 %>% arrange(yr) %>% pull(agb),
#                                       seq(0,tmax))[["y"]])

plot(raw.fit.ED2$yr,raw.fit.ED2$agb)
lines(df.fit.ED2$time,df.fit.ED2$agb, col = "red")

#######################################################################################################
# ORCHIDEE

h5file <- file.path('/home/femeunier/Documents/projects/Yoko.regrowth/data/ORCHIDEE','yoko_orc_10y_precut_stomate.nc')

# h5ls(h5file)
# h5readAttributes(h5file,"/time_centered")
# # Description of PFTs
# PFT_NAME__01=SoilBareGlobal
# PFT_NAME__02=BroadLeavedEvergreenTropical
# PFT_NAME__03=BroadLeavedRaingreenTropical
# PFT_NAME__04=NeedleleafEvergreenTemperate
# PFT_NAME__05=BroadLeavedEvergreenTemperate
# PFT_NAME__06=BroadLeavedSummergreenTemperate
# PFT_NAME__07=NeedleleafEvergreenBoreal
# PFT_NAME__08=BroadLeavedSummergreenBoreal
# PFT_NAME__09=LarixSpBoreal
# PFT_NAME__10=C3GrassTemperate
# PFT_NAME__11=C4GrassTemperate
# PFT_NAME__12=C3AgricultureTemperate
# PFT_NAME__13=C4AgricultureTemperate

# LEAF_M, SAP_M_AB, SAP_M_BE, HEART_M_AB, HEART_M_BE, ROOT_M, RESERVE_M, FRUIT_M

# Precut

time <- as.Date(h5read(h5file,"/time_centered")/86400,
                origin = "1991-01-01 00:00:00")

AGB <- h5read(h5file,"/SAP_M_AB")[1,1,,]/1000 +
  h5read(h5file,"/LEAF_M")[1,1,,]/1000 +
  h5read(h5file,"/HEART_M_AB")[1,1,,]/1000 +
  h5read(h5file,"/RESERVE_M")[1,1,,]/1000 +
  h5read(h5file,"/FRUIT_M")[1,1,,]/1000  # kgC/m²

# maxvegetfrac <-  h5read(h5file,"/maxvegetfrac")[1,1,,]
maxvegetfrac <- c( 0.0306,0.7394,0.0517,0.0055,0.0000,0.0000,0.0000,0.0000,0.0000,0.0695,0.0701,0.0271,0.0061)
AGB.tot <- sapply(1:dim(AGB)[2],function(i){
  forest <- c(1:11)
  return(weighted.mean(AGB[forest,i],maxvegetfrac[forest]))})

df_OP_ORCH_precut <- data.frame(times = time, agb = AGB.tot) %>%
  mutate(t = year(time)) %>%
  mutate(yr = t - min(t)) %>%
  group_by(yr) %>%
  summarise(agb = mean(agb),
            .groups = "keep")

# Postcut

h5files <- file.path('/home/femeunier/Documents/projects/Yoko.regrowth/data/ORCHIDEE',
                     c('yoko_orc_postcut_stomate.nc',
                       'yoko_orc_postcut_stomate2.nc',
                       'yoko_orc_postcut_stomate3.nc'))

df_OP_ORCH_postcut <- data.frame() ; cyr <- 0

for (i in seq(1,length(h5files))){

  h5file <- h5files[i]

  Origin <- h5readAttributes(h5file,"/time_centered")[["time_origin"]]
  time <- as.Date(h5read(h5file,"/time_centered")/86400,
                  origin = Origin)
  AGB <- h5read(h5file,"/SAP_M_AB")[1,1,,]/1000 +
    h5read(h5file,"/LEAF_M")[1,1,,]/1000 +
    h5read(h5file,"/HEART_M_AB")[1,1,,]/1000 +
    h5read(h5file,"/RESERVE_M")[1,1,,]/1000 +
    h5read(h5file,"/FRUIT_M")[1,1,,]/1000  # kgC/m²

  # maxvegetfrac <-  h5read(h5file,"/maxvegetfrac")[1,1,,]
  AGB.tot <- sapply(1:dim(AGB)[2],function(i){

    forest <- c(1:11)
    return(weighted.mean(AGB[forest,i],maxvegetfrac[forest]))})

  temp.df <- data.frame(times = time, agb = AGB.tot) %>%
    mutate(t = year(time),
           m = month(time),
           d = day(time))  %>%
    filter(m == 1 & d == 1) %>%
    mutate(yr = t)

  df_OP_ORCH_postcut <- bind_rows(list(df_OP_ORCH_postcut,
                                       temp.df %>% dplyr::select(yr,agb)))

  cyr <- cyr
}


df_OP_ORCH <- bind_rows(list(df_OP_ORCH_precut,
                             df_OP_ORCH_postcut %>%
                               mutate(yr = 1:nrow(df_OP_ORCH_postcut)) %>%
                               ungroup() %>%
                               dplyr::select(yr,agb) %>%
                               mutate(yr = yr + max(df_OP_ORCH_precut$yr))
)) %>%
  mutate(t.since =  yr - t.of.dist["ORCHIDEE"]) %>% ungroup()

raw.fit.ORCHIDEE <- bind_rows(list(df_OP_ORCH_postcut %>% ungroup() %>%
                                     mutate(yr = 1:nrow(df_OP_ORCH_postcut)),
                                   df_OP_ORCH_precut %>% mutate(yr = tmax)))

m.ORCHIDEE <- nlsLM(data = raw.fit.ORCHIDEE,
                    # agb ~ a*(1 - exp(-b*yr)),
                    agb ~ a*(1 - exp(-b*yr))**c,
                    start=list(a = 300,
                               b = 0.001,
                               c = 1),
                    lower = c(0,0,0),
                    upper = c(Inf,Inf,Inf),
                    control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                          printEval = TRUE, warnOnly = TRUE))

df.fit.ORCHIDEE <- data.frame(time = seq(0,tmax),
                              agb = coef(m.ORCHIDEE)[1]*(1 - exp(-coef(m.ORCHIDEE)[2]*seq(0,tmax)))**coef(m.ORCHIDEE)[3])

# df.fit.ORCHIDEE <- data.frame(time = seq(0,tmax),
#                               agb = approx(raw.fit.ORCHIDEE %>% arrange(yr) %>% pull(yr),
#                                            raw.fit.ORCHIDEE %>% arrange(yr) %>% pull(agb),
#                                            seq(0,tmax))[["y"]])

plot(raw.fit.ORCHIDEE$yr,raw.fit.ORCHIDEE$agb)
lines(df.fit.ORCHIDEE$time,df.fit.ORCHIDEE$agb, col = "red")



#######################################################################################################


OP.all <- bind_rows(list(
  cmass %>% dplyr::select(t.since,cmass) %>% rename(AGB = cmass) %>% mutate(model = "LPJ-GUESS"),
  df_OP_SA_Yoko.ref %>% dplyr::select(t.since,agb) %>% rename(AGB = agb) %>% mutate(model = "ED2"),
  df_OP_ORCH %>% dplyr::select(t.since, agb) %>% rename(AGB = agb) %>% mutate(model = "ORCHIDEE")))

########################################################################################################
pal <- wes_palette("Zissou1", 3, type = "continuous")

ggplot() +
  geom_point(data = data.wide,
             aes(x = t.since,
                 y = mean),
             color = "black",
             size = 2) +
  geom_errorbar(data = data.wide,
                aes(x = t.since,y = mean,ymin = low,ymax = up),
                width = 0,
                color = "black") +
  labs(x = "Time since disturbance (yr)", y = "AGB (kgC/m²)", color = "") +
  theme_bw() +
  scale_x_continuous(limits = c(-10,80)) +
  scale_y_continuous(limits = c(0,30)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(text = element_text(size = 22))

ggsave(last_plot(),filename = "./Figures/MI_regrowth_data.png",dpi = 300, width = 20, height = 10,unit = "cm")


ggplot(data = OP.all %>% filter(t.since >= -10,
                                t.since <= 80)) +
  geom_line(aes(x = t.since, y = AGB, color = model),show.legend = FALSE) +
  geom_point(data = data.wide,
             aes(x = t.since,
                 y = mean),
             color = "black",
             size = 2) +
  geom_errorbar(data = data.wide,
                aes(x = t.since,y = mean,ymin = low,ymax = up),
                width = 0,
                color = "black") +
  labs(x = "Time since disturbance (yr)", y = "AGB (kgC/m²)", color = "") +
  theme_bw() +
  scale_x_continuous(limits = c(-10,80)) +
  scale_y_continuous(limits = c(0,30)) +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(text = element_text(size = 22))

ggsave(last_plot(),filename = "./Figures/MI_regrowth_models.png",dpi = 300, width = 20, height = 10,unit = "cm")



MIP <- OP.all %>% filter(t.since >= -10,
                         t.since <= 80) %>%
  mutate(yr = floor(t.since)) %>%
  group_by(yr,model) %>%
  summarise(AGB = min(AGB),
            t.since = min(t.since),
            .groups = "keep") %>%
  filter(t.since %in% seq(-10,80,1)) %>%
  group_by(t.since) %>%
  summarise(AGB.m = mean(AGB),
            AGB.sd = sd(AGB),
            AGB.se = 1.96*sd(AGB)/sqrt(3),
            .groups = "keep")


data.select <- OP.all %>%
  left_join(data.wide %>% dplyr::select(t.since,mean),
            by = "t.since") %>%
  filter(!is.na(mean))


ggplot(data = MIP %>% filter(t.since <= 77),
       aes(x = t.since, y = AGB.m)) +
  geom_line(color = "darkgrey",
            show.legend = FALSE) +
  geom_ribbon(aes(ymin = pmax(0,AGB.m - AGB.se), ymax = AGB.m + AGB.se), alpha = 0.4, fill = "grey",color = NA) +
  geom_point(data = data.wide,
             aes(x = t.since,
                 y = mean),
             color = "black",
             size = 2) +
  geom_point(data = data.select,
             aes(x = t.since,
                 y = AGB,
                 shape = model),
             position = position_jitter(w = 2, h = 0),
             color = "black") +
  geom_errorbar(data = data.wide,
                aes(x = t.since,y = mean,ymin = low,ymax = up),
                width = 0.,
                color = "black") +
  labs(x = "Time since disturbance (yr)", y = "AGB (kgC/m²)", color = "") +
  scale_x_continuous(limits = c(-7,75),expand = c(0,1)) +
  scale_y_continuous(limits = c(0,31),expand = c(0,0)) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(0,1,2)) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(text = element_text(size = 22),
        panel.grid = element_blank(),
        axis.line =  element_line(arrow = grid::arrow(length = unit(0.3, "cm")),
                                  color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))


ggsave(plot = last_plot(),
       filename = "./Figures/Figure_chrono.png",dpi = 300,
       width = 30, height = 15, unit = "cm")


##################################################################################################################
# Fraction of recovery

OP.all.rel <- OP.all %>% filter(t.since >= -10,
                                t.since <= 80) %>%
  group_by(model) %>%
  mutate(AGB.OG = mean(AGB[t.since < 0])) %>%
  mutate(AGB.rel = AGB/AGB.OG)

data.wide.rel <- data.wide %>%
  rename(AGB = mean) %>%
  mutate(AGB.OG = mean(AGB[t.since < 0])) %>%
  mutate(AGB.rel = AGB/AGB.OG,
         low = low/AGB.OG,
         up = up/AGB.OG)


ggplot(data = OP.all.rel) +
  geom_line(aes(x = t.since, y = AGB.rel*100, color = model),show.legend = FALSE) +
  geom_point(data = data.wide.rel,
             aes(x = t.since,
                 y = AGB.rel*100),
             color = "black",
             size = 2) +
  geom_errorbar(data = data.wide.rel,
                aes(x = t.since,y = AGB.rel*100,ymin = low*100,ymax = up*100),
                width = 0,
                color = "black") +
  labs(x = "Time since disturbance (yr)", y = "Fraction of recovery (%)", color = "") +
  theme_bw() +
  scale_x_continuous(limits = c(-10,80)) +
  scale_y_continuous(limits = c(0,1.5*100)) +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme(text = element_text(size = 22),
        legend.position = c(0.8,0.8))

ggsave(last_plot(),filename = "./Figures/MI_regrowth_recovery_models.png",dpi = 300, width = 20, height = 10,unit = "cm")


######################################################################################################################
# Diff

# Data fit
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
#                                   seq(0,tmax)))

plot(data.wide$t,data.wide$mean,ylim = c(0,30))
lines(df.fit$time,df.fit$agb,col = "red")

# Model fit

df.all.fit <- df.fit.ORCHIDEE %>%
  rename(ORCHIDEE = agb) %>%
  left_join(df.fit.ED2 %>% rename(ED2 = agb),
            by = "time") %>%
  left_join(df.fit.LPJ %>% rename(LPJ = agb),
            by = "time") %>%
  pivot_longer(cols = c(ORCHIDEE,ED2,LPJ),
               values_to = "agb",
               names_to = "model")

df.ME <- df.all.fit  %>%
  group_by(time) %>%
  summarise(agb.m = mean(agb,na.rm = TRUE))

df.model.vs.data <- df.fit %>% rename(data = agb) %>%
  left_join(df.ME %>% rename(model = agb.m),
            by = "time") %>%
  mutate(diff = model - data,
         diff.rel = (model - data)/data)

ggplot() +
  geom_line(data = df.model.vs.data,
            aes(x = 1+time, y = data)) +
  geom_line(data = df.model.vs.data,
            aes(x = 1 + time, y = model),linetype = 2) +
  geom_line(data = df.all.fit,
            aes(x = 1 + time, y = agb, color = model)) +
  geom_point(data = data.wide.rel,
             aes(x = 1 + t,
                 y = AGB),
             color = "black",
             size = 2) +
  scale_x_log10() +
  theme_bw()

ggplot(data = df.model.vs.data %>% filter(time > 10)) +
  geom_line(aes(x = time, y = diff.rel)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Time since disturbance (yr)", y = "Model - data \r\n (kgC/m²)", color = "") +
  # scale_x_log10() +
  scale_x_continuous(breaks = seq(0,tmax,250)) +
  theme_bw() +
  theme(text = element_text(size = 22),
        panel.grid = element_blank(),
        axis.line =  element_line(arrow = grid::arrow(length = unit(0.3, "cm")),
                                  color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))


ggsave(last_plot(),
       filename = "./Figures/MI_regrowth_diff.png",
       dpi = 300,
       width = 11, height = 7,unit = "cm")









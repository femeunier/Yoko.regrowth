rm(list = ls())

library(dplyr)
library(ggplot2)
library(minpack.lm)
library(tidyr)
library(stringr)

data.Yoko <- read.csv("./data/yoko_inventory.csv",stringsAsFactors = FALSE)
data.Yoko.m <- data.Yoko %>%
  mutate(dbh = rowMeans(dplyr::select(data.Yoko, starts_with("DBH")),na.rm = TRUE)) %>%
  mutate(plot.num = as.numeric(substr(plots,nchar(stage)+1,nchar(stage)+1)))


sp.mapping <- read.csv("./data/sp_mapping.csv",stringsAsFactors = FALSE)

########################################################################################
# Height allom and AGB

href <- 61.7;b1Ht <- 0.0352;b2Ht <- 0.694
m0 <- nlsLM(data = data.Yoko.m,
            height ~ href*(1 -exp(-b1Ht*((dbh)**b2Ht))),
            start=list(href=href, b1Ht=b1Ht, b2Ht = b2Ht), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                                 printEval = TRUE, warnOnly = TRUE))
y <- data.Yoko.m$height
f <- predict(m0)
1-sum((y[!is.na(y)]-f)^2,na.rm = TRUE)/(length(y[!is.na(y)])*var(y[!is.na(y)]))

age.OG = 250

data.Yoko.h <- data.Yoko.m %>% mutate(h = case_when(!is.na(height) ~ height,
                                                    TRUE ~ coef(m0)[1]*(1 -exp(-coef(m0)[2]*((dbh)**coef(m0)[3]))))) %>%
  mutate(agb = 0.0673*(meanWD*(dbh**2)*h)**0.976) %>%
  mutate(age.num = case_when(Age == ">100" ~ as.character(age.OG),
                             TRUE ~ (Age))) %>%
  mutate(age.num = as.numeric(age.num)) %>%
  mutate(PFT = case_when(is.na(FunctionalG) ~ "Other",
                         TRUE ~ FunctionalG))

########################################################################################
# Aggregate

data.plot <-
  data.Yoko.h %>%
  group_by(age.num,plots,plot.num,PFT) %>%
  summarise(agb.m = sum(agb,na.rm = TRUE)/(40*40)/2*10, # MgC/ha
            .groups = "keep")

data.plot.pft.m <- data.plot %>%
  ungroup() %>%
  dplyr::select(age.num,PFT,plot.num,agb.m) %>%
  complete(age.num = c(5,12,20,60,age.OG),
           PFT = c("NPLD","PIO","SB","Other"),
           plot.num = c(1,2,3),
           fill = list(agb.m = 0)) %>%
  group_by(age.num,PFT) %>%
  summarise(agb.av = mean(agb.m),
            agb.sd = sd(agb.m),
            .groups = "keep") %>%
  ungroup()

data.plot.m <- data.plot %>%
  group_by(age.num,plots) %>%
  summarise(agb.m = sum(agb.m),
            .groups = "keep") %>%
  group_by(age.num) %>%
  summarise(agb.av = mean(agb.m),
            agb.sd = sd(agb.m),
            .groups = "keep") %>%
  ungroup()

############################################################################################
# AGB + AGB component dynamics

ggplot() +
  geom_bar(data = data.plot.pft.m,
           aes(x = age.num,y = agb.av,fill = PFT), stat = "identity") +
  geom_point(data = data.plot.m,
             aes(x = age.num, y = agb.av), size = 1) +
  geom_errorbar(data = data.plot.m,
                aes(x = age.num, y = agb.av,
                    ymin = agb.av - agb.sd, ymax = agb.av + agb.sd),width = 0) +
  theme_bw()

saveRDS(data.plot.pft.m,"./data/Yoko_AGB_dyn.RDS")
saveRDS(data.plot.m,"./data/Yoko_AGB.tot_dyn.RDS")

#############################################################################################
# Tree size distributions

data.Yoko.dbh <- data.Yoko.h %>%
  ungroup() %>%
  dplyr::select(age.num,PFT,plot.num,dbh) %>%
  mutate(dbh.cat = case_when(dbh < 10 ~ 1,
                             dbh < 20 ~ 2,
                             dbh < 30 ~ 3,
                             dbh < 40 ~ 4,
                             dbh < 50 ~ 5,
                             dbh < 60 ~ 6,
                             dbh < 70 ~ 7,
                             dbh < 80 ~ 8,
                             dbh < 90 ~ 9,
                             dbh < 100 ~ 10,
                             TRUE ~ 11),
         BA = dbh*dbh*pi/4)

data.Yoko.pft.dbh <- data.Yoko.dbh %>%
  group_by(age.num,plot.num,dbh.cat,PFT) %>%
  summarise(N = length(dbh[!is.na(dbh)])/(40*40)*10000, # Number/ha
            BA = sum(BA, na.rm = TRUE)/(40*40),         # cm²/m²
            .groups = "keep") %>%
  group_by(age.num,dbh.cat,PFT) %>%
  summarise(N.m = mean(N,na.rm = TRUE),
            BA.m = mean(BA,na.rm = TRUE),
            .groups = "keep")

data.Yoko.pft.dbh.tot <- data.Yoko.dbh %>%
  group_by(age.num,plot.num,dbh.cat) %>%
  summarise(N = length(dbh[!is.na(dbh)])/(40*40)*10000, # Number/ha
            BA = sum(BA, na.rm = TRUE)/(40*40),         # cm²/m²
            .groups = "keep") %>%
  ungroup() %>%
  complete(age.num = c(5,12,20,60,age.OG),
           plot.num = c(1,2,3),
           dbh.cat = seq(2,11),
           fill = list(N = 0,
                       BA = 0)) %>%
  group_by(age.num,dbh.cat) %>%
  summarise(N.m = mean(N,na.rm = TRUE),
            N.sd = sd(N,na.rm = TRUE),
            BA.m = mean(BA,na.rm = TRUE),
            .groups = "keep")

ggplot() +
  geom_bar(data = data.Yoko.pft.dbh %>% filter(dbh.cat > 1),
           aes(x = as.factor(dbh.cat),y = N.m,fill = PFT), stat = "identity") +
  facet_wrap(~as.factor(age.num), nrow = 1) +
  theme_bw() +
  labs(x = "",y = "Density (#/m²)")


ggplot() +
  geom_bar(data = data.Yoko.pft.dbh %>% filter(dbh.cat > 1),
           aes(x = as.factor(dbh.cat),y = BA.m,fill = PFT), stat = "identity") +
  facet_wrap(~as.factor(age.num), nrow = 1) +
  theme_bw() +
  labs(x = "",y = "BA (cm²/m²)")

saveRDS(data.Yoko.pft.dbh,"./data/Yoko_invent_dbh.RDS")
saveRDS(data.Yoko.pft.dbh.tot,"./data/Yoko_invent_dbh_tot.RDS")

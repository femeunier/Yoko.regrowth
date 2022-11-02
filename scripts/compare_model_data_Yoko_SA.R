rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(minpack.lm)

system2("scp",paste("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/df_OP_SA_Yoko.RDS","./outputs/"))
system2("scp",paste("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/df_param_SA_Yoko.RDS","./outputs/"))
system2("scp",paste("hpc:/kyukon/data/gent/vo/000/gvo00074/felicien/R/df_nplant_Yoko.RDS","./outputs/"))

############################################################################################################
# Data
tmax = 1998      # SA
tmax.OG = 250    # data


data.wide <- readRDS("./data/Yoko_AGB.tot_dyn.RDS") %>%
  mutate(agb.av= agb.av/10,
         agb.sd = agb.sd/10) %>%
  rename(mean = agb.av,
         t = age.num,
         sd = agb.sd) %>%
  mutate(low = mean - sd,
         up = mean + sd)


m0 <- nlsLM(data = data.wide,
            mean ~ a*(1 - exp(-b*t)),
            start=list(a = 300,
                       b = 0.001),
            lower = c(0,0),
            upper = c(Inf,Inf),
            control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                  printEval = TRUE, warnOnly = TRUE))


df.fit <- data.frame(time = seq(0,tmax.OG),
                     agb = coef(m0)[1]*(1 - exp(-coef(m0)[2]*seq(0,tmax.OG))))

############################################################################################################
# Model
df_OP_SA_Yoko <- readRDS("./outputs/df_OP_SA_Yoko.RDS")

df_OP_SA_Yoko.ref <- df_OP_SA_Yoko %>%   filter(param == "reference")
df_OP_SA_Yoko <- df_OP_SA_Yoko %>%   filter(param != "reference")

params <- unique(df_OP_SA_Yoko$param)

for (i in seq(1:length(params))){
  df_OP_SA_Yoko <- bind_rows(list(df_OP_SA_Yoko,
                                  df_OP_SA_Yoko.ref %>% mutate(param = params[i])))
}

df_OP_SA_Yoko <- df_OP_SA_Yoko %>%
  mutate(agb.early.frac = agb.early/agb) %>%
  pivot_longer(cols = -c(t,quantile,param),
               values_to = "val",
               names_to = "var")


df_OP_SA_Yoko.end <- df_OP_SA_Yoko %>%
  ungroup() %>%
  filter(t == tmax) %>%
  group_by(param,var) %>%
  mutate(val.rel = val/val[quantile == 0.5],
            .groups = "keep")

t.max.sim <- max(df_OP_SA_Yoko$t) - min(df_OP_SA_Yoko$t)

df.area <- df_OP_SA_Yoko %>%
  filter(quantile %in% c(0.025,0.5,0.975)) %>%
  mutate(quantile.str = case_when(quantile == 0.025 ~ "low",
                                  quantile == 0.5 ~ "median",
                                  quantile == 0.975 ~ "high")) %>%
  dplyr::select(-quantile) %>%
  pivot_wider(names_from = quantile.str,
              values_from = val) %>%
  rowwise() %>%
  mutate(temp.high = pmax(low,high),
         temp.low = pmin(low,high)) %>%
  mutate(low = temp.low,
         high = temp.high)


# Dynamics
cvar = "agb"
ggplot() +
  geom_ribbon(data = df.area %>% filter(var == cvar),
            aes(x = t - t[1],
                ymin = low,ymax = high, fill = param), color = NA, alpha = 0.5) +
  geom_line(data = df_OP_SA_Yoko %>% filter(quantile == 0.5,
                                            var == cvar),
            aes(x = t - t[1],
                y = val,
                linetype = as.factor(quantile)),
            color = "black",linetype = 1) +
  geom_point(data = data.wide,
             aes(x = t,
                 y = mean),
             color = "black",alpha = 0.7,
             size = 2) +
  geom_errorbar(data = data.wide,
                aes(x = t,y = mean,ymin = low,ymax = up),
                width = 1,
                color = "black",alpha = 0.7) +
  geom_line(data = df.fit,
            aes(x = time, y = agb),
            color = "black",linetype = 2) +
  facet_wrap(~ param) +
  scale_x_continuous(limits = c(0,t.max.sim)) +
  scale_linetype_manual(values = c(2,1,3)) +
  theme_bw()


# Dynamics, no data

cvar = "LAI"
ggplot() +
  geom_ribbon(data = df.area %>% filter(var == cvar),
              aes(x = t - t[1], ymin = low,
                  ymax = high, fill = param), color = NA, alpha = 0.5) +
  geom_line(data = df_OP_SA_Yoko %>% filter(quantile == 0.5,
                                            var == cvar),
            aes(x = t - t[1],
                y = val,
                linetype = as.factor(quantile)),
            color = "black",linetype = 1) +
  facet_wrap(~ param) +
  scale_x_continuous(limits = c(0,t.max.sim)) +
  scale_linetype_manual(values = c(2,1,3)) +
  theme_bw()


# SA
ggplot(data = df_OP_SA_Yoko.end,
       aes(x = quantile,
           y = val.rel,
           color = param)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 1,
             color = "black", linetype = 2) +
  facet_wrap(~ var, nrow = 1) +
  theme_bw()

############################################################################################################
# Density plot
df_nplant_Yoko <- readRDS("./outputs/df_nplant_Yoko.RDS") %>% mutate(t = year - 1700)

df_nplant_Yoko.ref <- df_nplant_Yoko %>% filter(quantile == 0.5)
df_nplant_Yoko <- df_nplant_Yoko %>% filter(quantile != 0.5)

for (i in seq(1:length(params))){
  df_nplant_Yoko <- bind_rows(list(df_nplant_Yoko,
                                   df_nplant_Yoko.ref %>% mutate(param = params[i])))
}

data.nplant <- readRDS("./data/Yoko_invent_dbh.RDS") %>% rename(t = age.num) %>%
  mutate(pft = case_when(PFT == "PIO" ~ 2,
                         PFT == "NPLD" ~ 3,
                         PFT == "SB" ~ 4,
                         PFT == "Other" ~ 3))

data.nplant.tot <-  readRDS("./data/Yoko_invent_dbh_tot.RDS") %>% rename(t = age.num)

ggplot() +
  geom_bar(data = df_nplant_Yoko.ref,
           aes(x = as.factor(dbh.cat - 1), y = nplant, fill = as.factor(pft)),
           alpha = 1,
           stat = "identity") +
  geom_point(data = data.nplant.tot %>% filter(dbh.cat > 1),
             aes(x = as.factor(dbh.cat - 1),
                 y = N.m), color = "black") +
  geom_errorbar(data = data.nplant.tot %>% filter(dbh.cat > 1),
                aes(x = as.factor(dbh.cat - 1),
                    y = N.m, ymin = N.m - N.sd, ymax = N.m + N.sd), color = "black",width = 0) +
  facet_grid(~ as.factor(t)) +
  theme_bw()

ggplot(data = df_nplant_Yoko) +
  geom_bar(aes(x = interaction(as.factor(quantile),as.factor(dbh.cat)), y = nplant, fill = as.factor(pft)), stat = "identity") +
  facet_grid(as.factor(t) ~ param) +
  theme_bw()


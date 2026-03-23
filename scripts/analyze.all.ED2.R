rm(list = ls())

library(ggplot2)
library(dplyr)
library(ncdf4)
library(ncdf4.helpers)
library(TrENDY.analyses)
library(lubridate)
library(dismo)
library(textrecipes)
library(polite)
library(dplyr)
library(tidyr)


raw <-  readRDS("./outputs/Compiled.model.outputs.R") %>%
  filter(!is.na(obs))
model.predictions <- raw %>%
  pivot_longer(cols = c(pred,obs),
               names_to = "source",
               values_to = "value") %>%
  filter((source == "pred" & model != "MEM") |
           (source == "obs" & model == "MEM"))


raw %>%
  filter(is.na(obs))

MEM <- model.predictions %>%
  dplyr::select(-c(lon,lat)) %>%
  ungroup() %>%
  distinct() %>%
  mutate(age = case_when(age < 0 ~ 300,
                         TRUE ~ age)) %>%
  filter(age > 0) %>%
  filter((type == "Primary") | (age < 300 & type == "Secondary"))


ggplot(data = MEM %>%
         ungroup() %>%
         filter(region == region[1]),
      aes(x = age,
          y = value,
          color = model)) +
  geom_line() +
  theme_bw()

MEM %>%
  ungroup() %>%
  filter(region == region[1],
         model == "MEM")

df.t <- MEM %>%
  filter(type == "Secondary") %>%
  group_by(region,model,site,country) %>%
  summarise(AGB.20 = approx(x = age,
                            y = value,
                   xout = 20)[["y"]],
            .groups = "keep")

df.eq <- MEM %>%
  filter(type == "Primary")

climate <- readRDS("./outputs/Climate.sites.RDS") %>%
  mutate(region =  gsub(" ", "_", chartr("éèêë", "eeee", tolower(sub(" ","_",gsub(" +", " ", region))))))

sites <- unique(climate$region)
all.climate <- data.frame()
for (csite in sites){

  cdf <- climate %>%
    filter(region == csite) %>%
    mutate(pre = pre*4*30.25)

  all.climate <- bind_rows(all.climate,
                           data.frame(region = csite,
                                      t(biovars(cdf$pre,
                                         (cdf$tmin - 273.15),
                                         (cdf$tmax - 273.15))[c(1:19)])))
}

climate.sum <- climate %>%
  group_by(region) %>%
  summarise(MAT = mean(tmp),
            MAP = sum(pre)*4*30.25,
            dswrf = mean(dswrf),
            spfh = mean(spfh),
            .groups = "keep") %>%
  left_join(all.climate,
            by = "region")


df.t.climate <- df.t %>%
  mutate(region = tolower(region)) %>%
  left_join(climate.sum,
            by = "region")

df.t.climate.long <- df.t.climate %>%
  pivot_longer(cols = c("MAT","MAP","dswrf","spfh",
                        starts_with("X")),
               names_to = "climate_var",
               values_to = "value")
DF2plot <- df.t.climate.long %>%
  filter(climate_var %in% c("X3","X12","X15"))  %>%
  mutate(climate_var = case_when(climate_var == "X3" ~ "Isothermality",
                                 climate_var == "X4" ~ "Temp. seasonality",
                                 climate_var == "X12" ~ "MAP",
                                 climate_var == "X15" ~ "Precip. seasonality"))


DF2plot %>%
  filter(climate_var == "Isothermality",
         value < 50)

ggplot(data = DF2plot %>%
         filter(climate_var %in% c("MAP","Precip. seasonality")) %>%
         na.omit(),
       aes(x = value, y = AGB.20,
           color = model, fill = model)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ poly(x,1),
              se = TRUE) +
  facet_wrap(~climate_var,
             scales = "free",
             nrow = 1) +
  scale_color_manual(values = c("#e7298a","#d95f02","#1b9e77","#7570b3",
                                "black")) +
  scale_fill_manual(values = c("#e7298a","#d95f02","#1b9e77","#7570b3",
                               "black")) +
  facet_wrap(~climate_var,
             scales = "free_x",
             nrow = 1) +
  theme_bw() +
  labs(x = "",y = "") +
  guides(fill = "none", color = "none") +
  theme(text = element_text(size = 20),
        legend.position = c(0.7,0.8),strip.background = element_blank(),
        strip.text = element_blank())

DF2plot %>%
  filter(climate_var %in% c("MAP","Precip. seasonality")) %>%
  group_by(climate_var,model) %>%
  summarise(R2 = summary(lm(AGB.20 ~ value))[["r.squared"]],
            .groups = "keep")

df.t.climate.long %>%
  group_by(model,climate_var) %>%
  summarise(rsq = summary(lm(AGB.20 ~ value))[["r.squared"]],
            .groups = "keep") %>%
  ungroup() %>%
  arrange(desc(rsq))

df.eq.climate <- df.eq %>%
  mutate(region = tolower(region)) %>%
  left_join(climate.sum,
            by = "region")

df.eq.climate.long <- df.eq.climate %>%
  pivot_longer(cols = c("MAT","MAP","dswrf","spfh",
                      starts_with("X")),
             names_to = "climate_var",
             values_to = "value.clim")

df.eq.climate.long %>%
  filter(climate_var == "MAP",
         model == "MEM") %>%
  pull(value.clim) %>%
  hist()

df.eq.climate.long %>%
  group_by(model,climate_var) %>%
  summarise(rsq = summary(lm(value ~ value.clim))[["r.squared"]],
            .groups = "keep") %>%
  group_by(climate_var) %>%
    summarise(R2m = mean(rsq)) %>%
  arrange(desc(R2m))


ggplot(data = df.eq.climate.long %>%
         filter(climate_var %in% c("MAP","X17"))  %>%
         mutate(climate_var = case_when(climate_var == "X10" ~ "Mean temp. of warmest quarter",
                                        climate_var == "X14" ~ "Precip. of driest month",
                                        climate_var == "X17" ~ "Precip. of driest quarter",
                                        climate_var == "X12" ~ "MAP",
                                        TRUE ~ climate_var)),
       aes(x = value.clim, y = value,
           color = model, fill = model)) +
  geom_point() +
  stat_smooth(method = "lm",formula = y ~ poly(x,1),
              se = TRUE) +
  scale_color_manual(values = c("#e7298a","#d95f02","#1b9e77","#7570b3",
                                "black")) +
  scale_fill_manual(values = c("#e7298a","#d95f02","#1b9e77","#7570b3",
                               "black")) +
  facet_wrap(~climate_var,
             scales = "free_x",
             nrow = 1) +
  theme_bw() +
  labs(x = "",y = "") +
  guides(fill = "none", color = "none") +
  theme(text = element_text(size = 20),
        legend.position = c(0.7,0.8),strip.background = element_blank(),
        strip.text = element_blank())

df.eq.climate.long %>%
  filter(climate_var %in% c("MAP","X17"))  %>%
  mutate(climate_var = case_when(climate_var == "X10" ~ "Mean temp. of warmest quarter",
                                 climate_var == "X14" ~ "Precip. of driest month",
                                 climate_var == "X17" ~ "Precip. of driest quarter",
                                 climate_var == "X12" ~ "MAP",
                                 TRUE ~ climate_var)) %>%
  filter(climate_var %in% c("MAP","Precip. of driest quarter")) %>%
  group_by(climate_var,model) %>%
  summarise(R2 = summary(lm(value ~ value.clim))[["r.squared"]],
            .groups = "keep")




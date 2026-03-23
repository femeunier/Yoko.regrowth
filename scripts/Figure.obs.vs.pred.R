rm(list = ls())

library(cowplot)
library(minpack.lm)
library(ggplot2)
library(dplyr)
library(tidyr)

df <- readRDS("./outputs/Compiled.model.outputs.R") %>%
  filter(var == "cAGB")
primary <- df %>%
  filter(type == "Primary") %>%
  pivot_longer(cols = c(obs,pred),
               names_to = "source",
               values_to = "value")
secondary <- df %>%
  filter(type != "Primary") %>%
  filter(age > 0)

# MEM
fit <- nlsLM(obs ~ a * ((1 - exp(-b * age))^ c),
           data = secondary %>%
             na.omit(),
           start = list(a = 10, b = 0.02, c = 1),
           lower = c(a = 0, b = 0, c = 0),   # enforce positivity
           upper = c(a = 20, b = 1, c = 2))

all.ages <- seq(0,70,length.out = 1000)
fit.pred <- data.frame(age = all.ages,
                       AGB = predict(fit,
                                     newdata = data.frame(age = seq(0,60,length.out = 1000))))

fit.model <- nlsLM(pred ~ a * ((1 - exp(-b * age))^ c),
             data = secondary %>%
               na.omit(),
             start = list(a = 10, b = 0.02, c = 1),
             lower = c(a = 0, b = 0, c = 0),   # enforce positivity
             upper = c(a = 20, b = 1, c = 2))

fit.pred.model <- data.frame(age = all.ages,
                             AGB = predict(fit.model,
                                     newdata = data.frame(age = seq(0,60,length.out = 1000))))


df.Viola <- data.frame(age = all.ages) %>%
  mutate(AGB = 113/10*((1- exp(-0.02*age))**0.706))


a <- ggplot() +
  geom_point(data = secondary,
             aes(x = age,
                 y = pred), size = 0.1,
             color = "darkblue", alpha = 0.7) +
  geom_point(data = secondary,
             aes(x = age,
                 y = obs), size = 0.1,
             color = "grey17", alpha = 0.7) +
  # geom_density_2d_filled(,
  #                        aes(
  #                            fill = source),
  #                        alpha = 0.7) +
  # stat_density_2d(data = secondary %>%
  #                   pivot_longer(cols = c(pred,obs),
  #                                names_to = "source",
  #                                values_to = "value"),
  #                 aes(x = age,
  #                     y = value,
  #                     fill = as.factor(source), alpha = ..level..),
  #                 geom = "polygon", contour = TRUE,
  #                 color = NA, bins = 10) +
  geom_line(data = df.Viola,
            aes(x = age, y = AGB), linetype = 2) +
  geom_line(data = fit.pred,
            aes(x = age, y = AGB), linetype = 1) +
  geom_line(data = fit.pred.model,
            aes(x = age, y = AGB), linetype = 1, color = "darkblue") +
  scale_fill_manual(values = c("grey17","darkblue")) +
  scale_y_continuous(limits = c(0,30)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(limits = c(0,70)) +
  labs(x = "", y = "") +
  guides(fill = "none", alpha = "none")
a
b <- ggplot() +
  geom_boxplot(data = primary,
               aes(x = 0,
                   y = value,
                   fill = source),
               width = 10, alpha = 0.6) +
  scale_fill_manual(values = c("grey17","darkblue")) +
  scale_y_continuous(limits = c(0,30),
                     breaks = c(0,10,20,30),
                     labels = c("","","","")) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks = c(0),
                 labels = c("")) +
  labs(x = "", y = "") +
  guides(fill = "none")
b

primary %>%
  group_by(source) %>%
  summarise(median(value,na.rm = TRUE))


summary(aov(data = primary,
            formula = value ~ source))


plot_grid(a,b,align = "hv",rel_widths = c(3,1))

df.all <- bind_rows(df.Viola %>% mutate(source = "Viola"),
                    fit.pred %>% mutate(source = "obs"),
                    fit.pred.model %>% mutate(source = "pred")) %>%
  filter(age <= 20) %>%
  group_by(source) %>%
  arrange(age) %>%
  mutate(growth.rate = c(NA,diff(AGB)/diff(age)))


df.all %>%
  filter(source != "Viola") %>%
  group_by(source) %>%
  summarise(median(growth.rate,na.rm = TRUE))

summary(aov(data = df.all %>%
      filter(source != "Viola"),
    formula = growth.rate ~ source))

df.all %>%
  filter(source != "Viola")

ggplot(data = df.all %>%
         filter(source != "Viola")) +
  geom_boxplot(aes(x = source, y = growth.rate,
                   fill = source),
               alpha = 0.6) +
  theme_bw() +
  scale_fill_manual(values = c("grey17","darkblue")) +
  theme(text = element_text(size = 20)) +
  scale_x_discrete(labels = c('','')) +
  scale_y_log10() +
  labs(x = "", y = "") +
  guides(fill = "none")


df %>%
  filter(region == "LIB_Sapo") %>%
  dplyr::select(age,obs)

ggplot(data = df %>%
         # filter(type != "Primary") %>%
         filter(region != "IC_Yaya") %>%
         filter(!(region == "LIB_Sapo" & (obs>10)))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(x = obs,
                 y = pred, color = model),
             size = 0.5, alpha = 0.8) +
  # stat_density_2d(aes(x = obs,
  #                     y = pred,
  #                     fill = as.factor(model), alpha = ..level..),
  #                 geom = "polygon", contour = TRUE,
  #                 color = NA, bins = 4) +
  stat_smooth(aes(x = obs,
                  y = pred,
                  color = model),
              method = "lm", se = FALSE) +
  facet_wrap(~type) +
  labs(x = "",
       y = "",
       color = "", shape = "") +
  scale_color_manual(values = c("#e7298a","#7570b3","#1b9e77","#66a61e","#d95f02",
                                "black")) +
  scale_fill_manual(values = c("#e7298a","#7570b3","#1b9e77","#66a61e","#d95f02",
                                "black")) +
  guides(color = "none", alpha = "none", fill = "none") +
  scale_x_continuous(limits = c(0,15)) +
  scale_y_continuous(limits = c(0,15)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = c(0.7,0.8),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  coord_equal()


df %>%
  # filter(type != "Primary") %>%
  filter(region != "IC_Yaya") %>%
  filter(!(region == "LIB_Sapo" & (obs>10))) %>%
  group_by(type,model) %>%
  summarise(R2 = summary(lm(formula = pred ~ obs))[["r.squared"]])

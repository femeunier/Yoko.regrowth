rm(list = ls())

library(ggplot2)
library(dplyr)
library(stringr)
library(stringi)
library(ncdf4)
library(ncdf4.helpers)
library(lubridate)
library(tidyr)
library(minpack.lm)


main.dir <- "/Users/felicien/Documents/projects/Yoko.regrowth/outputs/MIP"
models <- dir(main.dir)

all.models <- data.frame()

time.names <- c("time","time_counter")
vars2read <- "agb"

recompile <- FALSE

for (cmodel in models){

  print(cmodel)

  cdir <- file.path(main.dir,cmodel)

  model.output.file <- file.path(cdir,
                                 paste0(cmodel,".outputs.RDS"))

  if (!recompile & file.exists(model.output.file)){
    df.model <- readRDS(model.output.file)
    all.models <- bind_rows(all.models,
                            df.model)
    next()
  }


  files <- list.files(cdir,pattern = "*.nc$",recursive = TRUE,full.names = TRUE)
  files <- files[!grepl("pft|PFT",basename(files))]

  files <- files[!grepl("stomate",basename(files))]

  files.no.ext <- strsplit(tools::file_path_sans_ext(basename(files)),"\\_")

  file.model <- file.country <- file.site <-
    file.type <- file.age <- file.var <-
    c()

  for (ifile in seq(1,length(files.no.ext))){
    clist <- files.no.ext[[ifile]]

    num.max <- ifelse(tolower(clist[length(clist)]) %in% c("ecosystem","pft"),
                      7,6)

    if (length(clist) == num.max){
      file.model[ifile] <- clist[1]
      file.country[ifile] <- clist[2]
      file.site[ifile] <- clist[3]
      file.type[ifile] <- clist[4]
      file.age[ifile] <- clist[5]
      file.var[ifile] <- clist[6]
    } else {

      Ntot <- length(clist)

      Delta <- ifelse(tolower(clist[length(clist)]) %in% c("ecosystem","pft"),
                      1,0)

      if (cmodel %in% c("ORCHIDEE","ORCHIDEEv4.2")){

        file.model[ifile] <- clist[1]
        file.country[ifile] <- strsplit(clist[2],"\\.")[[1]][1]
        file.site[ifile] <- strsplit(clist[2],"\\.")[[1]][2]
        file.type[ifile] <- clist[3]
        file.age[ifile] <- clist[4]
        file.var[ifile] <- clist[5]
      } else {
        file.model[ifile] <- clist[1]
        file.country[ifile] <- clist[2]
        file.site[ifile] <- paste0(clist[3:(Ntot-3-Delta)],collapse = " ")
        file.type[ifile] <- clist[Ntot-2-Delta]
        file.age[ifile] <- clist[Ntot-1-Delta]
        file.var[ifile] <- clist[Ntot-Delta]
      }

    }
  }

  df.files <-
    data.frame(model = file.model,
               country = file.country,
               site = file.site,
               type = file.type,
               age = file.age,
               var = file.var,
               file = files)

  df.files.selected <- df.files %>%
    filter(grepl(tolower(paste0(c(vars2read,"allvariables"),collapse = "|")),
                 tolower(var)))

  df.model <- data.frame()
  for (ifile in seq(1,nrow(df.files.selected))){

    print(ifile/nrow(df.files.selected))

    cfile.model <- df.files.selected$model[ifile];  cfile.country <- df.files.selected$country[ifile];  cfile.site <- df.files.selected$site[ifile]
    cfile.type <- df.files.selected$type[ifile];   cfile.age <- as.numeric(df.files.selected$age[ifile]);   cfile.var <- df.files.selected$var[ifile]

    cfile <- df.files.selected$file[ifile]

    if (!file.exists(cfile)){
      next
    }

    nc <- nc_open(cfile)

    if (cfile.var == "allvariables"){

      if (cmodel == "ORCHIDEEv4.2"){
        cAGB <- (ncvar_get(nc,"LEAF_M_c") +
                   ncvar_get(nc,"SAP_M_AB_c") +
                   ncvar_get(nc,"HEART_M_AB_c") +
                   ncvar_get(nc,"FRUIT_M_c"))/1000

        Frac <- ncvar_get(nc,"VEGET_MAX")

        cdata <- apply(cAGB*Frac,2,sum)
      } else {
        cAGB <- (ncvar_get(nc,"LEAF_M") +
                   ncvar_get(nc,"SAP_M_AB") +
                   ncvar_get(nc,"HEART_M_AB") +
                   ncvar_get(nc,"FRUIT_M"))/1000

        Frac <- ncvar_get(nc,"VEGET_COV_MAX")

        cdata <- apply(cAGB*Frac,2,sum)
      }



    } else {
      cdata <- tryCatch(ncvar_get(nc,cfile.var),
                        error = function(e) NULL)

      if (is.null(cdata) & cmodel == "FATES"){
        cdata <- tryCatch(ncvar_get(nc,"cVegAbove"),
                          error = function(e) NA_real_)
      }
    }

    if (length(cdata) == 0){
      cdata <- NA_real_
    }

    clat <- tryCatch(ncvar_get(nc,"lat"),
                     error = function(e) NA_real_)
    clon <- tryCatch(ncvar_get(nc,"lon"),
                     error = function(e) NA_real_)

    all.times <- NULL ; i = 1
    while(is.null(all.times) & i <= length(time.names)){
      all.times <- tryCatch(suppressMessages(ncvar_get(nc,time.names[i])),
                       error = function(e) NULL)
      i = i +1
    }


    if (cmodel == "ED2.2"){

      start_date <- as.Date("1700-01-01")
      months_after <-all.times

      # Convert months after into dates
      ctimes <- start_date %m+% months(months_after)
    } else if (cmodel == "LPJ-GUESS") {

      start_date <- as.Date("1900-01-01")

      ctimes <- start_date + all.times
    } else if (cmodel == "ORCHIDEE"){
      start_date <- as.Date("1901-01-01")

      days_since_origin <- all.times/86400

      year <- round(year(start_date) + (days_since_origin / 365))
      days_remaining <- days_since_origin %% 365
      months <- month(as.Date(days_remaining,paste0(year(start_date),"/01/01")))
      days <- day(as.Date(days_remaining,paste0(year(start_date),"/01/01")))
      ctimes <- as.Date(paste0(year,"/",months,"/",days))

    } else if (cmodel == "ORCHIDEEv4.2"){

      start_date <-  strsplit(nc$dim$time_counter$units," ")[[1]][3]

      days_since_origin <- all.times/86400

      year <- round(year(start_date) + (days_since_origin / 365))
      days_remaining <- days_since_origin %% 365
      months <- month(as.Date(days_remaining,"1901/01/01"))
      days <- day(as.Date(days_remaining,"1901/01/01"))
      ctimes <- as.Date(paste0(year,"/",months,"/",days))

    } else if (cmodel == "FATES"){

      ctimes <- as.Date(all.times,origin = "1901-01-01 00:00:00")
    }

    if (length(ctimes) == 0){
      ctimes <- NA
    }


    df.model <- dplyr::bind_rows(df.model,
                          data.frame(model = cfile.model,
                                     country = cfile.country,
                                     site = cfile.site,
                                     lat = clat,
                                     lon = clon,
                                     type = cfile.type,
                                     age = cfile.age,
                                     var = cfile.var,
                                     time = ctimes,
                                     value = as.vector(cdata)))

    nc_close(nc)
  }

  all.models <- bind_rows(all.models,
                          df.model)

  if (recompile){
    saveRDS(df.model,
            model.output.file)
  }

}

all.models.corrected <- all.models
all.models.corrected$site <- str_replace_all(str_trim(all.models$site), regex("[-\\s\\u00A0]+"), "_")

# Correct variable names
all.models.corrected <- all.models.corrected %>%
  mutate(site = case_when(site == "Cusseque" ~ "Cusseque_Chutembo_Bie",
                          site == "Mpem" ~ "Mpem_et_Djim_NP",
                          site == "Ht" ~ "Ht_sassandra",
                          TRUE ~ site)) %>%
  mutate(var = case_when(var == "agb" ~ "cAGB",
                         var == "allvariables" ~ "cAGB",
                         TRUE ~ var)) %>%
  mutate(value = case_when(var %in% c("cAGB") & model == "LPJGUESSv4.1.1" ~ value/2,
                           TRUE ~ value)) %>%
  mutate(region = paste0(country,"_",site))

all.models.corrected %>%
  group_by(region) %>%
  summarise(N = length(unique(model)),
            models = paste0(unique(model),collapse = "|"),
            .groups = "keep")

df.final <- all.models.corrected %>%
  group_by(model,region,country,site,age,type,var) %>%
  filter(time == max(time)) %>%
  mutate(type = tolower(type)) %>%
  mutate(site)


df.final.all <- bind_rows(df.final %>%
                            filter(type != "oldgrowth"),
                          df.final %>%
                            filter(type == "oldgrowth") %>%
                            mutate(age = -10),
                          df.final %>%
                            filter(type == "oldgrowth") %>%
                            mutate(age = -0.1),
                          df.final %>%
                            filter(type == "oldgrowth") %>%
                            mutate(age = -10),
                          df.final %>%
                            filter(type != "Oldgrowth") %>%
                            group_by(model,region,country,site,lat,lon,var) %>%
                            summarise(value = 0,
                                      age = -0.,
                                      .groups = "keep")) %>%
  filter(!is.na(region))

raw.data <- readRDS("/Users/felicien/Documents/projects/Yoko.regrowth/data/Data.sum.RDS")
# TODO: Combine with year_inventory

data <- raw.data %>%
  mutate(region = stri_trans_general(region, "Latin-ASCII")) %>%
  mutate(region = str_replace_all(str_trim(region), regex("[-\\s\\u00A0]+"), "_")) %>%
  filter(region %in% unique(paste0(df.final.all$region))) %>%
  mutate(age = case_when(age > 119 ~ -5,
                         TRUE ~ age))


df.final.all %>%
  group_by(model,region) %>%
  summarise(N = n()) %>%
  pivot_wider(names_from = "model",
              values_from = "N" ) %>%
  filter(ED2 != `ELM-FATES`)

df.final.all %>%
  filter(age > 0) %>%
  filter(model %in% c("ED2","ELM-FATES"),
         region %in% c("DRC_Djolu","Gh_Abofour",
                       "MOZ_GileNatPark","TAN_Kilwa",
                       "ZAM_Chintumukulu")) %>%
  group_by(model,region) %>%
  summarise(ages = paste0(sort(unique(age)),collapse = "_")) %>%
  arrange(region,model)


model.vs.data <- df.final.all %>%
  left_join(data %>%
              mutate(age = case_when(age == -5 ~ - 10,
                                     TRUE ~ age)),
            by = c("age","region")) %>%
  mutate(type = case_when(age == -10 ~ "Primary",
                          TRUE ~ "Secondary"))

MEM <- model.vs.data %>%
  group_by(region,site,country,type,age,var) %>%
  summarise(value = mean(value),
            AGC_ha.m = mean(AGC_ha.m),
            .groups = "keep") %>%
  filter(var == "cAGB")


A <- df.final.all %>%
  filter(var == "cAGB",
         model == "ELM-FATES")

sort(unique(A$region))

ggplot() +
  geom_line(data = df.final.all %>%
              filter(var == "cAGB"),
            aes(x = age,
                 y = value, color = model)) +
  geom_point(data = df.final.all %>%
               filter(var == "cAGB"),
             aes(x = age,
                y = value, color = model)) +
  geom_point(data = data,
             aes(x = age, y = AGC_ha.m/10)) +
  geom_errorbar(data = data,
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0) +
  facet_wrap(~region, scales = "free_y") +
  labs(x = "time since disturbance (yr)",
       y = "AGB (kgC/m²)",
       color = "Model") +
  theme_bw() #+
  # theme(legend.position = c(0.15,0.85)) +

saveRDS(df.final.all,
        "./outputs/All.modeloutputs.R")


################################################################################

all.ages <- seq(0,60,length.out = 1000)

all.fits <- data.frame()
all.models <- c("ED2","ELM-FATES","LPJGUESSv4.1.1","ORCHIDEEv2","ORCHIDEEv4.2")
for (cmodel in all.models){
  df <- df.final.all %>%
    filter(var == "cAGB",
           model == cmodel)

  secondary <- df %>%
    filter(type != "Primary") %>%
    filter(age > 0) %>%
    filter(grepl("DRC",region))

  fit <- nlsLM(value ~ a * ((1 - exp(-b * age))^ c),
               data = secondary ,
               start = list(a = 10, b = 0.02, c = 1),
               lower = c(a = 0, b = 0, c = 0),   # enforce positivity
               upper = c(a = 20, b = 1, c = 2))

  max.age <- ifelse(cmodel == "ORCHIDEEv4.2",Inf,300)

  fit.pred <- data.frame(age = c(-5,-0.1,all.ages),
                         AGB = predict(fit,
                                       newdata = data.frame(age = c(max.age,max.age,
                                                                    seq(0,60,length.out = 1000)))))


  all.fits <- bind_rows(all.fits,
                        fit.pred %>%
                          mutate(model = cmodel))


}

# Add data


fit <- nlsLM(value ~ a * ((1 - exp(-b * age))^ c),
             data = data %>%
               # filter(grepl("DRC",region)) %>%
               filter(age > 0) %>%
               mutate(value = AGC_ha.m/10),
             start = list(a = 10, b = 0.02, c = 1),
             lower = c(a = 0, b = 0, c = 0),   # enforce positivity
             upper = c(a = 20, b = 1, c = 2))


all.fits <- bind_rows(all.fits,
                      data.frame(age = c(-5,-0.1,all.ages),
                                 AGB = predict(fit,
                                               newdata = data.frame(age = c(300,300,
                                                                            seq(0,60,length.out = 1000))))) %>%
                        mutate(model = "data"))


df.Viola <- data.frame(age = all.ages) %>%
  mutate(AGB = 113/10*((1- exp(-0.02*age))**0.706))

ggplot() +
  geom_point(data = data %>%
               filter(grepl("DRC",region)),
             aes(x = age, y = AGC_ha.m/10),
             alpha = 0.3) +
  geom_errorbar(data = data %>%
                  filter(grepl("DRC",region)),
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0,
                alpha = 0.3) +
  geom_line(data = all.fits,
            aes(x = age, y = AGB, color = model)) +
  geom_line(data = df.Viola,
            aes(x = age, y = AGB), linetype = 2) +
  scale_color_manual(values = c("black","#e7298a","#d95f02","#7570b3","#1b9e77","#66a61e")) +
  scale_fill_manual(values = c("black","#e7298a","#d95f02","#7570b3","#1b9e77","#66a61e")) +
  theme_bw() +
  labs(x = "",y = "") +
  scale_y_continuous(limits = c(0,20)) +
  theme(text = element_text(size = 20)) +
  guides(color = "none")



ggplot() +
  geom_point(data = data ,
             aes(x = age, y = AGC_ha.m/10),
             alpha = 0.3) +
  geom_errorbar(data = data,
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0,
                alpha = 0.3) +
  geom_line(data = all.fits,
            aes(x = age, y = AGB, color = model)) +
  geom_line(data = df.Viola,
            aes(x = age, y = AGB), linetype = 2) +
  scale_color_manual(values = c("black","#e7298a","#d95f02","#7570b3","#1b9e77","#66a61e")) +
  scale_fill_manual(values = c("black","#e7298a","#d95f02","#7570b3","#1b9e77","#66a61e")) +
  theme_bw() +
  labs(x = "",y = "") +
  scale_y_continuous(limits = c(0,20)) +
  theme(text = element_text(size = 20)) +
  guides(color = "none")

data %>%
  filter(grepl("DRC",region)) %>%
  pull(region) %>%
  unique()



################################################################################

ggplot() +
  geom_line(data = df.final.all %>%
              filter(var == "cAGB"),
            aes(x = age,
                y = value,
                group = interaction(model,region),
                color = model)) +

  geom_point(data = df.final.all %>%
               filter(var == "cAGB"),
             aes(x = age,
                 y = value, color = model)) +
  geom_point(data = data,
             aes(x = age, y = AGC_ha.m/10)) +
  geom_errorbar(data = data,
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0) +
  facet_wrap(~model, scales = "free_y") +
  labs(x = "time since disturbance (yr)",
       y = "AGB (kgC/m²)",
       color = "Model") +
  theme_bw() +
  theme(legend.position = c(0.85,0.15)) +
  scale_y_continuous(limits = c(0,20))



ggplot() +
  geom_line(data = df.final.all %>%
              filter(var == "cAGB"),
            aes(x = age,
                y = value,
                group = interaction(model,region),
                color = model)) +

  geom_point(data = df.final.all %>%
               filter(var == "cAGB"),
             aes(x = age,
                 y = value, color = model)) +
  geom_point(data = data,
             aes(x = age, y = AGC_ha.m/10)) +
  geom_errorbar(data = data,
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0) +
  facet_wrap(~model, scales = "free_y") +
  labs(x = "time since disturbance (yr)",
       y = "AGB (kgC/m²)",
       color = "Model") +
  theme_bw() +
  theme(legend.position = c(0.85,0.15)) +
  scale_y_continuous(limits = c(0,20))



ggplot() +

  geom_line(data = MEM %>%
              filter(var == "cAGB"),
            aes(x = age,
                y = value), color = "black") +

  geom_point(data = data,
             aes(x = age, y = AGC_ha.m/10)) +
  geom_errorbar(data = data,
                aes(x = age, y = AGC_ha.m/10,
                    ymin = (AGC_ha.m-AGC_ha.sd)/10, ymax = (AGC_ha.m+AGC_ha.sd)/10),
                width = 0) +
  facet_wrap(~region, scales = "free_y") +
  labs(x = "time since disturbance (yr)",
       y = "AGB (kgC/m²)",
       color = "Model") +
  theme_bw() +
  theme(legend.position = c(0.15,0.85))

data2plot <-
  bind_rows(model.vs.data %>%
              filter(var == "cAGB"),
            MEM %>% mutate(model = "MEM")) %>%
  mutate(model = factor(model,
                        levels = c("ED2",
                                   "LPJGUESSv4.1.1",
                                   'ORCHIDEEv2',
                                   "ORCHIDEEv4.2",
                                   "ELM-FATES",
                                   "MEM")))

ggplot(data = data2plot %>%
         filter(region != "IC_Yaya") %>%
         filter(!(region == "LIB_Sapo" & (AGC_ha.m/10)>10))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(x = AGC_ha.m/10,
                 y = value, color = model),
             size = 1) +
  stat_smooth(aes(x = AGC_ha.m/10,
                  y = value,
                  color = model),
              method = "lm", se = FALSE) +
  labs(x = "Observed AGB (kgC/m²)",
       y = "Modelled AGB (kgC/m²)",
       color = "Model", shape = "Site") +
  facet_wrap(~ type, scales = "free") +
  scale_color_manual(values = c(scales::hue_pal()(length(models)),
                       "black")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = c(0.7,0.8))



ggplot(data = data2plot %>%
         filter(type == "Secondary",
                model != "MEM",
                region != "IC_Yaya") %>%
         filter(!(region == "LIB_Sapo" & (AGC_ha.m/10)>10))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(x = AGC_ha.m,
                 y = value*10, color = model),
             size = 1) +
  stat_smooth(aes(x = AGC_ha.m,
                  y = value*10,
                  color = model),
              method = "lm", se = FALSE) +
  labs(x = "",
       y = "",
       color = "Model", shape = "Site") +
  scale_color_manual(values = c(scales::hue_pal()(length(models)),
                                "black")) +
  theme_bw() +
  coord_equal() +
  theme(text = element_text(size = 20),strip.background = element_blank(),strip.text = element_blank(),
        legend.position = c(0.4,0.82))


data2plot %>%
  filter(type == "Secondary",
         model != "MEM",
         region != "IC_Yaya") %>%
  filter(!(region == "LIB_Sapo" & (AGC_ha.m/10)>10)) %>%
  group_by(model) %>%
  summarise(r.sq = summary(lm(formula = value ~ AGC_ha.m))[["r.squared"]])

ggplot(data2plot %>%
         filter(AGC_ha.m/20 > 1)) +
  geom_density(aes(x = 100*(value - AGC_ha.m/10)/(AGC_ha.m/10),
                   fill = model), alpha = 0.5) +
  facet_wrap(~ type, scales = "free") +
  theme_bw()

saveRDS(data2plot %>%
          mutate(obs = AGC_ha.m/10) %>%
          rename(pred = value) %>%
          dplyr::select(-c(AGC_ha.m,time,AGC_ha.sd)),
        "./outputs/Compiled.model.outputs.R")


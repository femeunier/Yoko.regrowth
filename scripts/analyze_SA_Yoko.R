rm(list = ls())

library(dplyr)
library(LidarED)
library(purrr)
library(stringr)
library(reshape2)

rundir <- "/user/scratchkyukon/gent/gvo000/gvo00074/felicien/ED2_soil/ED2/ED/run/SA"
outdir <- "/kyukon/scratch/gent/vo/000/gvo00074/felicien/Yoko/SA"

df.all <- df.param.all <- df.nplant <- data.frame()
params <- c("stomatal_slope","Vm0","Delta_Vm0","stoma_psi_b","Delta_stoma_psi_b","clumping_factor")
Nparam <- length(params)

quantiles <- c(0.025,0.975)

years2select <- c(5,12,20,60,250)
pfts <- c(2,3,4)

for (iparam in seq(1,Nparam)){
  for (iquantile in seq(1,length(quantiles))){

    param_name <- params[iparam]
    run_name <- paste0("SA_",param_name,"_",quantiles[iquantile])

    run_ref <- file.path(rundir,run_name)
    out_ref <- file.path(outdir,run_name)

    datum.file <- file.path(out_ref,"analy","analysis.RData")

    if (file.exists(datum.file)){

      load(datum.file)

      yr <- datum$year
      AGB <- datum$emean$agb
      AGB.early <- datum$szpft$agb[,12,2]
      LAI <- datum$emean$lai

      sel.years <- which((datum$year - datum$year[1]) %in% years2select)
      years <- datum$year[sel.years]

      nplant <- melt(datum$szpft$nplant[sel.years,2:11,pfts]) %>%
        rename(year = Var1,
               dbh.cat = Var2,
               pft = Var3,
               nplant = value) %>%
        mutate(year = years[year],
               pft = pfts[pft],
               nplant = nplant*10000,
               dbh.cat = dbh.cat + 1)

      config.file <-  file.path(run_ref,"config.xml")

      if (params[iparam] %in% c("Vm0","stoma_psi_b")) {
        df.param <- data.frame(pft = pfts,
                               param = params[iparam],
                               value = unlist(map(pfts, function(x)
                                 get_ED_default_pft(config.file, params[iparam], pft_target = x))))
      } else if (params[iparam] %in% c("Delta_Vm0","Delta_stoma_psi_b")){

        cparam <- substring(params[iparam],7)
        raw.values <- unlist(map(pfts, function(x)
          get_ED_default_pft(config.file, cparam, pft_target = x)))
        df.param <- data.frame(pft = pfts,
                               param = params[iparam],
                               value = rep(mean(diff(rev(raw.values))),3))
      } else {
        df.param <-
          data.frame(pft = pfts,
                     param = params[iparam],
                     value = rep(get_ED_default_pft(config.file, params[iparam], pft_target = 2), 3))
      }
      df.param.all <- bind_rows(list(df.param.all,
                                     df.param %>% mutate(quantile = quantiles[iquantile],
                                                         param = params[iparam])))

      df.all <- bind_rows(list(df.all,
                               data.frame(t = yr,
                                          agb = AGB,
                                          agb.early = AGB.early,
                                          LAI = LAI,
                                          quantile = quantiles[iquantile],
                                          param = params[iparam])))

      df.nplant <- bind_rows(list(df.nplant,
                                  nplant %>% mutate(quantile = quantiles[iquantile],
                                                    param = params[iparam])))

    }
  }
}


# Reference run

run_name <- "SA_reference"
out_ref <- file.path(outdir,run_name)

datum.file <- file.path(out_ref,"analy","analysis.RData")

if (file.exists(datum.file)){
  load(datum.file)

  yr <- datum$year
  AGB <- datum$emean$agb
  AGB.early <- datum$szpft$agb[,12,2]
  LAI <- datum$emean$lai

  sel.years <- which((datum$year - datum$year[1]) %in% years2select)
  years <- datum$year[sel.years]

  nplant <- melt(datum$szpft$nplant[sel.years,2:11,pfts]) %>%
    rename(year = Var1,
           dbh.cat = Var2,
           pft = Var3,
           nplant = value) %>%
    mutate(year = years[year],
           pft = pfts[pft],
           nplant = nplant*10000,
           dbh.cat = dbh.cat + 1)

  config.file <-  file.path(run_ref,"config.xml")

  for (iparam in seq(1,Nparam)){
    if (params[iparam] %in% c("Vm0","stoma_psi_b")) {
      df.param <- data.frame(pft = pfts,
                             param = params[iparam],
                             value = unlist(map(pfts, function(x)
                               get_ED_default_pft(config.file, params[iparam], pft_target = x))))
    } else if (params[iparam] %in% c("Delta_Vm0","Delta_stoma_psi_b")){

      cparam <- substring(params[iparam],7)
      raw.values <- unlist(map(pfts, function(x)
        get_ED_default_pft(config.file, cparam, pft_target = x)))
      df.param <- data.frame(pft = pfts,
                             param = params[iparam],
                             value = rep(mean(diff(rev(raw.values))),3))
    } else {
      df.param <-
        data.frame(pft = pfts,
                   param = params[iparam],
                   value = rep(get_ED_default_pft(config.file, params[iparam], pft_target = 2), 3))
    }
    df.param.all <- bind_rows(list(df.param.all,
                                   df.param %>% mutate(quantile = 0.5,
                                                       param = params[iparam])))
  }

  df.all <- bind_rows(list(df.all,
                           data.frame(t = yr,
                                      agb = AGB,
                                      agb.early = AGB.early,
                                      LAI = LAI,
                                      quantile = 0.5,
                                      param = "reference")))

  df.nplant <- bind_rows(list(df.nplant,
                              nplant %>% mutate(quantile = 0.5,
                                                param = "reference")))

}


saveRDS(object = df.all,file = file.path('.',"df_OP_SA_Yoko.RDS"))
saveRDS(object = df.param.all,file = file.path('.',"df_param_SA_Yoko.RDS"))
saveRDS(object = df.nplant,file = file.path('.',"df_nplant_Yoko.RDS"))

# scp /home/femeunier/Documents/projects/YGB/scripts/analyze_SA_Yoko.R hpc:/data/gent/vo/000/gvo00074/felicien/R


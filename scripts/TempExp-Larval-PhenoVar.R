#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

#rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(fullfact)
library(parallel)


#### LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.lk <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LK-Larvae")
larval.ls <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LS-Larvae")
larval.lo <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LO-Larvae")

# Combine each population, temperature, and species
larval <- bind_rows(larval.lk, larval.ls, larval.lo) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 2.2, 4.0, 4.4, 6.9, 8, 8.9),
                              labels = c("2.0", "2.2", "4.0", "4.4", "6.9", "8.0", "8.9")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         # Create a variable with population and species combined
         population = factor(interaction(population, species), ordered = TRUE,
                             levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                             labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         group = interaction(population, temperature)) %>% 
  rename(sire = male, dam = female)

## Clean up environment
rm(larval.lo, larval.ls, larval.lk)


# FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter out missing lengths
larval.tl <- larval %>% filter(!is.na(length_mm), length_mm != 0)

## filter out missing yolks
larval.yolk <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0)


#### STATISTICAL ANALYSIS - GENERATE OBSERVED VARIANCES ------------------------------------------

## Length-at-Hatch
phenoVar.tl.obs <- do.call(rbind, lapply(as.character(unique(larval.tl$group)), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.tl %>% filter(group == grp) %>% 
      select(family, dam, sire, block, length_mm)
    
    obs.tl <- observLmer(observ = data.grp, dam = "dam", sire = "sire", response = "length_mm")
    
    obs.tl.df <- data.frame(group = substr(grp, 1, nchar(grp)-4),
                            temperature = as.numeric(substr(grp, nchar(grp)-2, nchar(grp))),
                            dam.var = obs.tl$random[3,2],
                            dam.p = obs.tl$random[3,7],
                            dam.perc = obs.tl$random[3,3],
                            sire.var = obs.tl$random[2,2],
                            sire.p = obs.tl$random[2,7],
                            sire.perc = obs.tl$random[2,3],
                            dam.sire.var = obs.tl$random[1,2],
                            dam.sire.p = obs.tl$random[1,7],
                            dam.sire.perc = obs.tl$random[1,3],
                            residual.var = obs.tl$other[1,2],
                            residual.perc = obs.tl$other[1,3]) %>% 
      mutate_if(is.numeric, round, 4)  
})) %>% mutate(trait = "LAH")

## Yolk-sac Volume
phenoVar.yolk.obs <- do.call(rbind, lapply(as.character(unique(larval.yolk$group)), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.yolk %>% filter(group == grp) %>% 
    select(family, dam, sire, block, y_vol_mm3)
  
  obs.yolk <- observLmer2(observ = data.grp, dam = "dam", sire = "sire", response = "y_vol_mm3", block = "block")
  
  obs.yolk.df <- data.frame(group = substr(grp, 1, nchar(grp)-4),
                            temperature = as.numeric(substr(grp, nchar(grp)-2, nchar(grp))),
                            dam.var = obs.yolk$random[3,2],
                            dam.p = obs.yolk$random[3,7],
                            dam.perc = obs.yolk$random[3,3],
                            sire.var = obs.yolk$random[2,2],
                            sire.p = obs.yolk$random[2,7],
                            sire.perc = obs.yolk$random[2,3],
                            dam.sire.var = obs.yolk$random[1,2],
                            dam.sire.p = obs.yolk$random[1,7],
                            dam.sire.perc = obs.yolk$random[1,3],
                            residual.var = obs.yolk$other[1,2],
                            residual.perc = obs.yolk$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
})) %>% mutate(trait = "YSV")


# CREATE TEMPERATURE TREATMENT DATAFRAME ------------------------------------------------------

temp <- data.frame(group = c("LK-Whitefish", "LK-Whitefish", "LK-Whitefish", "LK-Whitefish",
                             "LK-Vendace", "LK-Vendace", "LK-Vendace", "LK-Vendace",
                             "LS-Cisco", "LS-Cisco", "LS-Cisco", "LS-Cisco",
                             "LO-Cisco", "LO-Cisco", "LO-Cisco", "LO-Cisco"),
                   temperature = c(rep(c(2.2, 4.0, 6.9, 8.0),2), rep(c(2.0, 4.4, 6.9, 8.9),2)),
                   temp.treatment = factor(rep(c("Coldest", "Cold", "Warm", "Warmest"), 4), 
                                           ordered = TRUE, levels = c("Coldest", "Cold", "Warm", "Warmest")))


#### CALCUALTE MEAN VARIANCE ACROSS TEMPERATURES -------------------------------------------------

phenoVar.larval.mean <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.mean = mean(dam.perc),
            sire.perc.mean = mean(sire.perc),
            dam.sire.perc.mean = mean(dam.sire.perc),
            residual.perc.mean = mean(residual.perc)) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "variance") %>% 
  filter(group != "LK-Whitefish") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.mean", "sire.perc.mean", "dam.sire.perc.mean", "residual.perc.mean"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, 
                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")))


#### CALCUALTE ERROR ACROSS TEMPERATURES ---------------------------------------------------------

phenoVar.larval.error <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
  group_by(group, trait) %>% 
  summarize(dam.perc.se = sd(dam.perc)/sqrt(n()),
            sire.perc.se = sd(sire.perc)/sqrt(n()),
            dam.sire.perc.se = sd(dam.sire.perc)/sqrt(n()),
            residual.perc.se = sd(residual.perc)/sqrt(n())) %>% 
  pivot_longer(3:6, names_to = "component", values_to = "error") %>% 
  filter(group != "LK-Whitefish") %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("dam.perc.se", "sire.perc.se", "dam.sire.perc.se", "residual.perc.se"),
                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
         group = factor(group, ordered = TRUE, 
                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")))


# JOIN MEAN AND ERROR -------------------------------------------------------------------------

phenoVar.larval.all <- left_join(phenoVar.larval.mean, phenoVar.larval.error)


#### COMBINE ALL TRAITS --------------------------------------------------------------------------

#phenoVar.larval.all <- bind_rows(phenoVar.tl.obs, phenoVar.yolk.obs) %>% 
#  left_join(temp) %>% 
#  filter(group != "LK-Whitefish") %>% 
#  select(group, temp.treatment, trait, dam.perc, sire.perc, dam.sire.perc, residual.perc) %>% 
#  pivot_longer(4:7, names_to = "component", values_to = "variance") %>% 
#  mutate(component = factor(component, ordered = TRUE,
#                            levels = c("dam.perc", "sire.perc", "dam.sire.perc", "residual.perc"),
#                            labels = c("Dam", "Sire", "Dam.Sire", "Error")),
#         component.trt = factor(interaction(component, temp.treatment), ordered = TRUE,
#                                levels = c("Dam.Coldest", "Dam.Cold", "Dam.Warm", "Dam.Warmest",
#                                           "Sire.Coldest", "Sire.Cold", "Sire.Warm", "Sire.Warmest",
#                                           "Dam.Sire.Coldest", "Dam.Sire.Cold", "Dam.Sire.Warm", "Dam.Sire.Warmest",
#                                           "Error.Coldest", "Error.Cold", "Error.Warm", "Error.Warmest")),
#         group = factor(group, ordered = TRUE, 
#                        levels = c("LK-Vendace", "LS-Cisco", "LO-Cisco"))) %>% 
#  left_join(phenoVar.larval.boot)


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## LAH
phenoVar.tl.cor <- phenoVar.tl.obs %>% 
  filter(group != "LK-Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))

## YSV
phenoVar.yolk.cor <- phenoVar.yolk.obs %>% 
  filter(group != "LK-Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(dam.cor.2 = ifelse(dam.cor >= 0.7, "POSITIVE", ifelse(dam.cor <= -0.7, "NEGATIVE", "NC")),
         sire.cor.2 = ifelse(sire.cor >= 0.7, "POSITIVE", ifelse(sire.cor <= -0.7, "NEGATIVE", "NC")),
         dam.sire.cor.2 = ifelse(dam.sire.cor >= 0.7, "POSITIVE", ifelse(dam.sire.cor <= -0.7, "NEGATIVE", "NC")),
         error.cor.2 = ifelse(error.cor >= 0.7, "POSITIVE", ifelse(error.cor <= -0.7, "NEGATIVE", "NC")))


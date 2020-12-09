#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(fullfact)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)


#### LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.lk <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LK-Larvae")
larval.ls <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LS-Larvae")
larval.lo <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LO-Larvae")

# Combine each population, temperature, and species
larval <- bind_rows(larval.lk, larval.ls, larval.lo) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 2.2, 4.0, 4.4, 6.9, 8, 8.9),
                              labels = c("2.0°C", "2.2°C", "4.0°C", "4.4°C", "6.9°C", "8.0°C", "8.9°C")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         # Create a variable with population and species combined
         population = factor(interaction(population, species), ordered = TRUE,
                             levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                             labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         group = gsub("-", ".", interaction(population, temperature)),
         group = gsub("°", "", group),
         group = gsub("\\.", "_", group)) %>% 
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
phenoVar.tl.obs <- do.call(rbind, lapply(unique(larval.tl$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.tl %>% filter(group == grp) %>% 
      select(family, dam, sire, block, length_mm)
    
    obs.tl <- observLmer(observ = data.grp, dam = "dam", sire = "sire", response = "length_mm")
    
    obs.tl.df <- data.frame(group = substr(grp, 1, nchar(grp)-5),
                            temperature = as.numeric(gsub("C", "", gsub("_", ".", substr(grp, nchar(grp)-3, nchar(grp))))),
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
}))

## Yolk-sac Volume
phenoVar.yolk.obs <- do.call(rbind, lapply(unique(larval.yolk$group), function(grp) {
  ## Filter to only a single temperature treatment
  data.grp <- larval.yolk %>% filter(group == grp) %>% 
    select(family, dam, sire, block, y_vol_mm3)
  
  obs.yolk <- observLmer2(observ = data.grp, dam = "dam", sire = "sire", response = "y_vol_mm3", block = "block")
  
  obs.yolk.df <- data.frame(group = substr(grp, 1, nchar(grp)-5),
                            temperature = as.numeric(gsub("C", "", gsub("_", ".", substr(grp, nchar(grp)-3, nchar(grp))))),
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
}))


#### CALCULATE MEANS ACROSS TEMPERATURES ---------------------------------------------------------

## LAH
phenoVar.tl.mean <- phenoVar.tl.obs %>% 
  filter(group != "LK_Whitefish") %>% 
  group_by(group) %>% 
  summarize(mean.dam.perc = mean(dam.perc),
            mean.sire.perc = mean(sire.perc),
            mean.dam.sire.perc = mean(dam.sire.perc),
            mean.error.perc = mean(residual.perc),
            sd.dam.perc = sd(dam.perc),
            sd.sire.perc = sd(sire.perc),
            sd.dam.sire.perc = sd(dam.sire.perc),
            sd.error.perc = sd(residual.perc)) %>% 
  mutate(trait = "tl")

## YSV
phenoVar.yolk.mean <- phenoVar.yolk.obs %>% 
  filter(group != "LK_Whitefish") %>% 
  group_by(group) %>% 
  summarize(mean.dam.perc = mean(dam.perc),
            mean.sire.perc = mean(sire.perc),
            mean.dam.sire.perc = mean(dam.sire.perc),
            mean.error.perc = mean(residual.perc),
            sd.dam.perc = sd(dam.perc),
            sd.sire.perc = sd(sire.perc),
            sd.dam.sire.perc = sd(dam.sire.perc),
            sd.error.perc = sd(residual.perc)) %>% 
  mutate(trait = "yolk")
  

#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

phenoVar.tl.cor <- heritability.tl.obs %>% 
  filter(group != "LK_Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate(dam.r2 = abs(dam.cor^2),
         sire.r2 = abs(sire.cor^2),
         dam.sire.r2 = abs(dam.sire.cor^2),
         error.r2 = abs(error.cor^2)) %>% 
  mutate_if(is.numeric, round, 2)

phenoVar.yolk.cor <- heritability.yolk.obs %>% 
  filter(group != "LK_Whitefish") %>% 
  group_by(group) %>% 
  summarize(dam.cor = cor(dam.perc, temperature),
            sire.cor = cor(sire.perc, temperature),
            dam.sire.cor = cor(dam.sire.perc, temperature),
            error.cor = cor(residual.perc, temperature)) %>% 
  mutate(dam.r2 = abs(dam.cor^2),
         sire.r2 = abs(sire.cor^2),
         dam.sire.r2 = abs(dam.sire.cor^2),
         error.r2 = abs(error.cor^2)) %>% 
  mutate_if(is.numeric, round, 2)


# COMBINE ALL TRAITS --------------------------------------------------------------------------

phenoVar.mean <- bind_rows(phenoVar.tl.mean, phenoVar.yolk.mean) %>% 
  pivot_longer(2:5, names_to = "component", values_to = "mean.var") %>% 
  select(-sd.dam.perc, -sd.sire.perc, -sd.dam.sire.perc, -sd.error.perc) %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("mean.dam.perc", "mean.sire.perc", "mean.dam.sire.perc", "mean.error.perc"),
                            labels = c("Dam", "Sire", "Dam:Sire", "Error")))

phenoVar.sd <- bind_rows(phenoVar.tl.mean, phenoVar.yolk.mean) %>% 
  pivot_longer(6:9, names_to = "component", values_to = "sd.var") %>% 
  select(-mean.dam.perc, -mean.sire.perc, -mean.dam.sire.perc, -mean.error.perc) %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("sd.dam.perc", "sd.sire.perc", "sd.dam.sire.perc", "sd.error.perc"),
                            labels = c("Dam", "Sire", "Dam:Sire", "Error")))

phenoVar.mean.sd <- left_join(phenoVar.mean, phenoVar.sd) %>% 
  mutate(group = factor(group, ordered = TRUE, levels = c("LK_Vendace", "LS_Cisco", "LO_Cisco"),
                        labels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")),
         trait = factor(trait, ordered = TRUE, levels = c("tl", "yolk"),
                        labels = c("Length-at-Hatch", "Yolk-sac Volume")))


#### VISUALIZATION - HERITABILITY --------------------------------------------

ggplot(phenoVar.mean.sd, aes(x = group, y = mean.var, group = component, fill = component)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = ifelse(mean.var - sd.var < 0, 0, mean.var - sd.var), 
                    ymax = ifelse(mean.var + sd.var > 100, 100, mean.var + sd.var)), 
                position = position_dodge(0.9), size = 0.8, width = 0.4, show.legend = FALSE) +
  scale_y_continuous(limits = c(-0.5, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"),
                    labels = c("Dam  ", "Sire  ", "Dam:Sire  ", "Error")) +
  labs(y = "Mean % of Total Phenotypic Variation", x = "Study Group") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 16),
        strip.background = element_rect(color = "white", fill = "white"),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  facet_wrap(~trait)

ggsave("figures/2020-Larvae-PhenoVar.png", width = 10, height = 8, dpi = 300)


#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES & SET THREADS -----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(fullfact)
library(parallel)
library(gridExtra)
library(grid)
library(cowplot)

setDTthreads(threads = 0)  # 0 = all available


#### LOAD INCUBATION TEMPERATURE DATA ------------------------------------------------------------

ADD <- read.csv("data/Artedi-Temperature-ADD-2020.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA --------------------------------------------------------------------------

hatch.USA <- read_excel("data/Coregonine-Temperature-Experiment-NA-Hatch.xlsx", sheet = "2020HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

hatch.Finland <- read_excel("data/Coregonine-Temperature-Experiment-FI-Hatch.xlsx", sheet = "2019HatchingData") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA, hatch.Finland) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
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
rm(hatch.USA, hatch.Finland, ADD)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1)

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1)


#### STATISTICAL ANALYSIS - GENERATE OBSERVED VARIANCES ------------------------------------------

## Embryo Survival
phenoVar.survival.obs <- do.call(rbind, lapply(unique(hatch.survival$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.survival %>% filter(group == grp) %>% 
      select(family, dam, sire, block, hatch)
  
  obs.survival <- observGlmer(observ = data.group, dam = "dam", sire = "sire", response = "hatch",
                               fam_link = binomial(logit))
  
  obs.survival.df <- data.frame(group = substr(grp, 1, nchar(grp)-5),
                                temperature = as.numeric(gsub("C", "", gsub("_", ".", substr(grp, nchar(grp)-3, nchar(grp))))),
                                dam.var = obs.survival$random[3,2],
                                dam.p = obs.survival$random[3,7],
                                dam.perc = obs.survival$random[3,3],
                                sire.var = obs.survival$random[2,2],
                                sire.p = obs.survival$random[2,7],
                                sire.perc = obs.survival$random[2,3],
                                dam.sire.var = obs.survival$random[1,2],
                                dam.sire.p = obs.survival$random[1,7],
                                dam.sire.perc = obs.survival$random[1,3],
                                residual.var = obs.survival$other[1,2],
                                residual.perc = obs.survival$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
}))

## DPF
phenoVar.dpf.obs <- do.call(rbind, lapply(unique(hatch.dpf$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.dpf %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf)
    
  obs.dpf <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "dpf")
  
  obs.dpf.df <- data.frame(group = substr(grp, 1, nchar(grp)-5),
                           temperature = as.numeric(gsub("C", "", gsub("_", ".", substr(grp, nchar(grp)-3, nchar(grp))))),
                           dam.var = obs.dpf$random[3,2],
                           dam.p = obs.dpf$random[3,7],
                           dam.perc = obs.dpf$random[3,3],
                           sire.var = obs.dpf$random[2,2],
                           sire.p = obs.dpf$random[2,7],
                           sire.perc = obs.dpf$random[2,3],
                           dam.sire.var = obs.dpf$random[1,2],
                           dam.sire.p = obs.dpf$random[1,7],
                           dam.sire.perc = obs.dpf$random[1,3],
                           residual.var = obs.dpf$other[1,2],
                           residual.perc = obs.dpf$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)
  }))

## ADD
phenoVar.ADD.obs <- do.call(rbind, lapply(unique(hatch.ADD$group), function(grp) {
  ## Filter to only a single group
  data.group <- hatch.ADD %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD)
  
  obs.ADD <- observLmer(observ = data.group, dam = "dam", sire = "sire", response = "ADD")
  
  obs.ADD.df <- data.frame(group = substr(grp, 1, nchar(grp)-5),
                           temperature = as.numeric(gsub("C", "", gsub("_", ".", substr(grp, nchar(grp)-3, nchar(grp))))),
                           dam.var = obs.ADD$random[3,2],
                           dam.p = obs.ADD$random[3,7],
                           dam.perc = obs.ADD$random[3,3],
                           sire.var = obs.ADD$random[2,2],
                           sire.p = obs.ADD$random[2,7],
                           sire.perc = obs.ADD$random[2,3],
                           dam.sire.var = obs.ADD$random[1,2],
                           dam.sire.p = obs.ADD$random[1,7],
                           dam.sire.perc = obs.ADD$random[1,3],
                           residual.var = obs.ADD$other[1,2],
                           residual.perc = obs.ADD$other[1,3]) %>% 
    mutate_if(is.numeric, round, 4)  
  }))


#### CALCULATE MEANS ACROSS TEMPERATURES ---------------------------------------------------------

## Embryo Survival
phenoVar.survival.mean <- phenoVar.survival.obs %>% 
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
  mutate(trait = "survival")

## DPF
phenoVar.dpf.mean <- phenoVar.dpf.obs %>% 
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
  mutate(trait = "dpf")

## ADD
phenoVar.ADD.mean <- phenoVar.ADD.obs %>% 
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
  mutate(trait = "ADD")


#### CALCULATE CORRELATIONS ----------------------------------------------------------------------

## Embryo Survival
phenoVar.survival.cor <- phenoVar.survival.obs %>% 
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

## DPF
phenoVar.dpf.cor <- phenoVar.dpf.obs %>% 
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

## ADD
phenoVar.ADD.cor <- phenoVar.ADD.obs %>% 
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

phenoVar.mean <- bind_rows(phenoVar.survival.mean, phenoVar.dpf.mean, phenoVar.ADD.mean) %>% 
  pivot_longer(2:5, names_to = "component", values_to = "mean.var") %>% 
  select(-sd.dam.perc, -sd.sire.perc, -sd.dam.sire.perc, -sd.error.perc) %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("mean.dam.perc", "mean.sire.perc", "mean.dam.sire.perc", "mean.error.perc"),
                            labels = c("Dam", "Sire", "Dam:Sire", "Error")))

phenoVar.sd <- bind_rows(phenoVar.survival.mean, phenoVar.dpf.mean, phenoVar.ADD.mean) %>% 
  pivot_longer(6:9, names_to = "component", values_to = "sd.var") %>% 
  select(-mean.dam.perc, -mean.sire.perc, -mean.dam.sire.perc, -mean.error.perc) %>% 
  mutate(component = factor(component, ordered = TRUE,
                            levels = c("sd.dam.perc", "sd.sire.perc", "sd.dam.sire.perc", "sd.error.perc"),
                            labels = c("Dam", "Sire", "Dam:Sire", "Error")))

phenoVar.mean.sd <- left_join(phenoVar.mean, phenoVar.sd) %>% 
  mutate(group = factor(group, ordered = TRUE, levels = c("LK_Vendace", "LS_Cisco", "LO_Cisco"),
                        labels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")),
         trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")))


#### VISUALIZATION -------------------------------------------------------------------------------

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

ggsave("figures/2020-Embryo-PhenoVar.png", width = 14, height = 8, dpi = 300)



# CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


# SET RANDOM SEED FOR REPRODUCIBILITY -------------------------------------

set.seed(897231876)


# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggplot2)
library(lme4)
library(lmerTest)


# LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD.2020 <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


# LOAD HATCHING DATA ------------------------------------------------------

hatch.USA.2020 <- read_excel("data/2020-Artedi-Temperature-Experiment.xlsx", sheet = "2020HatchingData") %>% 
  mutate(year = 2020) %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD.2020) %>% 
  dplyr::select(year, population, species, family = family_herit, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.albula <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi vendace") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family = family_herit, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.lavaretus <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi whitefish") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family = family_herit, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 4.5, 7, 9),
                              labels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         # Create a variable with population and species combined
         population = factor(interaction(population, species), ordered = TRUE,
                             levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                             labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco"))) %>% 
  rename(sire = male, dam = female, tray = plate) %>% 
  group_by(year, population, species, temperature, family) %>% 
  mutate(rep = seq(1, n(), 1)) %>% ungroup()

## Clean up environment
rm(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus, ADD.2020)


# STATISTICAL ANALYSIS - SURVIVAL - HERITABILITY --------------------------

# filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.survival.summary <- do.call(rbind, lapply(unique(hatch.survival$temperature), function(j) {
  
  ## Filter to only a single temperature treatment
  data.temp <- hatch.survival %>% filter(temperature == j)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(k) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, population == k) %>% 
      select(dam, sire, block, tray, hatch)
    
    ## Fit standard random-effect model
    herit.survival.model <- observGlmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "hatch",
                                         block = "block", position = "tray", fam_link = binomial(logit))
    
    ## Create a vector of variances
    var.tray <- herit.survival.model$random[1,2]
    var.dom <- herit.survival.model$random[2,2]
    var.sire <- herit.survival.model$random[3,2]
    var.dam <- herit.survival.model$random[4,2]
    var.block <- herit.survival.model$random[5,2]
    var.resid <- herit.survival.model$other[1,2]
    var.total <- var.tray + var.dom + var.sire + var.dam + var.block + var.resid
    
    ## Create a vector of p-values
    p.tray <- herit.survival.model$random[1,7]
    p.dom <- herit.survival.model$random[2,7]
    p.sire <- herit.survival.model$random[3,7]
    p.dam <- herit.survival.model$random[4,7]
    p.block <- herit.survival.model$random[5,7]
    
    ## Create a data frame with variances and p-values
    variances.survival.df <- data.frame(population = k, temperature = j, trait = "surv", var.tray = var.tray, p.tray = p.tray,
                                   var.dom = var.dom, p.dom = p.dom, var.sire = var.sire, p.sire = p.sire,
                                   var.dam = var.dam, p.dam = p.dam, var.block = var.block, p.block = p.block,
                                   var.resid = var.resid, var.total = var.total, h2 = round(var.sire / var.total, 4))
  }))
}))


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - HERITABILITY -----------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf))

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.dpf.summary <- do.call(rbind, lapply(unique(hatch.dpf$temperature), function(j) {
  
  ## Filter to only a single temperature treatment
  data.temp <- hatch.dpf %>% filter(temperature == j)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(k) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, population == k) %>% 
      select(dam, sire, block, tray, dpf)
    
    ## Fit standard random-effect model
    herit.dpf.model <- observLmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "dpf",
                                   block = "block", position = "tray")
    
    ## Create a vector of variances
    var.tray <- herit.dpf.model$random[1,2]
    var.dom <- herit.dpf.model$random[2,2]
    var.sire <- herit.dpf.model$random[3,2]
    var.dam <- herit.dpf.model$random[4,2]
    var.block <- herit.dpf.model$random[5,2]
    var.resid <- herit.dpf.model$other[1,2]
    var.total <- var.tray + var.dom + var.sire + var.dam + var.block + var.resid
    
    ## Create a vector of p-values
    p.tray <- herit.dpf.model$random[1,7]
    p.dom <- herit.dpf.model$random[2,7]
    p.sire <- herit.dpf.model$random[3,7]
    p.dam <- herit.dpf.model$random[4,7]
    p.block <- herit.dpf.model$random[5,7]
    
    ## Create a data frame with variances and p-values
    variances.dpf.df <- data.frame(population = k, temperature = j, trait = "dpf", var.tray = var.tray, p.tray = p.tray,
                                   var.dom = var.dom, p.dom = p.dom, var.sire = var.sire, p.sire = p.sire,
                                   var.dam = var.dam, p.dam = p.dam, var.block = var.block, p.block = p.block,
                                   var.resid = var.resid, var.total = var.total, h2 = round(var.sire / var.total, 4))
  }))
}))


# STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - HERITABILITY -----------

# filter to only hatched embryos
hatch.add <- hatch %>% filter(!is.na(ADD))

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.add.summary <- do.call(rbind, lapply(unique(hatch.add$temperature), function(j) {
  
  ## Filter to only a single temperature treatment
  data.temp <- hatch.add %>% filter(temperature == j)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(k) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, population == k) %>% 
      select(dam, sire, block, tray, ADD)
    
    ## Fit standard random-effect model
    herit.add.model <- observLmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "ADD",
                                   block = "block", position = "tray")
    
    ## Create a vector of variances
    var.tray <- herit.add.model$random[1,2]
    var.dom <- herit.add.model$random[2,2]
    var.sire <- herit.add.model$random[3,2]
    var.dam <- herit.add.model$random[4,2]
    var.block <- herit.add.model$random[5,2]
    var.resid <- herit.add.model$other[1,2]
    var.total <- var.tray + var.dom + var.sire + var.dam + var.block + var.resid
    
    ## Create a vector of p-values
    p.tray <- herit.add.model$random[1,7]
    p.dom <- herit.add.model$random[2,7]
    p.sire <- herit.add.model$random[3,7]
    p.dam <- herit.add.model$random[4,7]
    p.block <- herit.add.model$random[5,7]
    
    ## Create a data frame with variances and p-values
    variances.add.df <- data.frame(population = k, temperature = j, trait = "ADD", var.tray = var.tray, p.tray = p.tray,
                                   var.dom = var.dom, p.dom = p.dom, var.sire = var.sire, p.sire = p.sire,
                                   var.dam = var.dam, p.dam = p.dam, var.block = var.block, p.block = p.block,
                                   var.resid = var.resid, var.total = var.total, h2 = round(var.sire / var.total, 4))
  }))
}))


# VISUALIZATION - HERITABILITY --------------------------------------------

heritability.all <- bind_rows(heritability.survival.summary, heritability.add.summary, heritability.dpf.summary) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("surv", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")),
         label = ifelse(h2 == 0, 0, NA))
write.csv(heritability.all, "data/embryo-heritability.csv", row.names = FALSE)

ggplot(heritability.all, aes(x = temperature, y = (h2 * 100), group = population, fill = population)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  geom_text(aes(x = temperature, y = 0.75, label = label), size = 5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymin = (mean.herit - se.herit) * 100, ymax = (mean.herit + se.herit) * 100), 
  #              position = position_dodge(0.9),
  #              size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_fill_grey("combine", start = 0.2, end = 0.9,
                  labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-0.5, 60), breaks = seq(0, 60, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (%)") +
  theme_bw() + 
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 15),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) + 
  facet_wrap(~trait, nrow = 1)

ggsave("figures/embryo/2020-Heritability-fullfact.png", width = 20, height = 12, dpi = 300)


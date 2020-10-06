# CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


# SET RANDOM SEED FOR REPRODUCIBILITY -------------------------------------

set.seed(897231876)


# LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggplot2)
library(fullfact)


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
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.albula <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi vendace") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.lavaretus <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi whitefish") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

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
  rename(sire = male, dam = female, tray = plate)

## Clean up environment
rm(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus, ADD.2020)


# LOAD PREVIOUS BOOTSTRAPPED DISTRIBUTIONS IF FOUND -----------------------
ifelse(file.exists("data/Heritability_Bootstrap/heritability_survival_boot.csv") == TRUE, 
       heritability.survival.boot <- read.csv("data/Heritability_Bootstrap/heritability_survival_boot.csv", header = TRUE), NA)
ifelse(file.exists("data/Heritability_Bootstrap/heritability_dpf_boot.csv") == TRUE, 
       heritability.dpf.boot <- read.csv("data/Heritability_Bootstrap/heritability_dpf_boot.csv", header = TRUE), NA)
ifelse(file.exists("data/Heritability_Bootstrap/heritability_ADD_boot.csv") == TRUE, 
       heritability.ADD.boot <- read.csv("data/Heritability_Bootstrap/heritability_ADD_boot.csv", header = TRUE), NA)


# STATISTICAL ANALYSIS - SURVIVAL - HERITABILITY --------------------------

# filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.survival.boot <- do.call(rbind, lapply(unique(hatch.survival$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.survival %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, tray, hatch.raw = hatch)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:2500, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.temp.group$family), function(fam) {
        data.family <- data.temp.group %>% filter(family == fam)
        hatch.boot <- sample(data.family$hatch.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, hatch.boot) %>% select(-hatch.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block),
               tray = factor(tray))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), 
                                 paste0("sire", length), paste0("block", length), 
                                 paste0("tray", length), paste0("hatch", length))
      data.family
    }))
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampGlmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "hatch",
                                          position = "tray", block = "block", fam_link = binomial(logit), start = 1, end = 2500)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(population = pop, temperature = temp) %>% 
      select(population, temperature, sire, dam, dam.sire, tray, block, residual = Residual, additive, nonadd, maternal)
  }))
})) %>% 
  mutate(trait = "survival")

## Save the bootstrapped data for future use
write.csv(heritability.survival.boot, "data/Heritability_Bootstrap/heritability_survival_boot.csv", row.names = FALSE)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - HERITABILITY -----------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf))

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.dpf.boot <- do.call(rbind, lapply(unique(hatch.dpf$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.dpf %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, tray, dpf.raw = dpf)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:2500, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.temp.group$family), function(fam) {
        data.family <- data.temp.group %>% filter(family == fam)
        dpf.boot <- sample(data.family$dpf.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, dpf.boot) %>% select(-dpf.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block),
               tray = factor(tray))
      
      colnames(data.family) <- c(paste0("temperature", length), paste0("dam", length), 
                                 paste0("sire", length), paste0("block", length), 
                                 paste0("tray", length), paste0("dpf", length))
      data.family
    }))
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.lmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "dpf",
                                        position = "tray", block = "block", start = 1, end = 2500)
    
    bootstrap.data.lmer2.h2 <- data.frame(bootstrap.data.lmer2) %>% 
      mutate(population = pop, temperature = temp) %>% 
      select(population, temperature, sire, dam, dam.sire, tray, block, residual = Residual, additive, nonadd, maternal)
  }))
})) %>% 
  mutate(trait = "dpf")
  
## Save the bootstrapped data for future use
write.csv(heritability.dpf.boot, "data/Heritability_Bootstrap/heritability_dpf_boot.csv", row.names = FALSE)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - HERITABILITY -----------

# filter to only hatched embryos
hatch.add <- hatch %>% filter(!is.na(ADD))

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.ADD.boot <- do.call(rbind, lapply(unique(hatch.add$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.add %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, tray, ADD.raw = ADD)
      
      ## Create a bootstrapped data set from each temperature and group
      bootstrap.data <- do.call(cbind, lapply(1:2500, function(length) {
        ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
        data.family <- do.call(rbind, lapply(unique(data.temp.group$family), function(fam) {
          data.family <- data.temp.group %>% filter(family == fam)
          ADD.boot <- sample(data.family$ADD.raw, replace = T, size = nrow(data.family))
          data.family.boot <- data.frame(data.family, ADD.boot) %>% select(-ADD.raw)
        })) %>% 
          mutate(family = factor(family),
                 dam = factor(dam),
                 sire = factor(sire),
                 block = factor(block),
                 tray = factor(tray))
        
        colnames(data.family) <- c(paste0("temperature", length), paste0("dam", length), 
                                   paste0("sire", length), paste0("block", length), 
                                   paste0("tray", length), paste0("ADD", length))
        data.family
      }))

      ## Calculate variance components from bootstrapped sample
      bootstrap.data.lmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "ADD",
                                          position = "tray", block = "block", start = 1, end = 2500)
      
      bootstrap.data.lmer2.h2 <- data.frame(bootstrap.data.lmer2) %>% 
        mutate(population = pop, temperature = temp) %>% 
        select(population, temperature, sire, dam, dam.sire, tray, block, residual = Residual, additive, nonadd, maternal)
  }))
})) %>% 
  mutate(trait = "ADD")

## Save the bootstrapped data for future use
write.csv(heritability.ADD.boot, "data/Heritability_Bootstrap/heritability_ADD_boot.csv", row.names = FALSE)


# CALCULATE THE VARIANCES FROM RAW DATA -----------------------------------

heritability.raw <- do.call(rbind, lapply(unique(hatch$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch %>% filter(temperature == temp)

  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, tray, eye, hatch, dpf, ADD)

    ## Filter to each trait
    hatch.survival <- data.temp.group %>% filter(eye != 0)
    hatch.dpf <- data.temp.group %>% filter(!is.na(dpf))
    hatch.add <- data.temp.group %>% filter(!is.na(ADD))
    
    ## Calculate variance components from bootstrapped sample
    data.hatch.glmer2 <- observGlmer2(observ = hatch.survival, dam = "dam", sire = "sire", response = "hatch",
                                      position = "tray", block = "block", fam_link = binomial(logit))
    data.dpf.lmer2 <- observLmer2(observ = hatch.dpf, dam = "dam", sire = "sire", response = "dpf",
                                  position = "tray", block = "block")
    data.add.lmer2 <- observLmer2(observ = hatch.add, dam = "dam", sire = "sire", response = "ADD",
                                  position = "tray", block = "block")
    
    ## Convert to data frame
    data.survival.h2 <- data.frame(population = pop,
                                   temperature = temp,
                                   trait = "survival",
                                   sire = data.hatch.glmer2$random[3,2],
                                   dam = data.hatch.glmer2$random[4,2],
                                   dam.sire = data.hatch.glmer2$random[2,2],
                                   tray = data.hatch.glmer2$random[1,2],
                                   block = data.hatch.glmer2$random[5,2],
                                   residual = data.hatch.glmer2$other[1,2],
                                   additive = data.hatch.glmer2$calculation[1,2],
                                   nonadd = data.hatch.glmer2$calculation[2,2],
                                   maternal = data.hatch.glmer2$calculation[3,2])
    
    data.dpf.h2 <- data.frame(population = pop,
                              temperature = temp,
                              trait = "dpf",
                              sire = data.dpf.lmer2$random[3,2],
                              dam = data.dpf.lmer2$random[4,2],
                              dam.sire = data.dpf.lmer2$random[2,2],
                              tray = data.dpf.lmer2$random[1,2],
                              block = data.dpf.lmer2$random[5,2],
                              residual = data.dpf.lmer2$other[1,2],
                              additive = data.dpf.lmer2$calculation[1,2],
                              nonadd = data.dpf.lmer2$calculation[2,2],
                              maternal = data.dpf.lmer2$calculation[3,2])

    data.ADD.h2 <- data.frame(population = pop,
                              temperature = temp,
                              trait = "ADD",
                              sire = data.add.lmer2$random[3,2],
                              dam = data.add.lmer2$random[4,2],
                              dam.sire = data.add.lmer2$random[2,2],
                              tray = data.add.lmer2$random[1,2],
                              block = data.add.lmer2$random[5,2],
                              residual = data.add.lmer2$other[1,2],
                              additive = data.add.lmer2$calculation[1,2],
                              nonadd = data.add.lmer2$calculation[2,2],
                              maternal = data.add.lmer2$calculation[3,2])
    
    data.h2 <- bind_rows(data.survival.h2, data.dpf.h2, data.ADD.h2)
  }))
})) 

heritability.raw.select <- heritability.raw %>% 
  mutate(pheno = additive + nonadd + residual + block + tray,
         h2.raw = additive / pheno) %>% 
  select(population, temperature, trait, h2.raw)


# CALCULATE THE MEAN, SE, AND CI FROM BOOTSTRAPPED DISTRIBUTIONS ----------

## Embryo Survival
heritability.survival.summary <- heritability.survival.boot %>% 
  mutate(pheno = additive + nonadd + residual + block + tray,
         h2 = additive / pheno) %>% 
  group_by(population, temperature, trait) %>% 
  summarize(mean.herit = mean(h2),
            se.herit = sd(h2),
            herit.upper.CI = quantile(h2, 0.975),
            herit.lower.CI = quantile(h2, 0.025),
            mean.var.add = mean(additive),
            se.var.add = sd(additive),
            var.add.upper.CI = quantile(additive, 0.975),
            var.add.lower.CI = quantile(additive, 0.025),
            mean.var.pheno = mean(pheno),
            se.var.pheno = sd(pheno),
            var.pheno.upper.CI = quantile(pheno, 0.975),
            var.pheno.lower.CI = quantile(pheno, 0.025))

## DPF
heritability.dpf.summary <- heritability.dpf.boot %>% 
  mutate(pheno = additive + nonadd + residual + block + tray,
         h2 = additive / pheno) %>% 
  group_by(population, temperature, trait) %>% 
  summarize(mean.herit = mean(h2),
            se.herit = sd(h2),
            herit.upper.CI = quantile(h2, 0.975),
            herit.lower.CI = quantile(h2, 0.025),
            mean.var.add = mean(additive),
            se.var.add = sd(additive),
            var.add.upper.CI = quantile(additive, 0.975),
            var.add.lower.CI = quantile(additive, 0.025),
            mean.var.pheno = mean(pheno),
            se.var.pheno = sd(pheno),
            var.pheno.upper.CI = quantile(pheno, 0.975),
            var.pheno.lower.CI = quantile(pheno, 0.025))

## ADD
heritability.ADD.summary <- heritability.ADD.boot %>% 
  mutate(pheno = additive + nonadd + residual + block + tray,
         h2 = additive / pheno) %>% 
  group_by(population, temperature, trait) %>% 
  summarize(mean.herit = mean(h2),
            se.herit = sd(h2),
            herit.upper.CI = quantile(h2, 0.975),
            herit.lower.CI = quantile(h2, 0.025),
            mean.var.add = mean(additive),
            se.var.add = sd(additive),
            var.add.upper.CI = quantile(additive, 0.975),
            var.add.lower.CI = quantile(additive, 0.025),
            mean.var.pheno = mean(pheno),
            se.var.pheno = sd(pheno),
            var.pheno.upper.CI = quantile(pheno, 0.975),
            var.pheno.lower.CI = quantile(pheno, 0.025))


# VISUALIZATION - HERITABILITY --------------------------------------------

heritability.all <- bind_rows(heritability.survival.summary, heritability.ADD.summary, heritability.dpf.summary) %>% 
  #left_join(heritability.raw.select) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")))
#write.csv(heritability.all, "data/Heritability_Bootstrap/Heritability_all.csv", row.names = FALSE)

ggplot(heritability.all, aes(x = temperature, y = (mean.herit * 100), group = population, fill = population)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  #geom_text(aes(x = temperature, y = 0.75, label = label), size = 5, position = position_dodge(width = 0.9)) +
  #geom_errorbar(aes(ymin = (mean.herit - se.herit) * 100, ymax = (mean.herit + se.herit) * 100), 
  #              position = position_dodge(0.9),
  #              size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  geom_errorbar(aes(ymin = (herit.upper.CI * 100), ymax = (herit.lower.CI * 100)), 
                position = position_dodge(0.9),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_fill_grey("combine", start = 0.2, end = 0.9,
                  labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-0.5, 100), breaks = seq(0, 100, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (% with 95% CI)") +
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

ggsave("figures/embryo/2020-Heritability-fullfact-2500.png", width = 20, height = 12, dpi = 300)


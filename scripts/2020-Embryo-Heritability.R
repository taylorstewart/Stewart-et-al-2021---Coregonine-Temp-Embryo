# CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


# SET RANDOM SEED FOR REPRODUCIBILITY -------------------------------------

set.seed(897231876)  ## 0 - 10,000


# LOAD PACKAGES & SET THREADS -----------------------------------------------------------------

library(tidyverse)
library(readxl)
library(data.table)
library(ggplot2)
library(fullfact)
library(parallel)

setDTthreads(threads = 0)  # 0 = all available


# LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


# LOAD HATCHING DATA ------------------------------------------------------

hatch.USA <- read_excel("data/2020-Artedi-Temperature-Experiment.xlsx", sheet = "2020HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

hatch.Finland <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, species, family, male, female, block, temperature, eye, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA, hatch.Finland) %>% 
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
                             labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         group = gsub("-", ".", interaction(population, temperature)),
         group = gsub("°", "", group),
         group = gsub("\\.", "_", group)) %>% 
  rename(sire = male, dam = female)

## Clean up environment
rm(hatch.USA, hatch.Finland, ADD)


# FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1)

## filter to only hatched embryos
hatch.add <- hatch %>% filter(!is.na(ADD), hatch == 1)


# STATISTICAL ANALYSIS - SURVIVAL - HERITABILITY --------------------------

## Run nested loop to create a bootstrapped dataset for each temperature and population
if(file.exists("data/Heritability_Bootstrap/heritability_survival_boot_fish.csv") == FALSE) {

  heritability.survival.boot.fish <- do.call(rbind, mclapply(unique(hatch.survival$group), mc.cores = 8, function(grp) {
    ## Filter to only a single temperature treatment
    data.filt <- hatch.survival %>% filter(group == grp) %>% 
        select(family, dam, sire, block, hatch.raw = hatch)
      
      ## Create a bootstrapped data set from each temperature and group
      bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
        ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
        data.family <- do.call(rbind, lapply(unique(data.filt$family), function(fam) {
          data.family <- data.filt %>% filter(family == fam)
          hatch.boot <- sample(data.family$hatch.raw, replace = T, size = nrow(data.family))
          data.family.boot <- data.frame(data.family, hatch.boot) %>% select(-hatch.raw)
        })) %>% 
          mutate(family = factor(family),
                 dam = factor(dam),
                 sire = factor(sire),
                 block = factor(block))
        
        colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                   paste0("block", length), paste0("hatch", length))
        data.family
      })) %>% 
        transmute(group = grp, !!!.)
      
  }))
  
  ## Save the bootstrapped fish for future use
  write.csv(heritability.survival.boot.fish, "data/Heritability_Bootstrap/heritability_survival_boot_fish.csv", row.names = FALSE)
} else {
  heritability.survival.boot.fish <- fread("data/Heritability_Bootstrap/heritability_survival_boot_fish.csv")
}

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
if(file.exists("data/Heritability_Bootstrap/heritability_survival_boot_10000.csv") == FALSE) {
  
  heritability.survival.boot <- do.call(rbind, mclapply(unique(heritability.survival.boot.fish$group), mc.cores = 8, function(grp) {
    ## Filter to only a single temperature treatment
    data.group <- heritability.survival.boot.fish %>% filter(group == grp)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.glmer2 <- resampGlmer2(resamp = data.group, dam = "dam", sire = "sire", response = "hatch",
                                          block = "block", fam_link = binomial(logit), start = 1, end = 10000)
    
    bootstrap.data.glmer2.h2 <- data.frame(bootstrap.data.glmer2) %>% 
      mutate(group = grp, rep = 1:n()) %>% 
      select(group, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  })) %>% 
    mutate(trait = "survival")
  
  ## Save the bootstrapped data for future use
  write.csv(heritability.survival.boot, "data/Heritability_Bootstrap/heritability_survival_boot_10000.csv", row.names = FALSE)
} else {
  heritability.survival.boot <- fread("data/Heritability_Bootstrap/heritability_survival_boot_10000.csv")
}
rm(heritability.survival.boot.fish)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - HERITABILITY -----------

## Run nested loop to create a bootstrapped dataset for each temperature and population
if(file.exists("data/Heritability_Bootstrap/heritability_dpf_boot_fish.csv") == FALSE) {
  
  heritability.dpf.boot.fish <- do.call(rbind, mclapply(unique(hatch.dpf$group), mc.cores = 8, function(grp) {
    ## Filter to only a single temperature treatment
    data.filt <- hatch.dpf %>% filter(group == grp) %>% 
      select(family, dam, sire, block, dpf.raw = dpf)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.filt$family), function(fam) {
        data.family <- data.filt %>% filter(family == fam)
        dpf.boot <- sample(data.family$dpf.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, dpf.boot) %>% select(-dpf.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("dpf", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
  }))
  
  ## Save the bootstrapped fish for future use
  write.csv(heritability.dpf.boot.fish, "data/Heritability_Bootstrap/heritability_dpf_boot_fish.csv", row.names = FALSE)
} else {
  heritability.dpf.boot.fish <- fread("data/Heritability_Bootstrap/heritability_dpf_boot_fish.csv")
}

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
if(file.exists("data/Heritability_Bootstrap/heritability_dpf_boot_10000.csv") == FALSE) {
  
  heritability.dpf.boot <- do.call(rbind, mclapply(unique(heritability.dpf.boot.fish$group), mc.cores = 8, function(temp) {
    ## Filter to only a single population and temperature treatment
    data.group <- heritability.dpf.boot.fish %>% filter(group == grp)
    
      ## Calculate variance components from bootstrapped sample
      bootstrap.data.lmer2 <- resampLmer2(resamp = data.group, dam = "dam", sire = "sire", response = "dpf",
                                          block = "block", start = 1, end = 10000)
      
      bootstrap.data.lmer2.h2 <- data.frame(bootstrap.data.lmer2) %>% 
        mutate(group = grp, rep = 1:n()) %>% 
        select(group, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  })) %>% 
    mutate(trait = "dpf")
  
  ## Save the bootstrapped data for future use
  write.csv(heritability.dpf.boot, "data/Heritability_Bootstrap/heritability_dpf_boot_10000.csv", row.names = FALSE)
} else {
  heritability.dpf.boot <- fread("data/Heritability_Bootstrap/heritability_dpf_boot_10000.csv")
}
rm(heritability.dpf.boot.fish)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - HERITABILITY -----------

## Run nested loop to create a bootstrapped dataset for each temperature and population
if(file.exists("data/Heritability_Bootstrap/heritability_ADD_boot_fish.csv") == FALSE) {
  
  heritability.ADD.boot.fish <- do.call(rbind, mclapply(unique(hatch.add$group), mc.cores = 8, function(grp) {
    ## Filter to only a single temperature treatment
    data.filt <- hatch.add %>% filter(group == grp) %>% 
      select(family, dam, sire, block, ADD.raw = ADD)
    
    ## Create a bootstrapped data set from each temperature and group
    bootstrap.data <- do.call(cbind, lapply(1:10000, function(length) {
      ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
      data.family <- do.call(rbind, lapply(unique(data.filt$family), function(fam) {
        data.family <- data.filt %>% filter(family == fam)
        ADD.boot <- sample(data.family$ADD.raw, replace = T, size = nrow(data.family))
        data.family.boot <- data.frame(data.family, ADD.boot) %>% select(-ADD.raw)
      })) %>% 
        mutate(family = factor(family),
               dam = factor(dam),
               sire = factor(sire),
               block = factor(block))
      
      colnames(data.family) <- c(paste0("family", length), paste0("dam", length), paste0("sire", length), 
                                 paste0("block", length), paste0("ADD", length))
      data.family
    })) %>% 
      transmute(group = grp, !!!.)
    
  }))
  
  ## Save the bootstrapped fish for future use
  write.csv(heritability.ADD.boot.fish, "data/Heritability_Bootstrap/heritability_ADD_boot_fish.csv", row.names = FALSE)
} else {
  heritability.ADD.boot.fish <- fread("data/Heritability_Bootstrap/heritability_ADD_boot_fish.csv")
}

if(file.exists("data/Heritability_Bootstrap/heritability_ADD_boot_10000.csv") == FALSE) {
  
  heritability.ADD.boot <- do.call(rbind, mclapply(unique(heritability.ADD.boot.fish$group), mc.cores = 8, function(temp) {
    ## Filter to only a single population and temperature treatment
    data.group <- heritability.ADD.boot.fish %>% filter(group == grp)
    
    ## Calculate variance components from bootstrapped sample
    bootstrap.data.lmer2 <- resampLmer2(resamp = data.group, dam = "dam", sire = "sire", response = "ADD",
                                        block = "block", start = 1, end = 10000)
    
    bootstrap.data.lmer2.h2 <- data.frame(bootstrap.data.lmer2) %>% 
      mutate(group = grp, rep = 1:n()) %>% 
      select(group, rep, sire, dam, dam.sire, block, residual = Residual, total = Total, additive, nonadd, maternal)
  })) %>% 
    mutate(trait = "ADD")
  
  ## Save the bootstrapped data for future use
  write.csv(heritability.ADD.boot, "data/Heritability_Bootstrap/heritability_ADD_boot_10000.csv", row.names = FALSE)
} else {
  ## If this has already been run, load previous data
  heritability.ADD.boot <- fread("data/Heritability_Bootstrap/heritability_ADD_boot_10000.csv")
}
rm(heritability.ADD.boot.fish)


# STATISTICAL ANALYSIS - GENERATE OBSERVED HERITABILITY ---------------------------------------

## Embryo Survival
heritability.survival.obs <- do.call(rbind, lapply(unique(hatch.survival$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.survival %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, hatch)
    
    obs.survival <- observGlmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "hatch",
                                 block = "block", fam_link = binomial(logit))
    
    obs.survival.df <- data.frame(population = pop, 
                                  temperature = temp,
                                  block = obs.survival$random[4,2],
                                  residual = obs.survival$other[1,2],
                                  additive.obs = obs.survival$calculation[1,2],
                                  nonadd = obs.survival$calculation[2,2]) %>% 
      mutate(pheno.obs = additive.obs + nonadd + residual + block,
             h2.obs = additive.obs / pheno.obs) %>% 
      select(population, temperature, additive.obs, pheno.obs, h2.obs)
  }))
}))

## DPF
heritability.dpf.obs <- do.call(rbind, lapply(unique(hatch.dpf$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.dpf %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, dpf)
    
    obs.dpf <- observLmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "dpf", block = "block")
    
    obs.dpf.df <- data.frame(population = pop, 
                             temperature = temp,
                             block = obs.dpf$random[4,2],
                             residual = obs.dpf$other[1,2],
                             additive.obs = obs.dpf$calculation[1,2],
                             nonadd = obs.dpf$calculation[2,2]) %>% 
      mutate(pheno.obs = additive.obs + nonadd + residual + block,
             h2.obs = additive.obs / pheno.obs) %>% 
      select(population, temperature, additive.obs, pheno.obs, h2.obs)
  }))
}))

## ADD
heritability.ADD.obs <- do.call(rbind, lapply(unique(hatch.add$temperature), function(temp) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.add %>% filter(temperature == temp)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  do.call(rbind, lapply(unique(data.temp$population), function(pop) {
    ## Filter to a single group
    data.temp.group <- data.temp %>% filter(population == pop) %>% 
      select(family, dam, sire, block, ADD)
    
    obs.ADD <- observLmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "ADD", block = "block")
    
    obs.ADD.df <- data.frame(population = pop, 
                             temperature = temp,
                             block = obs.ADD$random[4,2],
                             residual = obs.ADD$other[1,2],
                             additive.obs = obs.ADD$calculation[1,2],
                             nonadd = obs.ADD$calculation[2,2]) %>% 
      mutate(pheno.obs = additive.obs + nonadd + residual + block,
             h2.obs = additive.obs / pheno.obs) %>% 
      select(population, temperature, additive.obs, pheno.obs, h2.obs)
  }))
}))


# CALCULATE THE BIAS-CORRECTED MEAN AND SE FROM BOOTSTRAPPED DISTRIBUTIONS  -------------------

## Embryo Survival
heritability.survival.summary <- heritability.survival.boot %>% 
  mutate(pheno = additive + nonadd + residual + block,
         h2.boot = additive / pheno) %>% 
  left_join(heritability.survival.obs) %>% 
  group_by(population, temperature) %>% 
  summarize(h2.obs = mean(h2.obs),
            h2.obs.bias = h2.obs - (mean(h2.boot)-h2.obs),
            h2.se = sd(h2.boot),
            var.add.obs = mean(additive.obs),
            var.add.obs.bias = var.add.obs - (mean(additive.obs)-var.add.obs),
            var.add.se = sd(additive),
            var.pheno.obs = mean(pheno.obs),
            var.pheno.obs.bias = var.pheno.obs - (mean(pheno.obs)-var.pheno.obs),
            var.pheno.se = sd(pheno)) %>% 
  mutate(trait = "survival",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2)

## DPF
heritability.dpf.summary <- heritability.dpf.boot %>% 
  mutate(pheno = additive + nonadd + residual + block,
         h2.boot = additive / pheno) %>% 
  left_join(heritability.dpf.obs) %>% 
  group_by(population, temperature) %>% 
  summarize(h2.obs = mean(h2.obs),
            h2.obs.bias = h2.obs - (mean(h2.boot)-h2.obs),
            h2.se = sd(h2.boot),
            var.add.obs = mean(additive.obs),
            var.add.obs.bias = var.add.obs - (mean(additive.obs)-var.add.obs),
            var.add.se = sd(additive),
            var.pheno.obs = mean(pheno.obs),
            var.pheno.obs.bias = var.pheno.obs - (mean(pheno.obs)-var.pheno.obs),
            var.pheno.se = sd(pheno)) %>% 
  mutate(trait = "dpf",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2)

## ADD
heritability.ADD.summary <- heritability.ADD.boot %>% 
  mutate(pheno = additive + nonadd + residual + block,
         h2.boot = additive / pheno) %>% 
  left_join(heritability.ADD.obs) %>% 
  group_by(population, temperature) %>% 
  summarize(h2.obs = mean(h2.obs),
            h2.obs.bias = h2.obs - (mean(h2.boot)-h2.obs),
            h2.se = sd(h2.boot),
            var.add.obs = mean(additive.obs),
            var.add.obs.bias = var.add.obs - (mean(additive.obs)-var.add.obs),
            var.add.se = sd(additive),
            var.pheno.obs = mean(pheno.obs),
            var.pheno.obs.bias = var.pheno.obs - (mean(pheno.obs)-var.pheno.obs),
            var.pheno.se = sd(pheno)) %>% 
  mutate(trait = "ADD",
         h2.obs.bias = ifelse(h2.obs.bias < 0, 0, h2.obs.bias)) %>% 
  ## Round numeric columns
  mutate_if(is.numeric, round, 2)


# CALCULATE SAMPLE SIZES ----------------------------------------------------------------------

heritability.survival.n <- hatch.survival %>% group_by(population, temperature) %>% 
  summarise(n = n()) %>% mutate(trait = "survival")
heritability.dpf.n <- hatch.dpf %>% group_by(population, temperature) %>% 
  summarise(n = n()) %>% mutate(trait = "dpf")
heritability.ADD.n <- hatch.add %>% group_by(population, temperature) %>% 
  summarise(n = n()) %>% mutate(trait = "ADD")
heritability.n <- bind_rows(heritability.survival.n, heritability.dpf.n, heritability.ADD.n)


# COMBINE ALL TRAITS --------------------------------------------------------------------------

heritability.all <- bind_rows(heritability.survival.summary, heritability.ADD.summary, heritability.dpf.summary) %>% 
  left_join(heritability.n) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")),
         population = factor(population, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco"))) 

#%>% 
#  filter(population != "LK-Whitefish")
#write.csv(heritability.all, "data/Heritability_Bootstrap/heritability_all.csv", row.names = FALSE)


# VISUALIZATION - HERITABILITY --------------------------------------------

ggplot(heritability.all, aes(x = temperature, y = (h2.obs.bias * 100), group = population, fill = population)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  #geom_text(aes(x = temperature, y = (h2.upper.CI * 100) + 2, label = n), size = 5, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ifelse((h2.obs.bias - h2.se) * 100 < 0, 0, (h2.obs.bias - h2.se) * 100), ymax = (h2.obs.bias + h2.se) * 100), 
                position = position_dodge(0.9),
                size = 1.0, width = 0.25, linetype = "solid", show.legend = FALSE) +
  scale_fill_grey("combine", start = 0.2, end = 0.9,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 10), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (%)") +
  theme_bw() + 
  theme(axis.title.x = element_text(color = "Black", size = 24, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 24, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 18),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) + 
  facet_wrap(~trait, nrow = 1)

ggsave("figures/embryo/2020-Heritability-fullfact-10000.png", width = 20, height = 12, dpi = 300)


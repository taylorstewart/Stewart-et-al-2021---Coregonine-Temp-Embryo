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

# Define the number of bootstrap samples
boot_n <- 5000

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.survival.boot <- do.call(rbind, lapply(1:boot_n, function(i) {
  
  ## display current loop location
  print(i)

  do.call(rbind, lapply(unique(hatch.survival$temperature), function(j) {
    ## Filter to only a single temperature treatment
    data.temp <- hatch.survival %>% filter(temperature == j)
    
    ## Apply a nested loop to run each group (lake:species) within each temperature treatment
    temp <- do.call(rbind, lapply(unique(data.temp$group), function(k) {
      ## Filter to a single group
      data.temp.group <- filter(data.temp, group == k)
      
      ## Create a bootstrapped data set from each temperature and group
      bootstrap.data <- do.call(rbind, lapply(1:nrow(data.temp.group), function(na) {
        ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
        hatch <- sample(data.temp.group$hatch, replace = T, size = 1)
        male <- sample(data.temp.group$male, replace = T, size = 1)
        female <- sample(data.temp.group$female, replace = T, size = 1)
        block <- sample(data.temp.group$block, replace = T, size = 1)
        plate <- sample(data.temp.group$plate, replace = T, size = 1)
        
        ## combine each variable into a data frame
        single.random.obs.data <- data.frame(hatch, male, female, block, plate)
      }))
      
      ## Fit random-effect model
      herit.surv.model <- glmer(hatch ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                                family = binomial, 
                                data = bootstrap.data, 
                                control = glmerControl(optimizer = "Nelder_Mead" ,optCtrl = list(maxfun = 100000)),
                                nAGQ = 1)
      
      ## Calculate genetic variance 
      var.additive <- VarCorr(herit.surv.model)$male[1] * 4
      
      ## Calculate environmental sources of variance
      var.pheno <- var.additive +
        (VarCorr(herit.surv.model)$female[1]) +
        (VarCorr(herit.surv.model)$"male:female"[1]) +
        (VarCorr(herit.surv.model)$plate [1]) +
        (VarCorr(herit.surv.model)$block[1])
      
      ## Calculate heritability
      heritability <- var.additive / var.pheno
      
      herit.surv.df <- data.frame(group = k, temperature = j, var.add = round(var.additive, 2), 
                                  var.pheno = round(var.pheno, 2), herit = round(heritability, 2))
    }))
  }))
}))


## Calculate the mean and SE from bootstrapped distribution
heritability.survival.summary <- heritability.survival.boot %>% group_by(temperature, group) %>% 
  filter(!is.na(herit)) %>% 
  summarize(mean.herit = mean(herit),
            mean.var.add = mean(var.add),
            mean.var.pheno = mean(var.pheno),
            sd.herit = sd(herit),
            sd.var.add = sd(var.add),
            sd.var.pheno = sd(var.pheno),
            se.herit = sd.herit/sqrt(n()),
            se.var.add = sd.var.add/sqrt(n()),
            se.var.pheno = sd.var.pheno/sqrt(n()),
            herit.upper.CI = quantile(herit, 0.975),
            herit.lower.CI = quantile(herit, 0.025),
            herit.upper.CI = quantile(var.add, 0.975),
            herit.lower.CI = quantile(var.add, 0.025),
            herit.upper.CI = quantile(var.pheno, 0.975),
            herit.lower.CI = quantile(var.pheno, 0.025)) %>% 
  mutate(trait = "surv")


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - HERITABILITY -----------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf))

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.dpf.boot <- do.call(rbind, lapply(1:boot_n, function(i) {
  
  ## display current loop location
  print(i)
  
  do.call(rbind, lapply(unique(hatch.dpf$temperature), function(j) {
    ## Filter to only a single temperature treatment
    data.temp <- hatch.dpf %>% filter(temperature == j)
    
    ## Apply a nested loop to run each group (lake:species) within each temperature treatment
    do.call(rbind, lapply(unique(data.temp$group), function(k) {
      ## Filter to a single group
      data.temp.group <- filter(data.temp, group == k)
      
      ## Create a bootstrapped data set from each temperature and group
      bootstrap.data <- do.call(rbind, lapply(1:nrow(data.temp.group), function(na) {
        ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
        dpf <- sample(data.temp.group$dpf, replace = T, size = 1)
        male <- sample(data.temp.group$male, replace = T, size = 1)
        female <- sample(data.temp.group$female, replace = T, size = 1)
        block <- sample(data.temp.group$block, replace = T, size = 1)
        plate <- sample(data.temp.group$plate, replace = T, size = 1)
        
        ## combine each variable into a data frame
        single.random.obs.data <- data.frame(dpf, male, female, block, plate)
      }))
      
      ## Fit random-effect model
      herit.dpf.model <- lmer(dpf ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                              data = bootstrap.data,
                              REML = TRUE,
                              control = lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
      
      ## Calculate genetic variance 
      var.additive <- VarCorr(herit.dpf.model)$male[1] * 4
      
      ## Calculate environmental sources of variance
      var.pheno <- var.additive +
        (VarCorr(herit.dpf.model)$female[1]) +
        (VarCorr(herit.dpf.model)$"male:female"[1]) +
        (VarCorr(herit.dpf.model)$plate [1]) +
        (VarCorr(herit.dpf.model)$block[1])
      
      ## Calculate heritability
      heritability <- var.additive / var.pheno
      
      herit.dpf.df <- data.frame(group = k, temperature = j, var.add = round(var.additive, 2), 
                                 var.pheno = round(var.pheno, 2), herit = round(heritability, 2))
    }))
  }))
}))


## Calculate the mean and SE from bootstrapped distribution
heritability.dpf.summary <- heritability.dpf.boot %>% group_by(temperature, group) %>% 
  filter(!is.na(herit)) %>% 
  summarize(mean.herit = mean(herit),
            mean.var.add = mean(var.add),
            mean.var.pheno = mean(var.pheno),
            sd.herit = sd(herit),
            sd.var.add = sd(var.add),
            sd.var.pheno = sd(var.pheno),
            se.herit = sd.herit/sqrt(n()),
            se.var.add = sd.var.add/sqrt(n()),
            se.var.pheno = sd.var.pheno/sqrt(n()),
            herit.upper.CI = quantile(herit, 0.975),
            herit.lower.CI = quantile(herit, 0.025),
            herit.upper.CI = quantile(var.add, 0.975),
            herit.lower.CI = quantile(var.add, 0.025),
            herit.upper.CI = quantile(var.pheno, 0.975),
            herit.lower.CI = quantile(var.pheno, 0.025)) %>% 
  mutate(trait = "dpf")


# STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - HERITABILITY -----------

# filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD))

i <- 1
j <- "4.5°C"
k <- "LO-Cisco"

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.add.boot <- do.call(rbind, lapply(unique(hatch.ADD$temperature), function(j) {
    
  ## Filter to only a single temperature treatment
    data.temp <- hatch.ADD %>% filter(temperature == j)
    
    ## Apply a nested loop to run each group (lake:species) within each temperature treatment
    do.call(rbind, lapply(unique(data.temp$population), function(k) {
      ## Filter to a single group
      data.temp.group <- filter(data.temp, population == k) %>% 
        select(dam, sire, block, tray, ADD.raw = ADD)
      
      ## Create a bootstrapped data set from each temperature and group
      bootstrap.data <- do.call(cbind, lapply(1:1000, function(l) {
        ## create a vector of randomly generated data (randomly sample each variable and repeat by nrow(data.temp.group))
        ADD.boot <- sample(data.temp.group$ADD.raw, replace = T, size = nrow(data.temp.group))
        
        boot <- bind_cols(data.temp.group, ADD.boot) %>% select(-ADD.raw)
        colnames(boot) <- c(paste0("dam", l), paste0("sire", l), paste0("block", l), 
                            paste0("tray", l), paste0("ADD", l))
        boot
      }))
      
      ## Fit standard random-effect model
      herit.add.model <- observLmer2(observ = data.temp.group, dam = "dam", sire = "sire", response = "ADD.raw",
                                     block = "block", position = "tray")
      
      ## Create a vector of variances for bias correction
      var.tray <- herit.add.model$random[1,2]
      var.block <- herit.add.model$random[5,2]
      var.total <- herit.add.model$other[2,2]
      var.additive <- herit.add.model$calculation[1,2]
      var.non.additive <- herit.add.model$calculation[2,2]
      var.maternal <- herit.add.model$calculation[3,2]
      
      var.vector <- c(var.additive, var.non.additive, var.maternal, var.total, var.tray, var.block)
      
      ## Calculate variance components from bootstrapped sample
      bootstrap.data.lmer2 <- resampLmer2(resamp = bootstrap.data, dam = "dam", sire = "sire", response = "ADD",
                                          position = "tray", block = "block", start = 1, end = 1000)
      
      ## Calculate confidence intervals for variance components
      bootstrap.data.lmer2.ci <- ciMANA2(comp = bootstrap.data.lmer2, bias = var.vector,
                                         position = "tray", block = "block", trait = "ADD")
      
      variances.add.df <- data.frame(population = k, temperature = j, bootstrap.data.lmer2.ci)
  }))
}))

## Calculate the mean and SE from bootstrapped distribution
heritability.ADD.summary <- heritability.add.boot %>% group_by(temperature, group) %>% 
  filter(!is.na(herit)) %>% 
  summarize(mean.herit = mean(herit),
            mean.var.add = mean(var.add),
            mean.var.pheno = mean(var.pheno),
            sd.herit = sd(herit),
            sd.var.add = sd(var.add),
            sd.var.pheno = sd(var.pheno),
            se.herit = sd.herit/sqrt(n()),
            se.var.add = sd.var.add/sqrt(n()),
            se.var.pheno = sd.var.pheno/sqrt(n()),
            herit.upper.CI = quantile(herit, 0.975),
            herit.lower.CI = quantile(herit, 0.025),
            herit.upper.CI = quantile(var.add, 0.975),
            herit.lower.CI = quantile(var.add, 0.025),
            herit.upper.CI = quantile(var.pheno, 0.975),
            herit.lower.CI = quantile(var.pheno, 0.025)) %>% 
  mutate(trait = "ADD")


# VISUALIZATION - HERITABILITY --------------------------------------------

heritability.all <- bind_rows(heritability.survival.summary, heritability.ADD.summary, heritability.dpf.summary) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("surv", "dpf", "ADD"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")),
         label = ifelse(mean.herit == 0, 0, NA))

ggplot(heritability.ADD.summary, aes(x = temperature, y = (mean.herit * 100), group = group, fill = group)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  #geom_text(aes(x = temperature, y = 0.75, label = label), position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = (mean.herit - se.herit) * 100, ymax = (mean.herit + se.herit) * 100), 
                position = position_dodge(0.9),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_fill_grey("combine", start = 0.2, end = 0.9,
                  labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_y_continuous(limits = c(-0.5, 75), breaks = seq(0, 75, 25), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = "Incubation Temperature (°C)", y = "Narrow-sense Heritability (%)") +
  theme_bw() + 
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 15),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) #+ 
  #facet_wrap(~trait, nrow = 1)

ggsave("figures/embryo/2020-Heritability.png", width = 20, height = 12, dpi = 300)


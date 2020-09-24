# Clear the environment first ---------------------------------------------

rm(list = ls(all.names = TRUE))

# Load packages -----------------------------------------------------------

library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)
library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")

# Load incubation temperature data ----------------------------------------

ADD.2020 <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())

# Load hatching data ------------------------------------------------------

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
                       labels = c("M\n1", "M\n2", "M\n3", "M\n4", "M\n5", "M\n6", "M\n7", "M\n8", "M\n9", "M\n10", "M\n11", "M\n12", "M\n13", "M\n14", "M\n15", "M\n16")),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Clean up environment
rm(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus, ADD.2020)


# Summarize Data - Population Level ---------------------------------------

## Embryo Survival
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.survival = (mean(hatch)*100),
            sd.survival = (sd(hatch)*100),
            se.survival = sd.survival / sqrt(n()))

## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(hatch == 1, premature != 1, !is.na(dpf)) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf / sqrt(n()))

## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1, premature != 1) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.ADD = mean(ADD),
            sd.ADD = sd(ADD),
            se.ADD = sd.ADD / sqrt(n()))


# Summarize Data - Family Level -------------------------------------------

## Embryo Survival
hatch.survival.family.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, species, temperature, group, male, female, block) %>% 
  summarize(mean.survival = (mean(hatch)*100),
            sd.survival = (sd(hatch)*100),
            se.survival = sd.survival / sqrt(n()))

## Days Post Fertilization
hatch.dpf.family.summary <- hatch %>% filter(hatch == 1, premature != 1, !is.na(dpf)) %>% 
  group_by(population, species, temperature, group, male, female, block) %>% 
  summarize(mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf / sqrt(n()))

## Accumulated Degree-Days
hatch.ADD.family.summary <- hatch %>% filter(!is.na(ADD), hatch == 1, premature != 1) %>% 
  group_by(population, species, temperature, group, male, female, block) %>% 
  summarize(mean.ADD = mean(ADD),
            sd.ADD = sd(ADD),
            se.ADD = sd.ADD / sqrt(n()))


# Find Highest and Lowest SE for Each Male --------------------------------

## Embryo Survival
hatch.survival.family.summary <- do.call(rbind, lapply(unique(hatch.survival.family.summary$male), function(i) {
  male <- hatch.survival.family.summary %>% filter(male == i)
  
  male.se <- male %>% group_by(population, species, temperature, group, block) %>% 
    mutate(se.high = ifelse(mean.survival == max(mean.survival), se.survival, NA),
           se.low = ifelse(mean.survival == min(mean.survival), se.survival, NA),
           se.high = ifelse(!is.na(se.low) & is.na(se.high), 0, se.high),
           se.low = ifelse(!is.na(se.high) & is.na(se.low), 0, se.low))
}))

## Days Post Fertilization
hatch.dpf.family.summary <- do.call(rbind, lapply(unique(hatch.dpf.family.summary$male), function(i) {
  male <- hatch.dpf.family.summary %>% filter(male == i)
  
  male.se <- male %>% group_by(population, species, temperature, group, block) %>% 
    mutate(se.high = ifelse(mean.dpf == max(mean.dpf), se.dpf, NA),
           se.low = ifelse(mean.dpf == min(mean.dpf), se.dpf, NA),
           se.high = ifelse(!is.na(se.low) & is.na(se.high), 0, se.high),
           se.low = ifelse(!is.na(se.high) & is.na(se.low), 0, se.low))
}))

## Accumulated Degree-Days
hatch.ADD.family.summary <- do.call(rbind, lapply(unique(hatch.ADD.family.summary$male), function(i) {
  male <- hatch.ADD.family.summary %>% filter(male == i)
  
  male.se <- male %>% group_by(population, species, temperature, group, block) %>% 
    mutate(se.high = ifelse(mean.ADD == max(mean.ADD), se.ADD, NA),
           se.low = ifelse(mean.ADD == min(mean.ADD), se.ADD, NA),
           se.high = ifelse(!is.na(se.low) & is.na(se.high), 0, se.high),
           se.low = ifelse(!is.na(se.high) & is.na(se.low), 0, se.low))
}))


# Statistical Analysis - Survival - GLM -----------------------------------

# filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

# create linear mixed model
hatch.survival.glm <- glmer(hatch ~ temperature + group + temperature * group +        # fixed
                            (1|male) + (1|female),                                     # random
                            family = binomial("logit"),
                            data = hatch.survival,
                            control = glmerControl(optimizer = "bobyqa"))  

# to select all model based on AICc
hatch.survival.glm.AIC <- dredge(hatch.survival.glm)

# select best model based on AICc
hatch.survival.glm.AIC.best <- get.models(hatch.survival.glm.AIC, 1)[[1]]

# ANOVA
Anova(hatch.survival.glm.AIC.best)
Anova(hatch.survival.glm.AIC.best, type = "III")

# post-hoc test
hatch.survival.glm.emm <- emmeans(hatch.survival.glm.AIC.best, ~ temperature * group)
(hatch.survival.glm.emm.pair <- pairs(hatch.survival.glm.emm, simple = "each", adjust = "fdr"))


# Statistical Analysis - Survival - Heritability --------------------------

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.survival <- do.call(rbind, lapply(unique(hatch.survival$temperature), function(i) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.survival %>% filter(temperature == i)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  temp <- do.call(rbind, lapply(unique(data.temp$group), function(j) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, group == j)
    
    ## Fit random-effect model
    herit.surv.model <- lmer(hatch ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                             data = data.temp.group, REML = F)
    
    ## Calculate genetic variance 
    var.sire <- VarCorr(herit.surv.model)$male[1]
    
    ## Calculate environmental sources of variance
    var.resid <- (VarCorr(herit.surv.model)$male[1])+
      (VarCorr(herit.surv.model)$female[1])+
      (VarCorr(herit.surv.model)$"male:female"[1])+
      (VarCorr(herit.surv.model)$plate [1])+
      (VarCorr(herit.surv.model)$block[1])
    
    ## Calculate heritability
    heritability <- 4 * var.sire / (var.sire + var.resid)
    
    herit.surv.df <- data.frame(group = j, temperature = i, var.sire = round(var.sire, 2), var.resid = round(var.resid, 2), herit = round(heritability, 2))
  }))
}))


# Statistical Analysis - Incubation Period (ADD) - GLM --------------------

# filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD))

# create linear mixed model
hatch.ADD.glm <- lmer(ADD ~ temperature + group + temperature * group +       # Fixed
                      (1|male) + (1|female),                                  # Random
                      data = hatch.ADD, 
                      REML = FALSE)

# to select all model based on AICc
hatch.ADD.glm.AIC <- dredge(hatch.ADD.glm)

# select best model based on AICc
hatch.ADD.glm.AIC.best <- get.models(hatch.ADD.glm.AIC, 1)[[1]]

# ANOVA
Anova(hatch.ADD.glm.AIC.best)
Anova(hatch.ADD.glm.AIC.best, type = "III")

# post-hoc test
hatch.ADD.glm.emm <- emmeans(hatch.ADD.glm.AIC.best, ~ temperature * group)
pairs(hatch.ADD.glm.emm, simple = list("group", c("temperature")), adjust = "fdr") 


# Statistical Analysis - Incubation Period (ADD) - Heritability -----------

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.ADD <- do.call(rbind, lapply(unique(hatch.ADD$temperature), function(i) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.ADD %>% filter(temperature == i)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  temp <- do.call(rbind, lapply(unique(data.temp$group), function(j) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, group == j)
    
    ## Fit random-effect model
    herit.ADD.model <- lmer(ADD ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                            data = data.temp.group, REML = F)
    
    ## Calculate genetic variance 
    var.sire <- VarCorr(herit.ADD.model)$male[1]
    
    ## Calculate environmental sources of variance
    var.resid <- (VarCorr(herit.ADD.model)$male[1])+
      (VarCorr(herit.ADD.model)$female[1])+
      (VarCorr(herit.ADD.model)$"male:female"[1])+
      (VarCorr(herit.ADD.model)$plate [1])+
      (VarCorr(herit.ADD.model)$block[1])
    
    ## Calculate heritability
    heritability <- 4 * var.sire / (var.sire + var.resid)
    
    herit.ADD.df <- data.frame(group = j, temperature = i, var.sire = round(var.sire, 2), var.resid = round(var.resid, 2), herit = round(heritability, 2))
  }))
}))


# Statistical Analysis - Incubation Period (DPF) - GLM --------------------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf))

# create linear mixed model
hatch.dpf.glm <- lmer(dpf ~ temperature + group + temperature * group +       # Fixed
                      (1|male) + (1|female),                                  # Random
                      data = hatch.dpf, 
                      REML = FALSE)

# to select all model based on AICc
hatch.dpf.glm.AIC <- dredge(hatch.dpf.glm)

# select best model based on AICc
hatch.dpf.glm.AIC.best <- get.models(hatch.dpf.glm.AIC, 1)[[1]]

# ANOVA
Anova(hatch.dpf.glm.AIC.best)
Anova(hatch.dpf.glm.AIC.best, type = "III")

# post-hoc test
hatch.dpf.glm.emm <- emmeans(hatch.dpf.glm.AIC.best, ~ temperature * group)
pairs(hatch.dpf.glm.emm, simple = list("group", c("temperature")), adjust = "fdr") 


# Statistical Analysis - Incubation Period (DPF) - Heritability -----------

## Run nested loop to run narrow-sense heritability calculations for each temperature and population
heritability.dpf <- do.call(rbind, lapply(unique(hatch.dpf$temperature), function(i) {
  ## Filter to only a single temperature treatment
  data.temp <- hatch.dpf %>% filter(temperature == i)
  
  ## Apply a nested loop to run each group (lake:species) within each temperature treatment
  temp <- do.call(rbind, lapply(unique(data.temp$group), function(j) {
    ## Filter to a single group
    data.temp.group <- filter(data.temp, group == j)
    
    ## Fit random-effect model
    herit.dpf.model <- lmer(dpf ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                            data = data.temp.group, REML = F)
    
    ## Calculate genetic variance 
    var.sire <- VarCorr(herit.dpf.model)$male[1]
    
    ## Calculate environmental sources of variance
    var.resid <- (VarCorr(herit.dpf.model)$male[1])+
      (VarCorr(herit.dpf.model)$female[1])+
      (VarCorr(herit.dpf.model)$"male:female"[1])+
      (VarCorr(herit.dpf.model)$plate [1])+
      (VarCorr(herit.dpf.model)$block[1])
    
    ## Calculate heritability
    heritability <- 4 * var.sire / (var.sire + var.resid)
    
    herit.dpf.df <- data.frame(group = j, temperature = i, var.sire = round(var.sire, 2), var.resid = round(var.resid, 2), herit = round(heritability, 2))
  }))
}))


# Visualizations - Family Level -------------------------------------------

## Embryo Survival
ggplot(hatch.survival.family, aes(x = male, y = mean.survival, group = female, shape = female)) + 
  geom_point(size = 4.0) +
  geom_errorbar(aes(ymin = mean.survival-se.low, ymax = mean.survival+se.high), 
                size = 0.8, width = 0.3) +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), linetype = "dashed") +
  scale_x_discrete(expand = c(0, 0.3)) +
  scale_y_continuous(limits = c(-5, 110), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0)) +
  labs(y = "Embryo Survival (% ± SE)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        panel.spacing = unit(5, "mm"),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  facet_grid(group ~ temperature)

ggsave("figures/embryo/2020-Survival-Family.png", width = 20, height = 12, dpi = 300)

## Days Post Fertilization
ggplot(hatch.dpf.family, aes(x = male, y = mean.dpf, group = female, shape = female)) + 
  geom_point(size = 4.0) +
  geom_errorbar(aes(ymin = mean.dpf-se.low, ymax = mean.dpf+se.high), 
                size = 0.8, width = 0.3) +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), linetype = "dashed") +
  scale_x_discrete(expand = c(0, 0.35)) +
  scale_y_continuous(limits = c(25, 230), breaks = seq(25, 225, 50), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0)) +
  labs(y = "Incubation Period (No. Days ± SE)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        panel.spacing = unit(5, "mm"),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  facet_grid(group ~ temperature)

ggsave("figures/embryo/2020-DPF-Family.png", width = 20, height = 12, dpi = 300)

## Accumulated Degree-Days
ggplot(hatch.ADD.family, aes(x = male, y = mean.ADD, group = female, shape = female)) + 
  geom_point(size = 4.0) +
  geom_errorbar(aes(ymin = mean.ADD-se.low, ymax = mean.ADD+se.high),
                size = 0.8, width = 0.3) +
  geom_vline(xintercept = c(4.5, 8.5, 12.5), linetype = "dashed") +
  scale_x_discrete(expand = c(0, 0.35)) +
  scale_y_continuous(limits = c(150, 940), breaks = seq(200, 900, 100), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0)) +
  labs(y = "Incubation Period (ADD °C ± SE)") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        panel.spacing = unit(5, "mm"),
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) +
  facet_grid(group ~ temperature)

ggsave("figures/embryo/2020-ADD-Family.png", width = 20, height = 12, dpi = 300)


# Visualizations - Population Level ---------------------------------------

## Embryo Survival
ggplot(hatch.survival.summary, aes(x = temperature, y = mean.survival, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.25, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = mean.survival-se.survival, ymax = mean.survival+se.survival), 
                position = position_dodge(0.15),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(10, 101), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Embryo Survival (% ± SE)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        #legend.position = c(0.8, 0.86),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-Survival-BW.png", width = 8.5, height = 6, dpi = 300)

## Days Post Fertilization
ggplot(hatch.dpf.summary, aes(x = temperature, y = mean.dpf, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.25, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = mean.dpf-se.dpf, ymax = mean.dpf+se.dpf), 
                position = position_dodge(0.15),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(45, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (No. Days ± SE)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        #legend.position = c(0.8, 0.86),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-DPF-BW.png", width = 8.5, height = 6, dpi = 300)

## Accumulated Degree-Days
ggplot(hatch.ADD.summary, aes(x = temperature, y = mean.ADD, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.25, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = mean.ADD-sd.ADD, ymax = mean.ADD+sd.ADD), 
                position = position_dodge(0.15),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(250, 950), breaks = seq(250, 950, 100), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (ADD °C ± SE)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        #legend.position = c(0.17, 0.88),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-ADD-BW.png", width = 8.5, height = 6, dpi = 300)



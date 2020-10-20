#### CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(buildmer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

emm_options(pbkrtest.limit = 14000)

#### LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

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
  dplyr::select(population, latitude, species, male, female, block, temperature, eye, hatch, dpf, ADD)

hatch.Finland <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, latitude, species, male, female, block, temperature, eye, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA, hatch.Finland) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 4.5, 7, 9)),
                              #labels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         family = factor(paste0(female, male)),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Clean up environment
rm(hatch.USA, hatch.Finland, ADD)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival <- hatch %>% filter(eye != 0)

## filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1)

## filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD), hatch == 1)


#### STATISTICAL ANALYSIS - SURVIVAL - GLM -----------------------------------

## Find best fit model - be patient!
hatch.survival.model <- buildmer(hatch ~ temperature + group + temperature:group + (1|family) + (1|male) + (1|female) + (1|block), 
                                 direction = 'forward', data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.glm <- formula(hatch.survival.model@model))
## hatch ~ 1 + group + temperature + group:temperature + (1|family) + (1|female)

## Create generalized linear mixed models
# fit best (full) model
hatch.survival.glm.best <- lme4::glmer(hatch.survival.glm, data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Create generalized linear mixed models with fixed-effects removed for LRTs
# fit best model with temperature main effect removed
hatch.survival.glm.best.temp <- lme4::glmer(hatch ~ 1 + group + (1|family) + (1|female),
                                            data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))
# fit best model with population main effect removed
hatch.survival.glm.best.pop <- lme4::glmer(hatch ~ 1 + temperature + (1|family) + (1|female),
                                          data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Create generalized linear mixed models with random-effects removed for LRTs
# fit best model with family random-effect removed
hatch.survival.glm.best.family <- lme4::glmer(hatch ~ 1 + group + temperature + group:temperature + (1|female),
                                              data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))
# fit best model with female random-effect removed
hatch.survival.glm.best.female <- lme4::glmer(hatch ~ 1 + group + temperature + group:temperature + (1|family),
                                            data = hatch.survival, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Calculate LRT for both temperature and population fixed-effects and interaction
# temperature:population (significant; main effects irrelevant)
drop1(hatch.survival.glm.best, test = "Chisq")
# temperature
anova(hatch.survival.glm.best.temp, hatch.survival.glm.best, test = "Chisq")
# population
anova(hatch.survival.glm.best, hatch.survival.glm.best.pop, test = "Chisq")

## Calculate LRT for female, male, and plate random-effects
# family
anova(hatch.survival.glm.best.family, hatch.survival.glm.best, test = "Chisq")
# female
anova(hatch.survival.glm.best.female, hatch.survival.glm.best, test = "Chisq")

## Calculate estimated marginal means - be very patient!
if(file.exists("data/emmeans/hatch_survival_glm_emm.csv") == FALSE) {
  hatch.survival.glm.emm <- emmeans(hatch.survival.glm.best, ~ temperature*group)

  ## Pairwise, cld, confidence intervals
  pairs(hatch.survival.glm.emm, simple = list("temperature", "group"), adjust = "fdr", type = "response") 
  hatch.survival.glm.emm.confint <- multcomp::cld(hatch.survival.glm.emm, type = "response", adjust = "fdr",
                                                  sort = F, alpha = 0.05, Letters = LETTERS) %>% 
    mutate(.group = gsub("[[:space:]]", "", .group))
  
  ## Save output to prevent having to re-run time consuming models
  write.csv(hatch.survival.glm.emm.confint, "data/emmeans/hatch_survival_glm_emm.csv", row.names = FALSE)
} else {
  hatch.survival.glm.emm.confint <- read.csv("data/emmeans/hatch_survival_glm_emm.csv", header = TRUE) %>% 
    mutate(temperature = factor(temperature, ordered = TRUE, levels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
           group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))
}

#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - GLM --------------------

## Find best fit model - be patient!
hatch.dpf.model <- buildmer(dpf ~ temperature + group + temperature:group + (1|male) + (1|female) + (1|family) + (1|block), 
                            direction = 'forward', data = hatch.dpf, REML = TRUE)
( hatch.dpf.glm <- formula(hatch.dpf.model@model))
## hatch ~ 1 + group + temperature + group:temperature + (1|family) + (1|female) + (1|male)

## Create generalized linear mixed models
# fit best (full) model
hatch.dpf.glm.best <- lme4::lmer(hatch.dpf.glm, data = hatch.dpf)

## Create generalized linear mixed models with fixed-effects removed for LRTs
# fit best model with temperature main effect removed
hatch.dpf.glm.best.temp <- lme4::lmer(dpf ~ 1 + group + (1|family) + (1|female) + (1|male),
                                            data = hatch.dpf)
# fit best model with population main effect removed
hatch.dpf.glm.best.pop <- lme4::lmer(dpf ~ 1 + temperature + (1|family) + (1|female) + (1|male),
                                           data = hatch.dpf)

## Create generalized linear mixed models with random-effects removed for LRTs
# fit best model with family random-effect removed
hatch.dpf.glm.best.family <- lme4::lmer(dpf ~ 1 + group + temperature + group:temperature + (1|female) + (1|male),
                                       data = hatch.dpf)
# fit best model with female random-effect removed
hatch.dpf.glm.best.female <- lme4::lmer(dpf ~ 1 + group + temperature + group:temperature + (1|family) + (1|male),
                                              data = hatch.dpf)
# fit best model with male random-effect removed
hatch.dpf.glm.best.male <- lme4::lmer(dpf ~ 1 + group + temperature + group:temperature + (1|family) + (1|female),
                                            data = hatch.dpf)

## Calculate LRT for both temperature and population fixed-effects and interaction
# temperature:population (significant; main effects irrelevant)
drop1(hatch.dpf.glm.best, test = "Chisq")
# temperature
anova(hatch.dpf.glm.best.temp, hatch.dpf.glm.best, test = "Chisq")
# population
anova(hatch.dpf.glm.best, hatch.dpf.glm.best.pop, test = "Chisq")

## Calculate LRT for female, male, and plate random-effects
# family
anova(hatch.dpf.glm.best.family, hatch.dpf.glm.best, test = "Chisq")
# female
anova(hatch.dpf.glm.best.female, hatch.dpf.glm.best, test = "Chisq")
# male
anova(hatch.dpf.glm.best.male, hatch.dpf.glm.best, test = "Chisq")

## Calculate estimated marginal means - be very patient!
if(file.exists("data/emmeans/hatch_dpf_glm_emm.csv") == FALSE) {
  hatch.dpf.glm.emm <- emmeans(hatch.dpf.glm.best, ~ temperature*group)

  ## Pairwise, cld, confidence intervals
  pairs(hatch.dpf.glm.emm, simple = "temperature", adjust = "fdr", type = "response") 
  hatch.dpf.glm.emm.confint <- multcomp::cld(hatch.dpf.glm.emm, type = "response", adjust = "fdr", sort = F, alpha = 0.05, Letters = LETTERS) %>% 
    mutate(.group = gsub("[[:space:]]", "", .group))
  
  ## Save output to prevent having to re-run time consuming models
  write.csv(hatch.dpf.glm.emm.confint, "data/emmeans/hatch_dpf_glm_emm.csv", row.names = FALSE)
} else {
  hatch.dpf.glm.emm.confint <- read.csv("data/emmeans/hatch_dpf_glm_emm.csv", header = TRUE) %>% 
    mutate(temperature = factor(temperature, ordered = TRUE, levels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
           group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))
}

#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - GLM --------------------

## Find best fit model - be patient!
hatch.ADD.model <- buildmer(ADD ~ temperature + group + temperature:group + (1|family) + (1|male) + (1|female) + (1|block), 
                            direction = 'forward', data = hatch.ADD, REML = TRUE)
( hatch.ADD.glm <- formula(hatch.ADD.model@model))
## ADD ~ 1 + group + temperature + group:temperature + (1|family) + (1|female) + (1|male)

## Create generalized linear mixed models
# fit best (full) model
hatch.ADD.glm.best <- lme4::lmer(hatch.ADD.glm, data = hatch.ADD)

## Create generalized linear mixed models with fixed-effects removed for LRTs
# fit best model with temperature main effect removed
hatch.ADD.glm.best.temp <- lme4::lmer(ADD ~ 1 + group + (1|family) + (1|female) + (1|male),
                                      data = hatch.ADD)
# fit best model with population main effect removed
hatch.ADD.glm.best.pop <- lme4::lmer(ADD ~ 1 + temperature + (1|family) + (1|female) + (1|male),
                                     data = hatch.ADD)

## Create generalized linear mixed models with random-effects removed for LRTs
# fit best model with family random-effect removed
hatch.ADD.glm.best.family <- lme4::lmer(ADD ~ 1 + group + temperature + group:temperature + (1|female) + (1|male),
                                       data = hatch.ADD)
# fit best model with female random-effect removed
hatch.ADD.glm.best.female <- lme4::lmer(ADD ~ 1 + group + temperature + group:temperature + (1|family) + (1|male),
                                        data = hatch.ADD)
# fit best model with male random-effect removed
hatch.ADD.glm.best.male <- lme4::lmer(ADD ~ 1 + group + temperature + group:temperature + (1|family) + (1|female),
                                      data = hatch.ADD)

## Calculate LRT for both temperature and population fixed-effects and interaction
# temperature:population (significant; main effects irrelevant)
drop1(hatch.ADD.glm.best, test = "Chisq")
# temperature
anova(hatch.ADD.glm.best.temp, hatch.ADD.glm.best, test = "Chisq")
# population
anova(hatch.ADD.glm.best, hatch.ADD.glm.best.pop, test = "Chisq")

## Calculate LRT for female, male, and plate random-effects
# family
anova(hatch.ADD.glm.best.family, hatch.ADD.glm.best, test = "Chisq")
# female
anova(hatch.ADD.glm.best.female, hatch.ADD.glm.best, test = "Chisq")
# male
anova(hatch.ADD.glm.best.male, hatch.ADD.glm.best, test = "Chisq")

## Calculate estimated marginal means - be very patient!
if(file.exists("data/emmeans/hatch_ADD_glm_emm.csv") == FALSE) {
  hatch.ADD.glm.emm <- emmeans(hatch.ADD.glm.best, ~ temperature*group)
  
  ## Pairwise, cld, confidence intervals
  pairs(hatch.ADD.glm.emm, simple = "temperature", adjust = "fdr", type = "response") 
  hatch.ADD.glm.emm.confint <- multcomp::cld(hatch.ADD.glm.emm, type = "response", adjust = "fdr",
                                             sort = F, alpha = 0.05, Letters = LETTERS) %>% 
    mutate(.group = gsub("[[:space:]]", "", .group))
  
  ## Save output to prevent having to re-run time consuming models
  write.csv(hatch.ADD.glm.emm.confint, "data/emmeans/hatch_ADD_glm_emm.csv", row.names = FALSE)
} else {
  hatch.ADD.glm.emm.confint <- read.csv("data/emmeans/hatch_ADD_glm_emm.csv", header = TRUE) %>% 
    mutate(temperature = factor(temperature, ordered = TRUE, levels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
           group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))
}


#### CORRELATION WITH LATITUDE -------------------------------------------------------------------

## Create dataframe with latitudes to merge
latitude.df <- data.frame(group = c("LO-Cisco", "LS-Cisco", "LK-Whitefish", "LK-Vendace"),
                          latitude = c(44, 47, 62.5, 62.5))

## Embryo Survival
hatch.survival.corr <- hatch.survival.glm.emm.confint %>% filter(group != "LK-Whitefish") %>% 
  left_join(latitude.df) %>% 
  group_by(group, latitude) %>% 
  filter(prob == max(prob)) %>% ungroup() %>% 
  mutate(temperature = as.numeric(as.character(gsub("°C", "", temperature))),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "albula", temperature + 0.05, temperature),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "lavaretus", temperature - 0.05, temperature.dodge),
         scaled.survival = prob / max(prob))

lm.survival.temp <- lm(temperature ~ latitude, data = hatch.survival.corr)
summary(lm.survival.temp)
lm.survival.cgv <- lm(scaled.survival ~ latitude, data = hatch.survival.corr)
summary(lm.survival.cgv)

## Incubation Period (DPF)
hatch.dpf.corr <- hatch.dpf.glm.emm.confint %>% filter(group != "LK-Whitefish") %>% 
  left_join(latitude.df) %>% 
  group_by(group, latitude) %>% 
  filter(emmean == max(emmean)) %>% ungroup() %>% 
  mutate(temperature = as.numeric(as.character(gsub("°C", "", temperature))),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "albula", temperature + 0.05, temperature),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "lavaretus", temperature - 0.05, temperature.dodge),
         scaled.dpf = emmean / max(emmean))

lm.dpf.temp <- lm(temperature ~ latitude, data = hatch.dpf.corr)
summary(lm.dpf.temp)
lm.dpf.cgv <- lm(scaled.dpf ~ latitude, data = hatch.dpf.corr)
summary(lm.dpf.cgv)

## Incubation Period (ADD)
hatch.ADD.corr <- hatch.ADD.glm.emm.confint %>% filter(group != "LK-Whitefish") %>% 
  left_join(latitude.df) %>% 
  group_by(group, latitude) %>% 
  filter(emmean == max(emmean)) %>% ungroup() %>% 
  mutate(temperature = as.numeric(as.character(gsub("°C", "", temperature))),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "albula", temperature + 0.05, temperature),
         #temperature.dodge = ifelse(population == "konnevesi" & species == "lavaretus", temperature - 0.05, temperature.dodge),
         scaled.ADD = emmean / max(emmean))

lm.ADD.temp <- lm(temperature ~ latitude, data = hatch.ADD.corr)
summary(lm.dpf.temp)
lm.ADD.cgv <- lm(scaled.ADD ~ latitude, data = hatch.ADD.corr)
summary(lm.dpf.cgv)



# VISUALIZATIONS - CORRELATIONS ---------------------------------------------------------------

## Embryo Survival
plot.survival.corr.temp <- ggplot(hatch.survival.corr, aes(x = latitude, y = temperature, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = temperature), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1, 3), breaks = 2, expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Incubation Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

plot.survival.corr.cgv <- ggplot(hatch.survival.corr, aes(x = latitude, y = scaled.survival, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = scaled.survival), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.7, 1.02), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Relativized Embryo Survival") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Days Post Fertilization
plot.dpf.corr.temp <- ggplot(hatch.dpf.corr, aes(x = latitude, y = temperature, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = temperature), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1, 3), breaks = 2, expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Incubation Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

plot.dpf.corr.cgv <- ggplot(hatch.dpf.corr, aes(x = latitude, y = scaled.dpf, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = scaled.dpf), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.7, 1.02), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Relativized DPF") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Accumulated Degree-Days
plot.ADD.corr.temp <- ggplot(hatch.ADD.corr, aes(x = latitude, y = temperature, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = temperature), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(6, 8), breaks = 7, expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Incubation Temperature (°C)") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

plot.ADD.corr.cgv <- ggplot(hatch.ADD.corr, aes(x = latitude, y = scaled.ADD, fill = group, shape = group)) + 
  geom_smooth(aes(x = latitude, y = scaled.ADD), method = "lm", se = FALSE, inherit.aes = FALSE, color = "black") + 
  geom_point(size = 5, color = "black") +
  #annotate("text", x = 43.5, y = 0.575, hjust = 0, size = 6,
  #         label = "All Maximum Values at 7.0°C") +
  scale_x_continuous(limits = c(42.5, 63.5), breaks = seq(43, 63, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.6, 1.02), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
                  labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(24, 21, 22), 
                     labels = c("LK-Vendace   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Latitude (°N)", y = "Relativized ADD") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(0.75, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.all.corr <- grid.arrange(arrangeGrob(textGrob(""), 
                                          get_legend(plot.survival.corr.temp),
                                          nrow = 1,
                                          widths = c(0.09, 1)),
                              arrangeGrob(textGrob("A", x = 0.75, y = 0.95, gp = gpar(cex = 2, fontfamily = "Arial", fontface = "bold")),
                                          plot.survival.corr.temp + theme(legend.position = "none", axis.title.x = element_blank()), 
                                          plot.survival.corr.cgv + theme(legend.position = "none", axis.title.x = element_blank()), 
                                          textGrob("B", x = 0.75, y = 0.95, gp = gpar(cex = 2, fontfamily = "Arial", fontface = "bold")),
                                          plot.dpf.corr.temp + theme(legend.position = "none", axis.title.x = element_blank()),
                                          plot.dpf.corr.cgv + theme(legend.position = "none", axis.title.x = element_blank()),
                                          textGrob("C", x = 0.75, y = 0.95, gp = gpar(cex = 2, fontfamily = "Arial", fontface = "bold")),
                                          plot.ADD.corr.temp + theme(legend.position = "none", axis.title.x = element_blank()),
                                          plot.ADD.corr.cgv + theme(legend.position = "none", axis.title.x = element_blank()),
                                          nrow = 3,
                                          ncol = 3,
                                          widths = c(0.07, 1, 1),
                                          bottom = textGrob("Latitude (°N)", x = 0.545, gp = gpar(cex = 2, fontfamily = "Arial"))),
                              heights = c(0.025, 1)
                              )

ggsave("figures/embryo/2020-Embryo-LHT-Corr.png", plot = plot.all.corr, width = 12, height = 15, dpi = 300)


# VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
hatch.survival.glm.emm.confint <- hatch.survival.glm.emm.confint %>% 
  mutate(jit.y = ifelse(group == "LS-Cisco" & temperature %in% c("2.0°C", "4.5°C", "7.0°C"), (asymp.UCL * 100) - 2, (asymp.UCL * 100) + 2),
         jit.y = ifelse(group == "LO-Cisco" & temperature %in% c("2.0°C", "4.5°C"), (asymp.UCL * 100) + 3, jit.y),
         hjust = ifelse(group == "LS-Cisco" & temperature %in% c("2.0°C", "4.5°C", "7.0°C"), -0.4, 0.5))

plot.survival <- ggplot(hatch.survival.glm.emm.confint, aes(x = temperature, y = (prob * 100), group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.18)) +
  geom_point(size = 5, position = position_dodge(0.18)) +
  geom_errorbar(aes(ymin = (prob - SE) * 100, ymax = (prob + SE) * 100), 
                position = position_dodge(0.18),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_errorbar(aes(ymin = (asymp.LCL * 100), ymax = (asymp.UCL* 100)), 
  #              position = position_dodge(0.18),
  #              size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_text(aes(label = .group, y = jit.y, hjust = hjust), size = 4, 
  #          position = position_dodge(0.18), show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Mean Embryo Survival (%)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

#ggsave("figures/embryo/2020-Survival-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

## Days Post Fertilization
hatch.dpf.glm.emm.confint <- hatch.dpf.glm.emm.confint %>% 
  mutate(jit.y = ifelse(group == "LK-Whitefish" & temperature == "2.0°C" | group == "LS-Cisco" & temperature == "9.0°C",
                        lower.CL+2, upper.CL+3.5),
         jit.y = ifelse(group == "LS-Cisco" & temperature == "7.0°C", lower.CL+2, jit.y),
         hjust = ifelse(group == "LK-Whitefish" & temperature == "2.0°C", 1.9, 0.5),
         hjust = ifelse(group == "LS-Cisco" & temperature == "7.0°C", 2.3, hjust),
         hjust = ifelse(group == "LS-Cisco" & temperature == "9.0°C", 2.1, hjust))

plot.dpf <- ggplot(hatch.dpf.glm.emm.confint, aes(x = temperature, y = emmean, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.18)) +
  geom_point(size = 5, position = position_dodge(0.18)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.18),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
  #              position = position_dodge(0.18),
  #              size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_text(aes(label = .group, y = jit.y, hjust = hjust), size = 4, 
  #          position = position_dodge(0.18), show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(38, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Mean DPF", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

#ggsave("figures/embryo/2020-DPF-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

## Accumulated Degree-Days
hatch.ADD.glm.emm.confint <- hatch.ADD.glm.emm.confint %>% 
  mutate(jit.y = ifelse(group == "LO-Cisco" & temperature %in% c("2.0°C", "4.5°C", "7.0°C"), upper.CL-5, upper.CL+14),
         hjust = ifelse(group == "LO-Cisco" & temperature == "4.5°C", -0.8, 0.5),
         hjust = ifelse(group == "LO-Cisco" & temperature  == "2.0°C", -1.1, hjust),
         hjust = ifelse(group == "LO-Cisco" & temperature  == "7.0°C", -0.65, hjust))

plot.ADD <- ggplot(hatch.ADD.glm.emm.confint, aes(x = temperature, y = emmean, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.18)) +
  geom_point(size = 5, position = position_dodge(0.18)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                position = position_dodge(0.18),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
  #              position = position_dodge(0.18),
  #              size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  #geom_text(aes(label = .group, y = jit.y, hjust = hjust), size = 4, 
  #          position = position_dodge(0.18), show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(250, 850), breaks = seq(250, 850, 100), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Mean ADD (°C)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 15, 0, 0)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

#ggsave("figures/embryo/2020-ADD-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

## Combine all figures
plot.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                     get_legend(plot.survival),
                                     nrow = 1,
                                     widths = c(0.09, 1)),
                         arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()), 
                                     plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()),
                                     plot.ADD + theme(legend.position = "none", axis.title.x = element_blank()),
                                     nrow = 3,
                                     bottom = textGrob("Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
                         )

ggsave("figures/embryo/2020-Embryo-LHT-Confint.png", plot = plot.all, width = 12, height = 15, dpi = 300)

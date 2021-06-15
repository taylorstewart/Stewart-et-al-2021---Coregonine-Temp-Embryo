#### CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(emmeans)


#### LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD <- read_excel("data/Coregonine-Temperature-Experiment-NA-Hatch.xlsx", sheet = "temperature") %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA ------------------------------------------------------

hatch.NA <- read_excel("data/Coregonine-Temperature-Experiment-NA-Hatch.xlsx", sheet = "hatching") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, latitude, species, male, female, female_tl, female_fm, male_tl, male_fm, block, no, temperature, eye, hatch, dpf, ADD, include.incubation)

hatch.FI <- read_excel("data/Coregonine-Temperature-Experiment-FI-Hatch.xlsx", sheet = "hatching") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, latitude, species, male, female, female_tl, female_fm, male_tl, male_fm, block, no, temperature, eye, hatch, dpf, ADD, include.incubation)

## Combine all populations and years
hatch <- bind_rows(hatch.NA, hatch.FI) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         family = factor(paste0(female, male)),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")),
         trans.dpf = dpf^(1/3),
         trans.ADD = ADD^(1/3))

## Clean up environment
rm(hatch.NA, hatch.FI, ADD)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival.cisco <- hatch %>% filter(eye != 0, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.survival.finland <- hatch %>% filter(eye != 0, group %in% c( "LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## filter to only hatched embryos
hatch.dpf.cisco <- hatch %>% filter(!is.na(dpf), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels() %>% 
  filter(include.incubation == "y")
hatch.dpf.finland <- hatch %>% filter(!is.na(dpf), hatch == 1, group %in% c( "LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels() %>% 
  filter(include.incubation == "y")

## filter to only hatched embryos
hatch.ADD.cisco <- hatch %>% filter(!is.na(ADD), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco"))%>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels() %>% 
  filter(include.incubation == "y")
hatch.ADD.finland <- hatch %>% filter(!is.na(ADD), hatch == 1, group %in% c( "LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels() %>% 
  filter(include.incubation == "y")


#### STATISTICAL ANALYSIS - SURVIVAL - CISCO -----------------------------------------------------

## backward elimination to select best model
hatch.survival.cisco.glm <- buildmer(hatch ~ temperature + group + temperature:group + 
                                    (1|family) + (1|male) + (1|female) + (1|block), 
                                  direction = 'backward', data = hatch.survival.cisco, 
                                  family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.cisco.glm.formula <- formula(hatch.survival.cisco.glm@model))

## fit best model
hatch.survival.cisco.glm.final <- glmer(hatch.survival.cisco.glm.formula, data = hatch.survival.cisco, 
                                        family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.cisco.glm.formula, data = hatch.survival.cisco, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.cisco.glm.family <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                        (1 | female), data = hatch.survival.cisco, 
                                      family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.cisco.glm.female <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                        (1 | family), data = hatch.survival.cisco, 
                                      family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.cisco.glm.family, hatch.survival.cisco.glm.final)
# female
anova(hatch.survival.cisco.glm.female, hatch.survival.cisco.glm.final)


#### STATISTICAL ANALYSIS - SURVIVAL - FINLAND ---------------------------------------------------

## backward elimination to select best model
hatch.survival.finland.glm <- buildmer(hatch ~ temperature + group + temperature:group+
                                         (1|family) + (1|male) + (1|female) + (1|block), 
                                       direction = 'backward', data = hatch.survival.finland, 
                                       family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.finland.glm.formula <- formula(hatch.survival.finland.glm@model))

## fit best model
hatch.survival.finland.glm.final <- glmer(hatch.survival.finland.glm.formula, data = hatch.survival.finland, 
                                          family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.finland.glm.formula, data = hatch.survival.finland, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.finland.glm.family <- glmer(hatch ~ 1 + temperature + (1 | female), data = hatch.survival.finland, 
                                           family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.finland.glm.female <- glmer(hatch ~ 1 + temperature + (1 | family), data = hatch.survival.finland, 
                                           family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.finland.glm.family, hatch.survival.finland.glm.final)
# female
anova(hatch.survival.finland.glm.female, hatch.survival.finland.glm.final)

## Calculate estimated marginal means - be very patient!
hatch.survival.finland.glm.emm <- emmeans(hatch.survival.finland.glm.final, ~ temperature)

## Pairwise
pairs(hatch.survival.finland.glm.emm, simple = "temperature", adjust = "tukey", type = "response") 


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - CISCO --------------------------------------

## fit full model
hatch.dpf.cisco.glm.full <- lmer(trans.dpf ~ 1 + temperature + group + temperature:group + 
                                   (1|family) + (1|male) + (1|female) + (1|block), 
                                 data = hatch.dpf.cisco)

## backward elimination to select best model
hatch.dpf.cisco.glm <- step(hatch.dpf.cisco.glm.full)
( hatch.dpf.cisco.glm.formula <- get_model(hatch.dpf.cisco.glm)@call[["formula"]])

## fit best model
hatch.dpf.cisco.glm.final <- lmer(hatch.dpf.cisco.glm.formula, data = hatch.dpf.cisco)

## check residuals for normality
lattice::qqmath(hatch.dpf.cisco.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.dpf.cisco.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.cisco.glm.formula, data = hatch.dpf.cisco, method = "LRT")
rand(hatch.dpf.cisco.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - FINLAND ------------------------------------

## fit full model
hatch.dpf.finland.glm.full <- lmer(trans.dpf ~ 1 + temperature + group + temperature:group + 
                                     (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.dpf.finland)

## backward elimination to select best model
hatch.dpf.finland.glm <- step(hatch.dpf.finland.glm.full)
( hatch.dpf.finland.glm.formula <- get_model(hatch.dpf.finland.glm)@call[["formula"]])

## fit best model
hatch.dpf.finland.glm.final <- lmer(hatch.dpf.finland.glm.formula, data = hatch.dpf.finland)

## check residuals for normality
lattice::qqmath(hatch.dpf.finland.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.dpf.finland.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.finland.glm.formula, data = hatch.dpf.finland, method = "LRT")
rand(hatch.dpf.finland.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - CISCO --------------------------------------

## fit full model
hatch.ADD.cisco.glm.full <- lmer(trans.ADD ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.ADD.cisco)

## backward elimination to select best model
hatch.ADD.cisco.glm <- step(hatch.ADD.cisco.glm.full)
( hatch.ADD.cisco.glm.formula <- get_model(hatch.ADD.cisco.glm)@call[["formula"]])

## fit best model
hatch.ADD.cisco.glm.final <- lmer(hatch.ADD.cisco.glm.formula, data = hatch.ADD.cisco)

## check residuals for normality
lattice::qqmath(hatch.ADD.cisco.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.ADD.cisco.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.cisco.glm.formula, data = hatch.ADD.cisco, method = "LRT")
rand(hatch.ADD.cisco.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - FINLAND ------------------------------------

## fit full model
hatch.ADD.finland.glm.full <- lmer(trans.ADD ~ 1 + temperature + group + temperature:group + 
                                     (1|family) + (1|male) + (1|female) + (1|block), 
                                   data = hatch.ADD.finland)

## backward elimination to select best model
hatch.ADD.finland.glm <- step(hatch.ADD.finland.glm.full)
( hatch.ADD.finland.glm.formula <- get_model(hatch.ADD.finland.glm)@call[["formula"]])

## fit best model
hatch.ADD.finland.glm.final <- lmer(hatch.ADD.finland.glm.formula, data = hatch.ADD.finland)

## check residuals for normality
lattice::qqmath(hatch.ADD.finland.glm.final, id = 0.1, idLabels = ~.obs)
hist(rstudent(hatch.ADD.finland.glm.final))

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.finland.glm.formula, data = hatch.ADD.finland, method = "LRT")
rand(hatch.ADD.finland.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

temp <- data.frame(group = c("LK-Whitefish", "LK-Whitefish", "LK-Whitefish", "LK-Whitefish",
                             "LK-Vendace", "LK-Vendace", "LK-Vendace", "LK-Vendace",
                             "LS-Cisco", "LS-Cisco", "LS-Cisco", "LS-Cisco",
                             "LO-Cisco", "LO-Cisco", "LO-Cisco", "LO-Cisco"),
                   temperature = c(rep(c(2.2, 4.0, 6.9, 8.0),2), rep(c(2.0, 4.4, 6.9, 8.9),2)),
                   temp.treatment = factor(rep(c("Coldest", "Cold", "Warm", "Warmest"), 4), 
                                           ordered = TRUE, levels = c("Coldest", "Cold", "Warm", "Warmest")))

## Embryo Survival Overall
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.hatch = mean(hatch),
            se.hatch = sd(hatch)/sqrt(n())) %>% 
  group_by(temperature) %>% 
  mutate(width = 0.15 * n())

## Embryo Survival - Standardized Within Family
hatch.survival.summary.family <- hatch %>% filter(eye != 0) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.hatch = mean(hatch)) %>% ungroup()

hatch.survival.stand <- hatch.survival.summary.family %>% filter(temperature %in% c(2, 2.2)) %>% 
  select(group, family, local.survival = mean.hatch)

hatch.survival.summary.stand <- hatch.survival.summary.family %>% left_join(hatch.survival.stand) %>% 
  filter(group != "LK-Whitefish" | family != "F8M11") %>%  ## No data at 2C
  mutate(survival.diff = 100*(1+(mean.hatch-local.survival)/local.survival)) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.survival.diff = mean(survival.diff),
            se.survival.diff = sd(survival.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.survival.diff = ifelse(se.survival.diff == 0, NA, se.survival.diff),
         percent.loss = 100-mean.survival.diff,
         group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))


## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.dpf = mean(dpf),
            se.dpf = sd(dpf)/sqrt(n())) %>% ungroup() %>% 
  group_by(temperature) %>% 
  mutate(width = 0.15 * n())

## Days Post Fertilization - Standardized Within Family
hatch.dpf.summary.family <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.dpf = mean(dpf)) %>% ungroup()

hatch.dpf.stand <- hatch.dpf.summary.family %>% filter(temperature %in% c(2, 2.2)) %>% 
  select(group, family, local.dpf = mean.dpf)

hatch.dpf.summary.stand <- hatch.dpf.summary.family %>% left_join(hatch.dpf.stand) %>% 
  filter(group != "LK-Whitefish" | family != "F8M11") %>%  ## No data at 2C
  mutate(dpf.diff = 100*(1+(mean.dpf-local.dpf)/local.dpf)) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.dpf.diff = mean(dpf.diff),
            se.dpf.diff = sd(dpf.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.dpf.diff = ifelse(se.dpf.diff == 0, NA, se.dpf.diff),
         percent.loss = 100-mean.dpf.diff,
         group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.ADD = mean(ADD),
            se.ADD = sd(ADD)/sqrt(n())) %>% ungroup() %>% 
  group_by(temperature) %>% 
  mutate(width = 0.15 * n())

## Accumulated Degree-Days - Standardized Within Family
hatch.ADD.summary.family <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, temperature, group, family) %>% 
  summarize(mean.ADD = mean(ADD)) %>% ungroup()

hatch.ADD.stand <- hatch.ADD.summary.family %>% filter(temperature %in% c(2, 2.2)) %>% 
  select(group, family, local.ADD = mean.ADD)

hatch.ADD.summary.stand <- hatch.ADD.summary.family %>% left_join(hatch.ADD.stand) %>% 
  filter(group != "LK-Whitefish" | family != "F8M11") %>%  ## No data at 2C
  mutate(ADD.diff = 100*(1+(mean.ADD-local.ADD)/local.ADD)) %>%
  group_by(population, temperature, group) %>% 
  summarize(mean.ADD.diff = mean(ADD.diff),
            se.ADD.diff = sd(ADD.diff)/sqrt(n())) %>% 
  left_join(temp) %>% 
  mutate(se.ADD.diff = ifelse(se.ADD.diff == 0, NA, se.ADD.diff),
         percent.loss = 100-mean.ADD.diff,
         group = factor(group, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))


#### VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
plot.survival <- ggplot(hatch.survival.summary, aes(x = temperature, y = (mean.hatch * 100), 
                                                    group = group, color = group, shape = group, 
                                                    linetype = group, width = width)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = (mean.hatch - se.hatch) * 100, ymax = (mean.hatch + se.hatch) * 100), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean Embryo Survival (%)", x = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Plot Standardized Survival
plot.survival.stand <- ggplot(hatch.survival.summary.stand, aes(x = group, y = mean.survival.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.survival.diff - se.survival.diff), ymax = (mean.survival.diff + se.survival.diff)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 105), breaks = seq(0.0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized Survival (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 


## Days Post Fertilization
plot.dpf <- ggplot(hatch.dpf.summary, aes(x = temperature, y = mean.dpf, 
                                          group = group, color = group, shape = group, 
                                          linetype = group, width = width)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.dpf - se.dpf, ymax = mean.dpf + se.dpf), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(30, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean DPF") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Plot Standardized DPF
plot.dpf.stand <- ggplot(hatch.dpf.summary.stand, aes(x = group, y = mean.dpf.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.dpf.diff - se.dpf.diff), ymax = (mean.dpf.diff + se.dpf.diff)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 102), breaks = seq(0.0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized DPF (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 

## Accumulated Degree-Days
plot.ADD <- ggplot(hatch.ADD.summary, aes(x = temperature, y = mean.ADD, 
                                          group = group, color = group, shape = group, 
                                          linetype = group, width = width)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.ADD - se.ADD, ymax = mean.ADD + se.ADD), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(250, 850), breaks = seq(250, 850, 100), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean ADD") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Plot Standardized ADD
plot.ADD.stand <- ggplot(hatch.ADD.summary.stand, aes(x = group, y = mean.ADD.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = (mean.ADD.diff - se.ADD.diff), ymax = (mean.ADD.diff + se.ADD.diff)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 205), breaks = seq(80, 200, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  coord_cartesian(ylim = c(80, 205)) +
  labs(y = "Standardized ADD (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 


## Combine all figures
plot.all <- grid.arrange(
  arrangeGrob(
    arrangeGrob(textGrob(""),
                get_legend(plot.survival),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.survival.stand),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(1, 0.7)
    ),
  arrangeGrob(
    arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ADD + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 3,
                bottom = textGrob("Mean Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 0.8, fontfamily = "Arial"))),
    arrangeGrob(plot.survival.stand + theme(legend.position = "none", axis.title.x = element_blank()), 
                plot.dpf.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ADD.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 3,
                bottom = textGrob("Study Group", x = 0.55, gp = gpar(cex = 0.8, fontfamily = "Arial"))),
    ncol = 2,
    widths = c(1, 0.7)
    ),
  heights = c(0.035, 1.1)
)

ggsave("figures/Fig3.tiff", plot = plot.all, width = 6.9, height = 6.9, dpi = 600)


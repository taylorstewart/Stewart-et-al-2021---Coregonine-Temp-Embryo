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


#### LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD <- read.csv("data/Artedi-Temperature-ADD-2020.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


#### LOAD HATCHING DATA ------------------------------------------------------

hatch.NA <- read_excel("data/Coregonine-Temperature-Experiment-NA-Hatch.xlsx", sheet = "2020HatchingData") %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD) %>% 
  dplyr::select(population, latitude, species, male, female, block, temperature, eye, hatch, dpf, ADD)

hatch.FI <- read_excel("data/Coregonine-Temperature-Experiment-FI-Hatch.xlsx", sheet = "2019HatchingData") %>% 
  mutate(premature = 0) %>% 
  dplyr::select(population, latitude, species, male, female, block, temperature, eye, hatch, dpf, ADD)

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
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Clean up environment
rm(hatch.NA, hatch.FI, ADD)


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

## filter to only eyed embryos
hatch.survival.cisco <- hatch %>% filter(eye != 0, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.survival.vendace <- hatch %>% filter(eye != 0, group == "LK-Vendace") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()
hatch.survival.whitefish <- hatch %>% filter(eye != 0, group == "LK-Whitefish") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## filter to only hatched embryos
hatch.dpf.cisco <- hatch %>% filter(!is.na(dpf), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.dpf.vendace <- hatch %>% filter(!is.na(dpf), hatch == 1, group == "LK-Vendace") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()
hatch.dpf.whitefish <- hatch %>% filter(!is.na(dpf), hatch == 1, group == "LK-Whitefish") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## filter to only hatched embryos
hatch.ADD.cisco <- hatch %>% filter(!is.na(ADD), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco"))%>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.ADD.vendace <- hatch %>% filter(!is.na(ADD), hatch == 1, group == "LK-Vendace") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()
hatch.ADD.whitefish <- hatch %>% filter(!is.na(ADD), hatch == 1, group == "LK-Whitefish") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()


# STATISTICAL ANALYSIS - SURVIVAL - CISCO -----------------------------------------------------

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


#### STATISTICAL ANALYSIS - SURVIVAL - VENDACE ---------------------------------------------------

## backward elimination to select best model
hatch.survival.vendace.glm <- buildmer(hatch ~ temperature +
                                         (1|family) + (1|male) + (1|female) + (1|block), 
                                       direction = 'backward', data = hatch.survival.vendace, 
                                       family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.vendace.glm.formula <- formula(hatch.survival.vendace.glm@model))

## fit best model
hatch.survival.vendace.glm.final <- glmer(hatch.survival.vendace.glm.formula, data = hatch.survival.vendace, 
                                          family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.vendace.glm.formula, data = hatch.survival.vendace, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.vendace.glm.family <- glmer(hatch ~ 1 + temperature + (1 | female), data = hatch.survival.vendace, 
                                           family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.vendace.glm.female <- glmer(hatch ~ 1 + temperature + (1 | family), data = hatch.survival.vendace, 
                                           family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.vendace.glm.family, hatch.survival.vendace.glm.final)
# female
anova(hatch.survival.vendace.glm.female, hatch.survival.vendace.glm.final)


#### STATISTICAL ANALYSIS - SURVIVAL - WHITEFISH -------------------------------------------------

## backward elimination to select best model
hatch.survival.whitefish.glm <- buildmer(hatch ~ temperature +
                                         (1|family) + (1|male) + (1|female) + (1|block), 
                                       direction = 'backward', data = hatch.survival.whitefish, 
                                       family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.whitefish.glm.formula <- formula(hatch.survival.whitefish.glm@model))

## fit best model
hatch.survival.whitefish.glm.final <- glmer(hatch.survival.whitefish.glm.formula, data = hatch.survival.whitefish,
                                            family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.whitefish.glm.formula, data = hatch.survival.whitefish, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.whitefish.glm.family <- glmer(hatch ~ 1 + temperature + (1 | female) + (1 | male), data = hatch.survival.whitefish, 
                                             family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.whitefish.glm.female <- glmer(hatch ~ 1 + temperature + (1 | family) + (1 | male), data = hatch.survival.whitefish, 
                                             family = binomial, control = glmerControl(optimizer = "bobyqa"))
# male
hatch.survival.whitefish.glm.male <- glmer(hatch ~ 1 + temperature + (1 | family) + (1 | female), data = hatch.survival.whitefish, 
                                             family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.whitefish.glm.family, hatch.survival.whitefish.glm.final)
# female
anova(hatch.survival.whitefish.glm.female, hatch.survival.whitefish.glm.final)
# male
anova(hatch.survival.whitefish.glm.male, hatch.survival.whitefish.glm.final)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - CISCO --------------------------------------

## fit full model
hatch.dpf.cisco.glm.full <- lmer(dpf ~ 1 + temperature + group + temperature:group + 
                                   (1|family) + (1|male) + (1|female) + (1|block), 
                                 data = hatch.dpf.cisco)

## backward elimination to select best model
hatch.dpf.cisco.glm <- step(hatch.dpf.cisco.glm.full)
( hatch.dpf.cisco.glm.formula <- get_model(hatch.dpf.cisco.glm)@call[["formula"]])

## fit best model
hatch.dpf.cisco.glm.final <- lmer(hatch.dpf.cisco.glm.formula, data = hatch.dpf.cisco)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.cisco.glm.formula, data = hatch.dpf.cisco, method = "LRT")
rand(hatch.dpf.cisco.glm.final)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - VENDACE ------------------------------------

## fit full model
hatch.dpf.vendace.glm.full <- lmer(dpf ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.dpf.vendace)

## backward elimination to select best model
hatch.dpf.vendace.glm <- step(hatch.dpf.vendace.glm.full)
( hatch.dpf.vendace.glm.formula <- get_model(hatch.dpf.vendace.glm)@call[["formula"]])

## fit best model
hatch.dpf.vendace.glm.final <- lmer(hatch.dpf.vendace.glm.formula, data = hatch.dpf.vendace)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.vendace.glm.formula, data = hatch.dpf.vendace, method = "LRT")
rand(hatch.dpf.vendace.glm.final)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - WHITEFISH ----------------------------------

## fit full model
hatch.dpf.whitefish.glm.full <- lmer(dpf ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                   data = hatch.dpf.whitefish)

## backward elimination to select best model
hatch.dpf.whitefish.glm <- step(hatch.dpf.whitefish.glm.full)
( hatch.dpf.whitefish.glm.formula <- get_model(hatch.dpf.whitefish.glm)@call[["formula"]])

## fit best model
hatch.dpf.whitefish.glm.final <- lmer(hatch.dpf.whitefish.glm.formula, data = hatch.dpf.whitefish)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.whitefish.glm.formula, data = hatch.dpf.whitefish, method = "LRT")
rand(hatch.dpf.whitefish.glm.final)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - CISCO --------------------------------------

## fit full model
hatch.ADD.cisco.glm.full <- lmer(ADD ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.ADD.cisco)

## backward elimination to select best model
hatch.ADD.cisco.glm <- step(hatch.ADD.cisco.glm.full)
( hatch.ADD.cisco.glm.formula <- get_model(hatch.ADD.cisco.glm)@call[["formula"]])

## fit best model
hatch.ADD.cisco.glm.final <- lmer(hatch.ADD.cisco.glm.formula, data = hatch.ADD.cisco)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.cisco.glm.formula, data = hatch.ADD.cisco, method = "LRT")
rand(hatch.ADD.cisco.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - VENDACE ------------------------------------

## fit full model
hatch.ADD.vendace.glm.full <- lmer(ADD ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                 data = hatch.ADD.vendace)

## backward elimination to select best model
hatch.ADD.vendace.glm <- step(hatch.ADD.vendace.glm.full)
( hatch.ADD.vendace.glm.formula <- get_model(hatch.ADD.vendace.glm)@call[["formula"]])

## fit best model
hatch.ADD.vendace.glm.final <- lmer(hatch.ADD.vendace.glm.formula, data = hatch.ADD.vendace)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.vendace.glm.formula, data = hatch.ADD.vendace, method = "LRT")
rand(hatch.ADD.vendace.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - WHITEFISH ----------------------------------

## fit full model
hatch.ADD.whitefish.glm.full <- lmer(ADD ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                 data = hatch.ADD.whitefish)

## backward elimination to select best model
hatch.ADD.whitefish.glm <- step(hatch.ADD.whitefish.glm.full)
( hatch.ADD.whitefish.glm.formula <- get_model(hatch.ADD.whitefish.glm)@call[["formula"]])

## fit best model
hatch.ADD.whitefish.glm.final <- lmer(hatch.ADD.whitefish.glm.formula, data = hatch.ADD.whitefish)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.whitefish.glm.formula, data = hatch.ADD.whitefish, method = "LRT")
rand(hatch.ADD.whitefish.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

## Embryo Survival
hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.hatch = mean(hatch),
            se.hatch = sd(hatch)/sqrt(n()))

## Days Post Fertilization
hatch.dpf.summary <- hatch %>% filter(!is.na(dpf), hatch == 1) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.dpf = mean(dpf),
            se.dpf = sd(dpf)/sqrt(n()))

## Accumulated Degree-Days
hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.ADD = mean(ADD),
            se.ADD = sd(ADD)/sqrt(n()))


# VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
plot.survival <- ggplot(hatch.survival.summary, aes(x = temperature, y = (mean.hatch * 100), group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = (mean.hatch - se.hatch) * 100, ymax = (mean.hatch + se.hatch) * 100), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = "Mean ES (%)") +
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

## Days Post Fertilization
plot.dpf <- ggplot(hatch.dpf.summary, aes(x = temperature, y = mean.dpf, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.dpf - se.dpf, ymax = mean.dpf + se.dpf), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(30, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = "Mean DPF") +
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

## Accumulated Degree-Days
plot.ADD <- ggplot(hatch.ADD.summary, aes(x = temperature, y = mean.ADD, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 5, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.ADD - se.ADD, ymax = mean.ADD + se.ADD), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(250, 850), breaks = seq(250, 850, 100), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = "Mean ADD (°C)") +
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

## Combine all figures
plot.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                     get_legend(plot.survival),
                                     nrow = 1,
                                     widths = c(0.09, 1)),
                         arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()), 
                                     plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()),
                                     plot.ADD + theme(legend.position = "none", axis.title.x = element_blank()),
                                     nrow = 3,
                                     bottom = textGrob("Mean Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 1.75, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
)

ggsave("figures/2020-Embryo-LHT-SE.png", plot = plot.all, width = 11, height = 15, dpi = 300)

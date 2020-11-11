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

ADD <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
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
hatch.survival.na <- hatch %>% filter(eye != 0, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.survival.fi <- hatch %>% filter(eye != 0, group %in% c("LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## filter to only hatched embryos
hatch.dpf.na <- hatch %>% filter(!is.na(dpf), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.dpf.fi <- hatch %>% filter(!is.na(dpf), hatch == 1, group %in% c("LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## filter to only hatched embryos
hatch.ADD.na <- hatch %>% filter(!is.na(ADD), hatch == 1, group %in% c("LO-Cisco", "LS-Cisco"))%>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
hatch.ADD.fi <- hatch %>% filter(!is.na(ADD), hatch == 1, group %in% c("LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()


#### STATISTICAL ANALYSIS - SURVIVAL - NA ------------------------------------

## backward elimination to select best model
hatch.survival.na.glm <- buildmer(hatch ~ temperature + group + temperature:group + 
                                    (1|family) + (1|male) + (1|female) + (1|block), 
                                  direction = 'backward', data = hatch.survival.na, 
                                  family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.na.glm.formula <- formula(hatch.survival.na.glm@model))

## fit best model
hatch.survival.na.glm.final <- glmer(hatch.survival.na.glm.formula, data = hatch.survival.na, 
                                     family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.na.glm.formula, data = hatch.survival.na, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.na.glm.family <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                        (1 | female), data = hatch.survival.na, 
                                      family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.na.glm.female <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                        (1 | family), data = hatch.survival.na, 
                                      family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.na.glm.family, hatch.survival.na.glm.final)
# female
anova(hatch.survival.na.glm.female, hatch.survival.na.glm.final)


#### STATISTICAL ANALYSIS - SURVIVAL - FI ------------------------------------

## backward elimination to select best model
hatch.survival.fi.glm <- buildmer(hatch ~ temperature + group + temperature:group + 
                                    (1|family) + (1|male) + (1|female) + (1|block), 
                                  direction = 'backward', data = hatch.survival.fi, 
                                  family = binomial, control = glmerControl(optimizer = "bobyqa"))
( hatch.survival.fi.glm.formula <- formula(hatch.survival.fi.glm@model))

## fit best model
hatch.survival.fi.glm.final <- glmer(hatch.survival.fi.glm.formula, data = hatch.survival.fi, 
                                     family = binomial, control = glmerControl(optimizer = "bobyqa"))

## likelihood ratio test for fixed effects
mixed(hatch.survival.fi.glm.formula, data = hatch.survival.fi, method = "LRT")

## fit model without random effects for LRT
# family
hatch.survival.fi.glm.family <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                       (1 | female), data = hatch.survival.fi, 
                                     family = binomial, control = glmerControl(optimizer = "bobyqa"))
# female
hatch.survival.fi.glm.female <- glmer(hatch ~ 1 + temperature + group + temperature:group + 
                                        (1 | family), data = hatch.survival.fi, 
                                      family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Compare full to reduced models (LRT)
# family
anova(hatch.survival.fi.glm.family, hatch.survival.fi.glm.final)
# female
anova(hatch.survival.fi.glm.female, hatch.survival.fi.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - NA ---------------------

## fit full model
hatch.dpf.na.glm.full <- lmer(dpf ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.dpf.na)

## backward elimination to select best model
hatch.dpf.na.glm <- step(hatch.dpf.na.glm.full)
( hatch.dpf.na.glm.formula <- get_model(hatch.dpf.na.glm)@call[["formula"]])

## fit best model
hatch.dpf.na.glm.final <- lmer(hatch.dpf.na.glm.formula, data = hatch.dpf.na)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.na.glm.formula, data = hatch.dpf.na, method = "LRT")
rand(hatch.dpf.na.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - FI ---------------------

## fit full model
hatch.dpf.fi.glm.full <- lmer(dpf ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.dpf.fi)

## backward elimination to select best model
hatch.dpf.fi.glm <- step(hatch.dpf.fi.glm.full)
( hatch.dpf.fi.glm.formula <- get_model(hatch.dpf.fi.glm)@call[["formula"]])

## fit best model
hatch.dpf.fi.glm.final <- lmer(hatch.dpf.fi.glm.formula, data = hatch.dpf.fi)

## likelihood ratio test for fixed and random effects
mixed(hatch.dpf.fi.glm.formula, data = hatch.dpf.fi, method = "LRT")
rand(hatch.dpf.fi.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - NA ---------------------

## fit full model
hatch.ADD.na.glm.full <- lmer(ADD ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.ADD.na)

## backward elimination to select best model
hatch.ADD.na.glm <- step(hatch.ADD.na.glm.full)
( hatch.ADD.na.glm.formula <- get_model(hatch.ADD.na.glm)@call[["formula"]])

## fit best model
hatch.ADD.na.glm.final <- lmer(hatch.ADD.na.glm.formula, data = hatch.ADD.na)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.na.glm.formula, data = hatch.ADD.na, method = "LRT")
rand(hatch.ADD.na.glm.final)


#### STATISTICAL ANALYSIS - INCUBATION PERIOD (ADD) - FI ---------------------

## fit full model
hatch.ADD.fi.glm.full <- lmer(ADD ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = hatch.ADD.fi)

## backward elimination to select best model
hatch.ADD.fi.glm <- step(hatch.ADD.fi.glm.full)
( hatch.ADD.fi.glm.formula <- get_model(hatch.ADD.fi.glm)@call[["formula"]])

## fit best model
hatch.ADD.fi.glm.final <- lmer(hatch.ADD.fi.glm.formula, data = hatch.ADD.fi)

## likelihood ratio test for fixed and random effects
mixed(hatch.ADD.fi.glm.formula, data = hatch.ADD.fi, method = "LRT")
rand(hatch.ADD.fi.glm.final)


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

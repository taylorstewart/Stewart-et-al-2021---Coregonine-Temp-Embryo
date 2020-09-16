## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(tidyverse)
library(readxl)
library(lubridate)
library(magrittr)
library(ggplot2)
library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")


## ===========================================================
## Load and Initial Manipulations of the Larval Hatch Length Data
## ===========================================================
larvae.yolk.lo.2 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LO-2")
larvae.yolk.ls.2 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LS-2")
larvae.yolk.lk.2 <- read_excel("data/Finland2019-Measurements.xlsx", sheet = "LK-2")
larvae.yolk.lo.4.5 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LO-4.5")
larvae.yolk.ls.4.5 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LS-4.5")
larvae.yolk.lk.4.5 <- read_excel("data/Finland2019-Measurements.xlsx", sheet = "LK-4.5")
larvae.yolk.lo.7.0 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LO-7")
larvae.yolk.ls.7.0 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LS-7")
larvae.yolk.lk.7.0 <- read_excel("data/Finland2019-Measurements.xlsx", sheet = "LK-7")
larvae.yolk.lo.9.0 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LO-9")
larvae.yolk.ls.9.0 <- read_excel("data/Artedi-Temperature-Experiments-LarvaeMeasurements-2020.xlsx", sheet = "LS-9")
larvae.yolk.lk.9.0 <- read_excel("data/Finland2019-Measurements.xlsx", sheet = "LK-9")

larvae.yolk <- bind_rows(larvae.yolk.lo.2, larvae.yolk.ls.2, larvae.yolk.lk.2,
                         larvae.yolk.lo.4.5, larvae.yolk.ls.4.5, larvae.yolk.lk.4.5,
                         larvae.yolk.lo.7.0, larvae.yolk.ls.7.0, larvae.yolk.lk.7.0,
                         larvae.yolk.lo.9.0, larvae.yolk.ls.9.0, larvae.yolk.lk.9.0) %>% 
  filter(!is.na(y_vol_mm3)) %>% 
  mutate(population = factor(population, ordered = TRUE, levels = c("konnevesi", "superior", "ontario")),
         temperature = factor(temperature),
         species = factor(species),
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi")))

## -----------------------------------------------------------
## Calculate means and std. error for each treatment
## -----------------------------------------------------------
larvae.yolk.summary <- larvae.yolk %>% group_by(population, temperature, species, group) %>% 
  summarize(mean.yolk = mean(y_vol_mm3),
            sd.yolk = sd(y_vol_mm3),
            n = n(),
            se.yolk = sd.yolk/sqrt(n)) %>% 
  arrange(population, temperature)


ggplot(filter(larvae.yolk.summary, species %in% c("albula", "artedi")), aes(x = temperature, y = mean.yolk, group = population, color = population, shape = population)) + 
  geom_point(size = 3.5, position = position_dodge(0.09)) + 
  geom_line(size = 1, position = position_dodge(0.09)) +
  geom_errorbar(aes(ymin = mean.yolk - se.yolk, ymax = mean.yolk + se.yolk), linetype = "solid", width = 0.3, size = 0.85, position = position_dodge(0.09)) +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, 0.2), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.1)) +
  scale_color_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c("#33a02c", "#1f78b4", "#a6cee3")) +
  scale_shape_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c(17, 16, 15)) +
  labs(y = expression("Mean Yolk-sac Volume (mm"^3*" ± SE)"), x = 'Incubation Temperature (°C)') +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/larvae/2020-Yolk.png", width = 6.5, height = 7, dpi = 300)


ggplot(larvae.yolk.summary, aes(x = temperature, y = mean.yolk, group = group, color = group, shape = group)) + #, linetype = group)) + 
  geom_point(size = 3.5, position = position_dodge(0.09)) + 
  geom_line(size = 1, position = position_dodge(0.09)) +
  geom_errorbar(aes(ymin = mean.yolk - se.yolk, ymax = mean.yolk + se.yolk), linetype = "solid", width = 0.3, size = 0.85, position = position_dodge(0.09)) +
  scale_y_continuous(limits = c(0.0, 1.4), breaks = seq(0.0, 1.4, 0.2), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.1)) +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual(values = c(17, 18, 16, 15), 
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  #scale_linetype_manual(values = c("dashed", "dotdash", "solid", "solid"),
  #                      labels = c("L. Konnevesi (C. albula)    ", "L. Konnevesi (C. lavaretus)    ", "L. Superior (C. artedi)    ", "L. Ontario (C. artedi)")) +
  labs(y = expression("Mean Yolk-sac Volume (mm"^3*" ± SE)"), x = 'Incubation Temperature (°C)') +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/larvae/2020-Yolk-wWhitefish.png", width = 6.5, height = 7, dpi = 300)


##############################################################
## ANALYSIS
##############################################################
## -----------------------------------------------------------
## Fit model
## -----------------------------------------------------------
glm <- lmer(length_mm ~ temperature + group + temperature * group + (1|male) + (1|female), 
            data = larvae.yolk, REML = FALSE)

dg1 <- dredge(glm)                    # to select all model based on AICc
dg1

best <- get.models(dg1, 1)[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ temperature * group)
(best.emm.pair <- pairs(best.emm, simple = list("temperature", "group"), adjust = "bonferroni"))


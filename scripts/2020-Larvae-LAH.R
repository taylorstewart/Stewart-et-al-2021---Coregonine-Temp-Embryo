# Clear the environment first ---------------------------------------------

rm(list = ls(all.names = TRUE))

# Load packages -----------------------------------------------------------

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

# Load and Initial Manipulations of the Larval Hatch Length Data ----------

larvae.tl.lo.2 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-2")
larvae.tl.ls.2 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-2")
larvae.tl.lk.2 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-2")
larvae.tl.lo.4.5 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-4.5")
larvae.tl.ls.4.5 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-4.5")
larvae.tl.lk.4.5 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-4.5")
larvae.tl.lo.7.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-7")
larvae.tl.ls.7.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-7")
larvae.tl.lk.7.0 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-7")
larvae.tl.lo.9.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-9")
larvae.tl.ls.9.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-9")
larvae.tl.lk.9.0 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-9")

## Combine each population, temperature, and species
larvae.tl <- bind_rows(larvae.tl.lo.2, larvae.tl.ls.2, larvae.tl.lk.2,
                       larvae.tl.lo.4.5, larvae.tl.ls.4.5, larvae.tl.lk.4.5,
                       larvae.tl.lo.7.0, larvae.tl.ls.7.0, larvae.tl.lk.7.0,
                       larvae.tl.lo.9.0, larvae.tl.ls.9.0, larvae.tl.lk.9.0) %>% 
  filter(!is.na(length_mm)) %>% 
  mutate(population = factor(population, ordered = TRUE, levels = c("konnevesi", "superior", "ontario")),
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

# Calculate summary statistics --------------------------------------------

larvae.tl.summary <- larvae.tl %>% group_by(population, temperature, species, group) %>% 
  summarize(mean.tl = mean(length_mm),
            sd.tl = sd(length_mm),
            n = n(),
            se.tl = sd.tl/sqrt(n)) %>% 
  arrange(population, temperature)

# Visualizations ----------------------------------------------------------

ggplot(filter(larvae.tl.summary, species %in% c("albula", "artedi")), aes(x = temperature, y = mean.tl, group = population, color = population, shape = population)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean.tl - se.tl, ymax = mean.tl + se.tl), width = 0.1, size = 0.85) +
  scale_y_continuous(limits = c(6.5, 11.5), breaks = seq(6.5, 11.5, 1.0), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.1)) +
  scale_color_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c("#33a02c", "#1f78b4", "#a6cee3")) +
  scale_shape_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c(17, 16, 15)) +
  labs(y = "Mean Length-at-Hatch (mm ± SE)", x = 'Incubation Temperature (°C)') +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/larvae/2020-LAH.png", width = 6.5, height = 7, dpi = 300)



ggplot(larvae.tl.summary, aes(x = temperature, y = mean.tl, group = group, color = group, shape = group)) + 
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean.tl - se.tl, ymax = mean.tl + se.tl), width = 0.1, size = 0.85) +
  scale_y_continuous(limits = c(6.5, 11.5), breaks = seq(6.5, 11.5, 1.0), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.0, 0.1)) +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual(values = c(17, 18, 16, 15), 
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(y = "Mean Length-at-Hatch (mm ± SE)", x = 'Incubation Temperature (°C)') +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/larvae/2020-LAH-wWhitefish.png", width = 6.5, height = 7, dpi = 300)

# Statistical Analyses ----------------------------------------------------

## Fit model
glm <- lmer(length_mm ~ temperature  + group + temperature * group +       # Fixed
            (1|male) + (1|female),                                         # Random
            data = larvae.tl, 
            REML = FALSE)

dg1 <- dredge(glm)                    # to select all model based on AICc
dg1

best <- get.models(dg1, 1)[[1]]      # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ temperature * group)
pairs(best.emm, simple = list("group", "temperature"), adjust = "fdr") 

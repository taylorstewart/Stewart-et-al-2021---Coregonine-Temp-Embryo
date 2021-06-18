#### Clear the environment first ---------------------------------------------
rm(list = ls(all.names = TRUE))


#### Load packages -----------------------------------------------------------
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)
library(emmeans)

#### Load embryo diameter data -----------------------------------------------
data <- read_excel("data/Coregonine-Temperature-Experiment-EggMeasurements.xlsx", sheet = "egg_diameter_mm", skip = 35) %>% 
  mutate(female = factor(female),
         population = factor(population, levels = c("konnevesi", "superior", "ontario"))) %>% 
  filter(fert_success == "y")


data.summary <- data %>% group_by(population) %>% 
  summarize(mean.diam = mean(dia_mm),
            sd.diam = sd(dia_mm),
            n = n())

## Fit linear model and run ANOVA
lm.embryo <- lm(dia_mm ~ 0 + population, data = data)
summary(lm.embryo)
anova(lm.embryo)

## Calculate estimated margin means
embryo.emm <- emmeans(lm.embryo, ~ population)

## Pairwise
pairs(embryo.emm, simple = "population", adjust = "tukey", type = "response") 


#### Visualizations ----------------------------------------------------------
ggplot(data, aes(x = population, y = dia_mm)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1.25, 2.5), breaks = seq(1.25, 2.5, 0.25), expand = c(0, 0)) +
  scale_x_discrete(labels = c("LK-Vendace", "LS-Cisco", "LO-Cisco")) +
  #scale_fill_grey("combine", start = 0.0, end = 0.8,
  #                labels = c("Konnevesi", "Superior", "Ontario")) +
  labs(x = "", y = "Egg Diameter (mm)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))


ggsave("figures/adult/2020-EmbryoDiameter.png", width = 5.5, height = 7, dpi = 300)

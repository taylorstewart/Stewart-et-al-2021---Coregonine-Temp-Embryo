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
  dplyr::select(year, population, species, family, male, female, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.albula <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi vendace") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.lavaretus <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi whitefish") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, temperature, eye, premature, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, levels = c(2, 4.5, 7, 9), ordered = TRUE),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi")))

# Summarize Data ----------------------------------------------------------

hatch.survival.summary <- hatch %>% filter(eye != 0) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.survival = (mean(hatch)*100),
            sd.survival = (sd(hatch)*100),
            se.survival = sd.survival / sqrt(n()))

hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1, premature != 1) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.ADD = mean(ADD),
            sd.ADD = sd(ADD),
            se.ADD = sd.ADD / sqrt(n()))

hatch.dpf.summary <- hatch %>% filter(hatch == 1, premature != 1, !is.na(dpf)) %>% 
  group_by(population, species, temperature, group) %>% 
  summarize(mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf / sqrt(n()))

# Visualizations ----------------------------------------------------------

ggplot(hatch.survival.summary, aes(x = temperature, y = mean.survival, group = group, color = group, shape = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.survival-se.survival, ymax = mean.survival+se.survival), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(10, 100), breaks = seq(10, 100, 10), expand = c(0, 0)) +
  scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(17, 18, 16, 15), 
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(x = "Incubation Temperature (°C)", y = "Embryo Survival (% ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-Survival-wWhitefish.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.dpf.summary, aes(x = temperature, y = mean.dpf, group = group, color = group, shape = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.dpf-se.dpf, ymax = mean.dpf+se.dpf), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(45, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(17, 18, 16, 15), 
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (No. Days ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        #legend.position = c(0.8, 0.86),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-DPF-wWhitefish.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.ADD.summary, aes(x = temperature, y = mean.ADD, group = group, color = group, shape = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.ADD-sd.ADD, ymax = mean.ADD+sd.ADD), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(250, 950), breaks = seq(250, 950, 100), expand = c(0, 0)) +
  scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_shape_manual("combine", values = c(17, 18, 16, 15), 
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (ADD °C ± SE)", color = "Populations") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        #legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        #legend.title = element_text(face = "bold", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        #legend.position = c(0.17, 0.88),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-ADD-wWhitefish.png", width = 12, height = 7, dpi = 300)


# Statistical Analyses ----------------------------------------------------

## survival 
hatch.survival <- hatch %>% filter(eye != 0)

c <- glmer(hatch ~ temperature + group + temperature * group +        # fixed
          (1|male) + (1|female),                                      # random
          family = binomial("logit"),
          data = hatch.survival,
          control = glmerControl(optimizer = "bobyqa"))  

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1, 1)[[1]]     # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ temperature * group)
(best.emm.pair <- pairs(best.emm, simple = list("group", c("temperature")), adjust = "fdr"))


##################### incubation  ########################

hatch.ADD <- hatch %>% filter(!is.na(ADD))

c <- lmer(ADD ~ temperature + group + temperature * group +           # Fixed
         (1|male) + (1|female),                                       # Random
         data = hatch.ADD, 
         REML = FALSE)

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1, 1)[[1]]     # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post hoc test :
best.emm <- emmeans(best, ~ temperature * group)
pairs(best.emm, simple = list("group", c("temperature")), adjust = "fdr") 



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
emm_options(pbkrtest.limit = 14000)


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
hatch.survival.glm.AIC.best <- emmeans(hatch.survival.glm.AIC.best, ~ temperature * group)

## Pairwise, cld, confidence intervals
pairs(hatch.survival.glm.AIC.best, simple = list("group", c("temperature")), adjust = "bonferroni", type = "response") 
hatch.survival.glm.emm.confint <- multcomp::cld(hatch.survival.glm.AIC.best, type = "response", sort = F, adjust = "tukey") %>% 
  mutate(.group = gsub("[[:space:]]", "", .group))


# Statistical Analysis - Incubation Period (DPF) - GLM --------------------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf))

# create linear mixed model
hatch.dpf.glm <- lmer(dpf ~ temperature + group + temperature * group +       # Fixed
                      (1|male) + (1|female),                                  # Random
                      data = hatch.dpf, 
                      REML = TRUE)

# to select all model based on AICc
hatch.dpf.glm.AIC <- dredge(hatch.dpf.glm)

# select best model based on AICc
hatch.dpf.glm.AIC.best <- get.models(hatch.dpf.glm.AIC, 1)[[1]]

# ANOVA
Anova(hatch.dpf.glm.AIC.best)
Anova(hatch.dpf.glm.AIC.best, type = "III")

# post-hoc test
hatch.dpf.glm.emm <- emmeans(hatch.dpf.glm.AIC.best, ~ temperature * group)

## Pairwise, cld, confidence intervals
pairs(hatch.dpf.glm.emm, simple = list("group", c("temperature")), adjust = "bonferroni", type = "response") 
hatch.dpf.glm.emm.confint <- multcomp::cld(hatch.dpf.glm.emm, type = "response", sort = F, adjust = "tukey") %>% 
  mutate(.group = gsub("[[:space:]]", "", .group))


# Statistical Analysis - Incubation Period (ADD) - GLM --------------------

# filter to only hatched embryos
hatch.ADD <- hatch %>% filter(!is.na(ADD))

# create linear mixed model
hatch.ADD.glm <- lmer(ADD ~ temperature + group + temperature * group +       # Fixed
                        (1|male) + (1|female),                                  # Random
                      data = hatch.ADD, 
                      REML = TRUE)

# to select all model based on AICc
hatch.ADD.glm.AIC <- dredge(hatch.ADD.glm)

# select best model based on AICc
hatch.ADD.glm.AIC.best <- get.models(hatch.ADD.glm.AIC, 1)[[1]]

# ANOVA
Anova(hatch.ADD.glm.AIC.best)
Anova(hatch.ADD.glm.AIC.best, type = "III")

# post-hoc test
hatch.ADD.glm.emm <- emmeans(hatch.ADD.glm.AIC.best, ~ temperature * group)

## Pairwise, cld, confidence intervals
pairs(hatch.ADD.glm.emm, simple = list("group", c("temperature")), adjust = "bonferroni", type = "response") 
hatch.ADD.glm.emm.confint <- multcomp::cld(hatch.ADD.glm.emm, type = "response", sort = F, adjust = "tukey") %>% 
  mutate(.group = gsub("[[:space:]]", "", .group))


# Visualizations - Survival Heritability ----------------------------------

heritability.all <- bind_rows(heritability.survival, heritability.ADD, heritability.dpf) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("surv", "dpf", "add"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")),
         label = ifelse(herit == 0, 0, NA))

ggplot(heritability.all, aes(x = temperature, y = (herit * 100), group = group, fill = group)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  #geom_text(aes(x = temperature, y = 0.75, label = label), position = position_dodge(width = 0.9)) +
  scale_fill_grey("combine", start = 0.0, end = 0.8,
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
        plot.margin = unit(c(5, 5, 5, 5), 'mm')) + 
  facet_wrap(~trait, nrow = 1)

ggsave("figures/embryo/2020-Heritability.png", width = 20, height = 12, dpi = 300)


# Visualizations - Population Level ---------------------------------------

## Embryo Survival
ggplot(hatch.survival.glm.emm.confint, aes(x = temperature, y = (prob * 100), group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.22)) +
  geom_point(size = 3.25, position = position_dodge(0.22)) +
  geom_errorbar(aes(ymin = (asymp.LCL * 100), ymax = (asymp.UCL* 100)), 
                position = position_dodge(0.22),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  geom_text(aes(label = .group, y = (asymp.UCL * 100) + 2), size = 3, 
            position = position_dodge(0.22), 
            show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(10, 103), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Embryo Survival (% ± 95% CI)", color = "Populations") +
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

ggsave("figures/embryo/2020-Survival-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

## Days Post Fertilization
ggplot(hatch.dpf.glm.emm.confint, aes(x = temperature, y = emmean, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.25, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(0.15),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  geom_text(aes(label = .group, y = upper.CL+3.5), size = 3, 
            position = position_dodge(0.15), 
            show.legend = FALSE, color = "black") +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(40, 230), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (No. Days ± 95% CI)", color = "Populations") +
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

ggsave("figures/embryo/2020-DPF-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

## Accumulated Degree-Days
ggplot(hatch.ADD.glm.emm.confint, aes(x = temperature, y = emmean, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.15)) +
  geom_point(size = 3.25, position = position_dodge(0.15)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(0.15),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  geom_text(aes(label = .group, y = upper.CL+14), size = 3, 
            position = position_dodge(0.15), 
            show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.15)) +
  scale_y_continuous(limits = c(250, 850), breaks = seq(250, 850, 100), expand = c(0, 0)) +
  #scale_color_manual("combine", values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
  #                   labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(x = "Incubation Temperature (°C)", y = "Incubation Period (ADD °C ± 95% CI)", color = "Populations") +
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

ggsave("figures/embryo/2020-ADD-BW-Confint.png", width = 8.5, height = 6, dpi = 300)



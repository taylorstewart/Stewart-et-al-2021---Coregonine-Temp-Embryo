## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)


## ===========================================================
## Load incubation temperature data
## ===========================================================
ADD.2020 <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())

ADD.2019 <- read.csv("data/2019-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(temperature) %>% 
  mutate(dpf = 1:n())


## ===========================================================
## Load hatching data
## ===========================================================
hatch.USA.2020 <- read_excel("data/Artedi-Temperature-Experiments-2020.xlsx", sheet = "2020HatchingData") %>% 
  mutate(year = 2020) %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD.2020) %>% 
  dplyr::select(year, population, family, male, female, temperature, eye, premature, hatch, dpf, ADD)

hatch.USA.2019 <- read_excel("data/Artedi-Temperature-Experiments-2020.xlsx", sheet = "2019HatchingData") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  left_join(ADD.2019) %>% 
  dplyr::select(year, population, family, male, female, temperature, eye, premature, hatch, dpf, ADD)
  
hatch.Finland <- read_excel("data/Finland2019.xlsx", sheet = "L. Konnevesi vendace") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, family, male, female, temperature, eye, premature, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA.2020, hatch.Finland) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, levels = c(2, 4.5, 7, 9), ordered = TRUE))


#----------------------------------------------------------------------------------------------#
################################    Results     #############################
#----------------------------------------------------------------------------------------------#

hatch.survival.summary <- hatch %>% filter(eye != 0) %>% group_by(population, temperature) %>% 
  summarize(mean.survival = (mean(hatch)*100),
            sd.survival = (sd(hatch)*100),
            se.survival = sd.survival / sqrt(n()))

hatch.ADD.summary <- hatch %>% filter(!is.na(ADD), hatch == 1, premature != 1) %>% group_by(population, temperature) %>% 
  summarize(mean.ADD = mean(ADD),
            sd.ADD = sd(ADD),
            se.ADD = sd.ADD / sqrt(n()))

hatch.dpf.summary <- hatch %>% filter(hatch == 1, premature != 1, !is.na(dpf)) %>% group_by(population, temperature) %>% 
  summarize(mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf / sqrt(n()))


ggplot(hatch.survival.summary, aes(x = temperature, y = mean.survival, group = population, color = population, shape = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.survival-se.survival, ymax = mean.survival+se.survival), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(65, 100), breaks = seq(70, 100, 10), expand = c(0, 0)) +
  scale_color_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c("#33a02c", "#1f78b4", "#a6cee3")) +
  scale_shape_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c(17, 16, 15)) +
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

ggsave("figures/embryo/2020-Survival.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.dpf.summary, aes(x = temperature, y = mean.dpf, group = population, color = population, shape = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.dpf-se.dpf, ymax = mean.dpf+se.dpf), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(45, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c("#33a02c", "#1f78b4", "#a6cee3")) +
  scale_shape_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c(17, 16, 15)) +
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

ggsave("figures/embryo/2020-DPF.png", width = 12, height = 7, dpi = 300)


ggplot(hatch.ADD.summary, aes(x = temperature, y = mean.ADD, group = population, color = population, shape = population)) + 
  geom_line(size = 1.0, position = position_dodge(0.2)) +
  geom_point(size = 3.0, position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = mean.ADD-sd.ADD, ymax = mean.ADD+sd.ADD), size = 1.0, width = 0.3, position = position_dodge(0.2)) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(250, 950), breaks = seq(250, 950, 100), expand = c(0, 0)) +
  scale_color_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c("#33a02c", "#1f78b4", "#a6cee3")) +
  scale_shape_manual("combine", labels = c("LK-V   ", "LS-C   ", "LO-C"),
                     values = c(17, 16, 15)) +
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

ggsave("figures/embryo/2020-ADD.png", width = 12, height = 7, dpi = 300)


################################## statistique #########################

library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")


################### survival #####################

hatch.survival <- hatch %>% filter(eye != 0)

c <- glmer(hatch ~ temperature + population + temperature * population +        # fixed
          (1|male) + (1|female),                                                # random
          family = binomial("logit"),
          data = hatch.survival,
          control = glmerControl(optimizer = "bobyqa"))  

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1,"8")[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post-hoc test:
best.emm <- emmeans(best, ~ temperature * population)
(best.emm.pair <- pairs(best.emm, simple = list("population", c("temperature")), adjust = "fdr"))


##################### incubation  ########################

hatch.ADD <- hatch %>% filter(!is.na(ADD))

c <- lmer(ADD ~ population + temperature + population * temperature +
         (1|male) + (1|female), data = hatch.ADD, REML = FALSE)

dg1 <- dredge(c)                    # to select all model based on AICc
dg1

best <- get.models(dg1,"8")[[1]]    # select best model based on AICc
summary(best)

Anova(best)
Anova(best, type = "III")

# Post hoc test :
library(emmeans)
best.emm <- emmeans(best, ~ temperature * population)
pairs(best.emm, simple = list("population", c("temperature")), adjust = "fdr") 



## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(dplyr)
library(magrittr)
library(readxl)
library(ggplot2)
library(multcomp)
library(tidyr)
library(lubridate)
library(lme4)


## ===========================================================
## Load and Initial Manipulations of the Larval Hatch 
##   Number of eggs hatched
## ===========================================================
hatch.data <- read_excel("data/LakeOntario-Cisco-2019.xlsx", sheet = "HatchingData") %>% 
  dplyr::select(family, female, male, block, treatment, plate, plate_group, survival, hatch, shrinkage, mutated, rearing, dpf) %>% 
  ## remove poor survival families
  filter(family != "F01M02", family != "F03M01", family != "F03M02", family != "F03M03", family != "F03M04") %>% 
  mutate(female = factor(substr(family, 1, 3))) %>% 
  filter(hatch == 1, survival == 1)


## ===========================================================
## 
## ===========================================================
climate.temp.2_0 <- read_excel("data/LakeOntario-Cisco-IncubationTemps.xlsx", sheet = "2.0") %>% 
  mutate(date = as.POSIXct(paste0(substr(timestamp, 1, 4), "-", substr(timestamp, 6, 7), "-", substr(timestamp, 9, 10))))

climate.temp.4_5 <- read_excel("data/LakeOntario-Cisco-IncubationTemps.xlsx", sheet = "4.5") %>% 
  mutate(date = as.POSIXct(paste0(substr(timestamp, 1, 4), "-", substr(timestamp, 6, 7), "-", substr(timestamp, 9, 10))))
  
climate.temp <- bind_rows(climate.temp.2_0, climate.temp.4_5)

climate.temp.summary <- climate.temp %>% group_by(trt) %>% 
  summarize(mean.temp = mean(temp_c),
            sd.temp = sd(temp_c)) %>% ungroup()

climate.temp.2_0.daily <- climate.temp.2_0 %>% group_by(date) %>% 
  summarize(daily.incubation.temp = mean(temp_c)) %>%
  mutate(dd = cumsum(daily.incubation.temp),
         dpf = seq(0, n()-1, 1)) %>% ungroup()

climate.temp.4_5.daily <- climate.temp.4_5 %>% group_by(date) %>% 
  summarize(daily.incubation.temp = mean(temp_c)) %>%
  mutate(dd = cumsum(daily.incubation.temp),
         dpf = seq(0, n()-1, 1)) %>% ungroup()


## ===========================================================
## Calculate cumulative total, total number, and proportion for each treatment
## ===========================================================
larvae.hatch.2_0 <- hatch.data %>% 
  filter(treatment == 1) %>% 
  group_by(treatment, female, plate_group, dpf) %>% 
  summarize(quantity = sum(hatch)) %>% 
  mutate(prop = cumsum(quantity)/sum(quantity),
         trt = 2) %>% ungroup() %>% 
  dplyr::select(-treatment)

larvae.hatch.4_5 <- hatch.data %>% 
  filter(treatment == 2) %>% 
  group_by(treatment, female, plate_group, dpf) %>% 
  summarize(quantity = sum(hatch)) %>% 
  mutate(prop = cumsum(quantity)/sum(quantity),
         trt = 4.5) %>% ungroup() %>% 
  dplyr::select(-treatment)


## ===========================================================
## 
## ===========================================================
hatch.dd.2_0 <- left_join(larvae.hatch.2_0, climate.temp.2_0.daily) %>% 
  filter(female != "F03", female != "F01", female != "F09") %>% 
  mutate(plate_group = factor(plate_group))

hatch.dd.4_5 <- left_join(larvae.hatch.4_5, climate.temp.4_5.daily) %>% 
  filter(female != "F03", female != "F01", female != "F09") %>% 
  mutate(plate_group = factor(plate_group))


#lme <- glmer(prop ~ dd + female + (1|plate_group), data = hatch.dd, family = binomial)
#summary(lme)
#mc.hatch <- glht(lme, linfct = mcp(female = "Tukey"))
#summary(mc.hatch)
#cld(mc.hatch)



## Logistic regression with the 50% hatched line highlighted
ggplot(hatch.dd.2_0, aes(x = dd, y = prop, colour = female)) + 
  geom_point(show.legend = FALSE, size = 3.5) + 
  stat_smooth(aes(y = prop, group = female),  method = "glm", 
              method.args = list(family = "binomial"), se = FALSE, size = 2) +
  #geom_hline(yintercept = 0.50, colour = "#4d4d4d", size = 0.65) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(225, 410), breaks = seq(225, 400, 25), expand = c(0, 0)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#fdbf6f", "#a65628", "#f781bf", "#999999")) +
  labs(x = "Degree Days (°C)", y = "Proportion of Hatched Larvae", colour = "Female\nFamilies") +
  annotate("text", x = 236.5, y = 0.95, label = "2.0°C", size = 10) +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 30, colour = "black"),
        axis.text.x = element_text(size = 30, colour = "black"),
        axis.title.y = element_text(size = 35, margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(size = 35, margin = margin(20, 0, 0, 0)),
        legend.title = element_text(size = 30), legend.text = element_text(size = 30), 
        legend.key.size = unit(1.3, 'cm'),
        plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatching-logistic-regression-trt2_0.png", width = 17, height = 10, dpi = 300)


## Logistic regression with the 50% hatched line highlighted
ggplot(hatch.dd.4_5, aes(x = dd, y = prop, colour = female)) + 
  geom_point(show.legend = FALSE, size = 3.5) + 
  stat_smooth(aes(y = prop, group = female),  method = "glm", 
              method.args = list(family = "binomial"), se = FALSE, size = 2) +
  #geom_hline(yintercept = 0.50, colour = "#4d4d4d", size = 0.65) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(400, 700), breaks = seq(400, 700, 50), expand = c(0, 0)) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#fdbf6f", "#a65628", "#f781bf", "#999999")) +
  labs(x = "Degree Days (°C)", y = "Proportion of Hatched Larvae", colour = "Female\nFamilies") +
  annotate("text", x = 418.5, y = 0.95, label = "4.5°C", size = 10) +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 30, colour = "black"),
        axis.text.x = element_text(size = 30, colour = "black"),
        axis.title.y = element_text(size = 35, margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(size = 35, margin = margin(20, 0, 0, 0)),
        legend.title = element_text(size = 30), legend.text = element_text(size = 30), 
        legend.key.size = unit(1.3, 'cm'),
        plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatching-logistic-regression-trt4_5.png", width = 17, height = 10, dpi = 300)




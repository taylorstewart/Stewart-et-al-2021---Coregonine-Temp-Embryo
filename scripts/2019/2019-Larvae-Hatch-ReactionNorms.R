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
  mutate(female = factor(substr(family, 1, 3)),
         trt = factor(ifelse(treatment == 1, "2.0", "4.5"))) %>% 
  filter(hatch == 1, survival == 1) %>% 
  dplyr::select(-treatment)


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
larvae.hatch.2_0 <- hatch.data %>% filter(trt == "2.0")

larvae.hatch.4_5 <- hatch.data %>% filter(trt == "4.5")


## ===========================================================
## 
## ===========================================================
hatch.dd.2_0 <- left_join(larvae.hatch.2_0, climate.temp.2_0.daily) %>% 
  filter(female != "F03", female != "F01", female != "F09") %>% 
  mutate(plate_group = factor(plate_group))

hatch.dd.4_5 <- left_join(larvae.hatch.4_5, climate.temp.4_5.daily) %>% 
  filter(female != "F03", female != "F01", female != "F09") %>% 
  mutate(plate_group = factor(plate_group))

hatch.dd <- bind_rows(hatch.dd.2_0, hatch.dd.4_5)


## ===========================================================
## 
## ===========================================================
hatch.dd %<>% group_by(trt) %>% 
  summarize(n = n(),
            mean.dd = mean(dd),
            sd.dd = sd(dd),
            se.dd = sd.dd/sqrt(n),
            mean.dpf = mean(dpf),
            sd.dpf = sd(dpf),
            se.dpf = sd.dpf/sqrt(n))

data <- data.frame(lake = c(rep("Ontario", 4), rep("Superior", 4)),
                   trt = rep(c(2, 4.5, 7, 9), 2),
                   mean.dd = c(357, 493, 628, 764, 407, 518, 603, 714))

ggplot(hatch.dd, aes(x = factor(trt), y = mean.dd, group = 1)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean.dd + sd.dd, ymax = mean.dd - sd.dd), width = 0.1, size = 0.85) +
  scale_y_continuous(limits = c(325, 550), breaks = seq(350, 550, 50), expand = c(0, 0)) +
  labs(x = "Temperature (°C)", y = "Incubation Time (Degree Days)") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        text = element_text(family = 'Helvetica'), plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dd-reaction-norms.png", width = 10, height = 8, dpi = 300)



ggplot(hatch.dd, aes(x = factor(trt), y = mean.dpf, group = 1)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = mean.dpf + sd.dpf, ymax = mean.dpf - sd.dpf), width = 0.1, size = 0.85) +
  scale_y_continuous(limits = c(100, 175), breaks = seq(100, 175, 25), expand = c(0, 0)) +
  labs(x = "Temperature (°C)", y = "Incubation Time (Days Post Fertilization)") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        text = element_text(family = 'Helvetica'), plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dpf-reaction-norms.png", width = 10, height = 8, dpi = 300)



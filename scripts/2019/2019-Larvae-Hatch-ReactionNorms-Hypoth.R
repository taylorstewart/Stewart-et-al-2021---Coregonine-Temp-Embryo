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


data <- data.frame(lake = c(rep("Ontario", 4), rep("Superior", 4)),
                   trt = rep(c(2, 4.5, 7, 9), 2),
                   mean.dd = c(357, 493, 628, 764, 407, 518, 603, 714))

ggplot(data, aes(x = factor(trt), y = mean.dd, color = lake, group = lake)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  #geom_errorbar(aes(ymin = mean.dd + sd.dd, ymax = mean.dd - sd.dd), width = 0.1, size = 0.85) +
  #scale_y_continuous(limits = c(325, 550), breaks = seq(350, 550, 50), expand = c(0, 0)) +
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
  labs(x = "Temperature (°C)", y = "Incubation Time (Days)") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        text = element_text(family = 'Helvetica'), plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dpf-reaction-norms.png", width = 10, height = 8, dpi = 300)



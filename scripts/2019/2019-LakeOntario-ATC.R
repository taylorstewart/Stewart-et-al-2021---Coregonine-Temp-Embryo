#################################################################
##
##
##
#################################################################

## CLEAR THE ENVIRONMENT FIRST ==================================

rm(list = ls(all.names = TRUE))


## LOAD PACKAGES ================================================

library(dplyr)           # data manipulation
library(magrittr)        # for %<>%
library(ggplot2)         # visualization
library(lubridate)       # for easy data manipulation
library(readxl)          # reading data


## LOAD DATA ====================================================

atc.data <- read_excel("data/ATC/2019-LakeOntario-Cisco-ATC.xlsx", sheet = "Data") %>% 
  filter(include.utt == "y")


## LOAD DATA ====================================================

atc.summary <- atc.data %>% 
  filter(tank.id != "ATC 09") %>% 
  group_by(tank.id) %>% 
  summarize(n = n(),
            LD50.adm = median(utt.adm),
            LD50.temp = median(lethal.temp),
            utt.adm = max(utt.adm),
            utt.temp = max(lethal.temp))

ggplot(atc.data, aes(x = tank.id, y = lethal.temp)) +
  geom_boxplot(width = 0.5) + 
  scale_x_discrete(labels = c("Tank 1", "Tank 2", "Tank 3", "Tank 4")) + 
  scale_y_continuous(limits = c(14, 25), breaks = seq(15, 25, 2.5), expand = c(0, 0)) +
  labs(x = "Replicate", y = "Loss of Equilibrium Temperature (°C)") +
  theme(axis.text = element_text(size = 18), axis.line.x = element_line(), 
        axis.line.y = element_line(), 
        axis.title.x = element_text(size = 22, margin = margin(15, 0, 0, 0)),
        axis.title.y = element_text(size = 22, margin = margin(0, 15, 0, 0)),
        axis.ticks.length = unit(2.5, 'mm'), 
        legend.position = "none",
        panel.background = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(5, 7.5, 7.5, 7.5), "mm"))

ggsave("figures/2019-lake-ontario-4.5-atc.png", dpi = 300, width = 9, height = 7)





  

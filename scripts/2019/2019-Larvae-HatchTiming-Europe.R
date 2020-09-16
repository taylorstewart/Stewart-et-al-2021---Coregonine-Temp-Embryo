## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(dplyr)
library(readxl)
library(ggplot2)


## ===========================================================
## Load and Initial Manipulations of the Larval Hatch 
##   Number of eggs hatched
## ===========================================================
hatch.data <- read_excel("data/Europe2019-Data.xlsx", sheet = "Data") %>% 
  mutate(Lake = factor(Lake, levels = c("Konnevesi", "Constance Pelagic", "Constance Littoral", "Ontario",
                                        "Bourget", "Geneva")),
         Survival = Survival * 100)


ggplot(filter(hatch.data, Species == "lavaretus", Treatment %in% c(7, 9)), aes(x = factor(Treatment), y = ADD, color = Lake, group = Lake)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(350, 650), breaks = seq(350, 650, 50), expand = c(0, 0)) +
  scale_color_manual(values = c("grey40", "#bdc9e1", "#67a9cf", "#1c9099", "#016c59")) +
  labs(x = "Temperature (°C)", y = "Incubation Time (Degree Days)", color = "Lake Name") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks.length = unit(3, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15),
        plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dd-reaction-norms-Europe.png", width = 10, height = 6, dpi = 300)


ggplot(filter(hatch.data, Species == "lavaretus", Lake %in% c("Constance Pelagic", "Constance Littoral", "Bourget", "Geneva")),
       aes(x = factor(Treatment), y = ADD, color = Lake, group = Lake)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(400, 450), breaks = seq(400, 450, 10), expand = c(0, 0)) +
  scale_color_manual(values = c("#bdc9e1", "#67a9cf", "#1c9099", "#016c59")) +
  labs(x = "Temperature (°C)", y = "Incubation Time (Degree Days)", color = "Lake Name") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks.length = unit(3, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15),
        plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dd-reaction-norms-Europe-noFinland.png", width = 10, height = 6, dpi = 300)


ggplot(filter(hatch.data, Species %in% c("albula", "artedi")), aes(x = factor(Treatment), y = ADD, color = Lake, group = Lake)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(350, 800), breaks = seq(350, 800, 50), expand = c(0, 0)) +
  scale_color_manual(values = c("grey40", "#67a9cf")) +
  labs(x = "Temperature (°C)", y = "Incubation Time (Degree Days)", color = "Lake Name") +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks.length = unit(3, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 15), legend.title = element_text(size = 15),
        plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-hatch-dd-reaction-norms-Europe-cisco.png", width = 10, height = 6, dpi = 300)

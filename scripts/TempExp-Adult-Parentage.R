# Clear the environment first ---------------------------------------------
rm(list = ls(all.names = TRUE))

# Load packages -----------------------------------------------------------
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)

# Load parentage data -----------------------------------------------------
parentage <- read_excel("data/Coregonine-Temperature-Experiment-AdultMeasurements.xlsx", sheet = "Adult") %>% 
  mutate(Population = factor(Population, ordered = TRUE, levels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

parentage.summary <- parentage %>% group_by(Population, Sex) %>% 
  summarize(mean.TL = mean(TL),
            sd.TL = sd(TL),
            mean.FM = mean(FM),
            sd.FM = sd(FM))

# Visualization -----------------------------------------------------------
ggplot(parentage, aes(x = TL, y = FM, group = Population, color = Population)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 1000, 250), expand = c(0, 0)) +
  scale_x_continuous(limits = c(95, 550), breaks = seq(100, 550, 150), expand = c(0, 0)) +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(x = "Total Length (mm)", y = 'Weight (g)', color = "Population", shape = "Sex") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(3, 7, 5, 5), 'mm'))

ggsave("figures/adult/2020-LW.png", width = 7, height = 7, dpi = 300)


ggplot(parentage, aes(x = log(total.length.mm), y = log(weight.g), group = group, color = group)) +
  geom_point(size = 2, show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_continuous(limits = c(1.93, 7.0), breaks = seq(2, 7, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4.67, 6.3), breaks = seq(4.7, 6.3, 0.4), expand = c(0, 0)) +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#1f78b4", "#a6cee3"),
                     labels = c("LK-V   ", "LK-W   ", "LS-C   ", "LO-C")) +
  labs(x = "Log Total Length", y = 'Log Weight', color = "Population", shape = "Sex") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(3, 7, 5, 5), 'mm'))

ggsave("figures/adult/2020-LW-Log.png", width = 7, height = 7, dpi = 300)




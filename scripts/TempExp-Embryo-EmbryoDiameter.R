# Clear the environment first ---------------------------------------------
rm(list = ls(all.names = TRUE))

# Load packages -----------------------------------------------------------
library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)

# Load embryo diameter data -----------------------------------------------
data <- read_excel("data/Coregonine-Temperature-Experiment-EmbryoMeasurements.xlsx", sheet = "Data") %>% 
  mutate(female = factor(female)) %>% 
  filter(fert_success == "y")


data.summary <- data %>% group_by(population) %>% 
  summarize(mean.diam = mean(dia_mm),
            sd.diam = sd(dia_mm))

# Visualizations ----------------------------------------------------------
ggplot(data, aes(x = population, y = dia_mm, fill = population)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1.9, 2.5), breaks = seq(1.9, 2.5, 0.1), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Ontario", "Superior")) +
  scale_fill_manual(labels = c("Ontario    ", "Superior"), 
                    values = c("#a6cee3", "#1f78b4")) +
  labs(x = "", y = "Egg Diameter (mm)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/adult/2020-EmbryoDiameter.png", width = 5, height = 7, dpi = 300)

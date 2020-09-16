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
## Load embryo diameter data
## ===========================================================
data <- read_excel("data/Artedi-Temperature-Experiments-EmbryoMeasurements-2020.xlsx", sheet = "Data") %>% 
  mutate(female = factor(female))


ggplot(data, aes(x = population, y = dia_mm)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_boxplot() +
  scale_y_continuous(limits = c(1.8, 2.5), breaks = seq(1.8, 2.5, 0.1), expand = c(0, 0)) +
  scale_x_discrete(labels = c("Ontario", "Superior")) +
  labs(x = "", y = "Embryo Diameter (mm)") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/adult/2020-EmbryoDiameter.png", width = 5, height = 7, dpi = 300)


ggplot(data, aes(x = f_wt_g, y = dia_mm, group = population, color = population)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()


ggplot(data, aes(x = f_tl_mm, y = dia_mm, group = population, color = population)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw()


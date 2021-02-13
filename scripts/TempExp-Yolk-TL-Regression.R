#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(afex)
library(buildmer)
library(emmeans)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)


#### LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.lk <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LK-Larvae")
larval.ls <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LS-Larvae")
larval.lo <- read_excel("data/Coregonine-Temperature-Experiment-LarvalMeasurements.xlsx", sheet = "LO-Larvae")

# Combine each population, temperature, and species
larval <- bind_rows(larval.lk, larval.ls, larval.lo) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         block = factor(block),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Clean up environment
rm(larval.lo, larval.ls, larval.lk)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

temp <- data.frame(group = c("LK-Whitefish", "LK-Whitefish", "LK-Whitefish", "LK-Whitefish",
                             "LK-Vendace", "LK-Vendace", "LK-Vendace", "LK-Vendace",
                             "LS-Cisco", "LS-Cisco", "LS-Cisco", "LS-Cisco",
                             "LO-Cisco", "LO-Cisco", "LO-Cisco", "LO-Cisco"),
                   temperature = c(rep(c(2.2, 4.0, 6.9, 8.0),2), rep(c(2.0, 4.4, 6.9, 8.9),2)),
                   temp.treatment = factor(rep(c("Coldest ", "Cold ", "Warm ", "Warmest"), 4), 
                                           ordered = TRUE, levels = c("Coldest ", "Cold ", "Warm ", "Warmest")))

## Length-at-Hatch - Overall
larval.tl.summary <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.tl = mean(length_mm),
            se.tl = sd(length_mm)/sqrt(n())) %>% ungroup()

## Yolk-sac Volume - Overall
larval.yolk.summary <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.yolk = mean(y_vol_mm3),
            se.yolk = sd(y_vol_mm3)/sqrt(n())) %>% ungroup()

larval.all <- left_join(larval.tl.summary, larval.yolk.summary) %>% 
  left_join(temp) %>% 
  mutate(trans.tl = log(mean.tl),
         trans.yolk = log(mean.yolk))

larval.tl.max <- larval.all %>% group_by(group) %>% 
  summarize(max.tl = max(trans.tl))

lk.v.slope <- round(coef(summary(lm(trans.tl ~ trans.yolk, data = filter(larval.all, group == "LK-Vendace"))))[2,1], 2)
lk.w.slope <- round(coef(summary(lm(trans.tl ~ trans.yolk, data = filter(larval.all, group == "LK-Whitefish"))))[2,1], 2)
ls.c.slope <- round(coef(summary(lm(trans.tl ~ trans.yolk, data = filter(larval.all, group == "LS-Cisco"))))[2,1], 2)
lo.c.slope <- round(coef(summary(lm(trans.tl ~ trans.yolk, data = filter(larval.all, group == "LO-Cisco"))))[2,1], 2)

slope.data <- data.frame(group = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco"),
                         slope = c(lk.v.slope, lk.w.slope, ls.c.slope, lo.c.slope)) %>% 
  mutate(slope.comb = paste0(group, ": ", slope)) %>% 
  select(slope.comb)


ggplot(larval.all, aes(x = trans.yolk, y = trans.tl, group = group, color = group, shape = temp.treatment)) +
  geom_point(size = 3) +
  #geom_hline(data = larval.tl.max, aes(yintercept = max.tl, color = group), linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(name = "Study Groups", values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")) +
  scale_shape_manual(name = "Temperature Treatment", values = c(15, 16, 17, 3)) +
  scale_y_continuous(limits = c(1.9, 2.5), breaks = seq(1.9, 2.5, 0.1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-3, 0.5), breaks = seq(-3, 0.5, 0.5), expand = c(0, 0)) +
  #annotate(geom = "table", x = -2.9, y = 2.475, label = list(slope.data), size = 5,
  #         vjust = 1, hjust = 0, table.colnames = FALSE, table.theme = ttheme_gtminimal) +
  labs(x = "Log Yolk-sac Volume", y = "Log Total Length") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.length = unit(1.0, 'mm'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "right",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/yolk-conversion.png", width = 12, height = 7, dpi = 300)



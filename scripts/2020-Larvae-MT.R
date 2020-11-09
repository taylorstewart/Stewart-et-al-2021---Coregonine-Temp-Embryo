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


#### FILTER TO EACH TRAITS' DATASET --------------------------------------------------------------

larval.tl.na <- larval %>% filter(!is.na(length_mm), group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
larval.tl.fi <- larval %>% filter(!is.na(length_mm), group %in% c("LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()
larval.yolk.na <- larval %>% filter(!is.na(y_vol_mm3), group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
larval.yolk.fi <- larval %>% filter(!is.na(y_vol_mm3), group %in% c("LK-Vendace", "LK-Whitefish")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()

## Clean up environment
rm(larval.lo, larval.ls, larval.lk)


# STATISTICAL ANALYSIS - LENGTH-AT-HATCH - NA -------------------------------------------------

## fit full model
larval.tl.na.glm.full <- lmer(length_mm ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.tl.na)

## backward elimination to select best model
larval.tl.na.glm <- step(larval.tl.na.glm.full)
( larval.tl.na.glm.formula <- get_model(larval.tl.na.glm)@call[["formula"]])

## fit best model
larval.tl.na.glm.final <- lmer(larval.tl.na.glm.formula, data = larval.tl.na)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.na.glm.formula, data = larval.tl.na, method = "LRT")
rand(larval.tl.na.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.glm.emm <- emmeans(larval.tl.na.glm.final, ~ temperature | group)

## Pairwise
pairs(larval.tl.glm.emm, simple = list("temperature", "group"), adjust = "tukey", type = "response") 


# STATISTICAL ANALYSIS - LENGTH-AT-HATCH - FI -------------------------------------------------

## fit full model
larval.tl.fi.glm.full <- lmer(length_mm ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.tl.fi)

## backward elimination to select best model
larval.tl.fi.glm <- step(larval.tl.fi.glm.full)
( larval.tl.fi.glm.formula <- get_model(larval.tl.fi.glm)@call[["formula"]])

## fit best model
larval.tl.fi.glm.final <- lmer(larval.tl.fi.glm.formula, data = larval.tl.fi)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.fi.glm.formula, data = larval.tl.fi, method = "LRT")
rand(larval.tl.fi.glm.final)


# STATISTICAL ANALYSIS - YOLK-SAC VOLUME - NA -------------------------------------------------

## fit full model
larval.yolk.na.glm.full <- lmer(y_vol_mm3 ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.yolk.na)

## backward elimination to select best model
larval.yolk.na.glm <- step(larval.yolk.na.glm.full)
( larval.yolk.na.glm.formula <- get_model(larval.yolk.na.glm)@call[["formula"]])

## fit best model
larval.yolk.na.glm.final <- lmer(larval.yolk.na.glm.formula, data = larval.yolk.na)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.na.glm.formula, data = larval.yolk.na, method = "LRT")
rand(larval.yolk.na.glm.final)


# STATISTICAL ANALYSIS - YOLK-SAC VOLUME - FI -------------------------------------------------

## fit full model
larval.yolk.fi.glm.full <- lmer(y_vol_mm3 ~ 1 + temperature + group + temperature:group + 
                                  (1|family) + (1|male) + (1|female) + (1|block), 
                                data = larval.yolk.fi)

## backward elimination to select best model
larval.yolk.fi.glm <- step(larval.yolk.fi.glm.full)
( larval.yolk.fi.glm.formula <- get_model(larval.yolk.fi.glm)@call[["formula"]])

## fit best model
larval.yolk.fi.glm.final <- lmer(larval.yolk.fi.glm.formula, data = larval.yolk.fi)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.fi.glm.formula, data = larval.yolk.fi, method = "LRT")
rand(larval.yolk.fi.glm.final)


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

## Length-at-Hatch
larval.tl.summary <- larval %>% filter(!is.na(length_mm)) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.tl = mean(length_mm),
            se.tl = sd(length_mm)/sqrt(n()),
            n = n())

## Yolk-sac Volume
larval.yolk.summary <- larval %>% filter(!is.na(y_vol_mm3)) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.yolk = mean(y_vol_mm3),
            se.yolk = sd(y_vol_mm3)/sqrt(n()),
            n = n())


#### VISUALIZATIONS ----------------------------------------------------------

## Length-at-Hatch
plot.tl <- ggplot(larval.tl.summary, aes(x = temperature, y = mean.tl, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 3.25, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.tl - se.tl, ymax = mean.tl + se.tl), 
                position = position_dodge(0.1),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(6.5, 12), breaks = seq(7, 12, 1), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = "Mean LAH (mm)", color = "Populations") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 18, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 18, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Yolk-sac Volume
plot.yolk <- ggplot(larval.yolk.summary, aes(x = temperature, y = mean.yolk, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.13)) +
  geom_point(size = 3.25, position = position_dodge(0.13)) +
  geom_errorbar(aes(ymin = mean.yolk - se.yolk, ymax = mean.yolk + se.yolk), 
                position = position_dodge(0.13),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 1.45), breaks = seq(0, 1.4, 0.2), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = expression("Mean YSV (mm"^3*")")) +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 18, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 18, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

## Combine all figures
plot.all <- grid.arrange(arrangeGrob(textGrob(""), 
                                     get_legend(plot.tl),
                                     nrow = 1,
                                     widths = c(0.03, 1)),
                         arrangeGrob(plot.tl + theme(legend.position = "none", axis.title.x = element_blank()), 
                                     plot.yolk + theme(legend.position = "none", axis.title.x = element_blank()),
                                     nrow = 2,
                                     bottom = textGrob("Mean Incubation Temperature (°C)", x = 0.545, gp = gpar(cex = 1.5, fontfamily = "Arial"))),
                         heights = c(0.025, 1)
                         )

ggsave("figures/2020-Larval-MT-SE.png", plot = plot.all, width = 8, height = 10, dpi = 300)


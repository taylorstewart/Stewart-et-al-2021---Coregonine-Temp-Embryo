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


#### FILTER TO EACH SPECIES' DATASET --------------------------------------------------------------

## filter to only length
larval.cisco <- larval %>% filter(!is.na(length_mm), length_mm != 0, !is.na(y_vol_mm3), y_vol_mm3 != 0, group %in% c("LO-Cisco", "LS-Cisco")) %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2, 4.4, 6.9, 8.9))) %>% droplevels()
larval.vendace <- larval %>% filter(!is.na(length_mm), length_mm != 0, !is.na(y_vol_mm3), y_vol_mm3 != 0, group == "LK-Vendace") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()
larval.whitefish <- larval %>% filter(!is.na(length_mm), length_mm != 0, !is.na(y_vol_mm3), y_vol_mm3 != 0, group == "LK-Whitefish") %>% 
  mutate(temperature = factor(temperature, ordered = TRUE, levels = c(2.2, 4.0, 6.9, 8))) %>% droplevels()


## Clean up environment
rm(larval.lo, larval.ls, larval.lk)


#### STATISTICAL ANALYSIS - LENGTH-AT-HATCH - CISCO ----------------------------------------------

## fit full model
larval.tl.cisco.glm.full <- lmer(length_mm ~ 1 + temperature + group + temperature:group + 
                                (1|family) + (1|male) + (1|female) + (1|block), 
                              data = larval.cisco)

## backward elimination to select best model
larval.tl.cisco.glm <- step(larval.tl.cisco.glm.full)
( larval.tl.cisco.glm.formula <- get_model(larval.tl.cisco.glm)@call[["formula"]])

## fit best model
larval.tl.cisco.glm.final <- lmer(larval.tl.cisco.glm.formula, data = larval.cisco)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.cisco.glm.formula, data = larval.cisco, method = "LRT")
rand(larval.tl.cisco.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.cisco.glm.emm <- emmeans(larval.tl.cisco.glm.final, ~ temperature | group)

## Pairwise
pairs(larval.tl.cisco.glm.emm, simple = list("temperature", "group"), adjust = "tukey", type = "response") 


#### STATISTICAL ANALYSIS - LENGTH-AT-HATCH - VENDACE --------------------------------------------

## fit full model
larval.tl.vendace.glm.full <- lmer(length_mm ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                 data = larval.vendace)

## backward elimination to select best model
larval.tl.vendace.glm <- step(larval.tl.vendace.glm.full)
( larval.tl.vendace.glm.formula <- get_model(larval.tl.vendace.glm)@call[["formula"]])

## fit best model
larval.tl.vendace.glm.final <- lmer(larval.tl.vendace.glm.formula, data = larval.vendace)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.vendace.glm.formula, data = larval.vendace, method = "LRT")
rand(larval.tl.vendace.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.vendace.glm.emm <- emmeans(larval.tl.vendace.glm.final, ~ temperature)

## Pairwise
pairs(larval.tl.vendace.glm.emm, simple = "temperature", adjust = "tukey", type = "response") 


#### STATISTICAL ANALYSIS - LENGTH-AT-HATCH - WHITEFISH ------------------------------------------

## fit full model
larval.tl.whitefish.glm.full <- lmer(length_mm ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                   data = larval.whitefish)

## backward elimination to select best model
larval.tl.whitefish.glm <- step(larval.tl.whitefish.glm.full)
( larval.tl.whitefish.glm.formula <- get_model(larval.tl.whitefish.glm)@call[["formula"]])

## fit best model
larval.tl.whitefish.glm.final <- lmer(larval.tl.whitefish.glm.formula, data = larval.whitefish)

## likelihood ratio test for fixed and random effects
mixed(larval.tl.whitefish.glm.formula, data = larval.whitefish, method = "LRT")
rand(larval.tl.whitefish.glm.final)

## Calculate estimated marginal means - be very patient!
larval.tl.whitefish.glm.emm <- emmeans(larval.tl.whitefish.glm.final, ~ temperature)

## Pairwise
pairs(larval.tl.whitefish.glm.emm, simple = "temperature", adjust = "tukey", type = "response") 


# STATISTICAL ANALYSIS - YOLK-SAC VOLUME - CISCO ----------------------------------------------

## fit full model
larval.yolk.cisco.glm.full <- lmer(y_vol_mm3 ~ 1 + temperature + group + temperature:group + 
                                  (1|family) + (1|male) + (1|female) + (1|block), 
                                data = larval.cisco)

## backward elimination to select best model
larval.yolk.cisco.glm <- step(larval.yolk.cisco.glm.full)
( larval.yolk.cisco.glm.formula <- get_model(larval.yolk.cisco.glm)@call[["formula"]])

## fit best model
larval.yolk.cisco.glm.final <- lmer(larval.yolk.cisco.glm.formula, data = larval.cisco)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.cisco.glm.formula, data = larval.cisco, method = "LRT")
rand(larval.yolk.cisco.glm.final)


# STATISTICAL ANALYSIS - YOLK-SAC VOLUME - VENDACE --------------------------------------------

## fit full model
larval.yolk.vendace.glm.full <- lmer(y_vol_mm3 ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                   data = larval.vendace)

## backward elimination to select best model
larval.yolk.vendace.glm <- step(larval.yolk.vendace.glm.full)
( larval.yolk.vendace.glm.formula <- get_model(larval.yolk.vendace.glm)@call[["formula"]])

## fit best model
larval.yolk.vendace.glm.final <- lmer(larval.yolk.vendace.glm.formula, data = larval.vendace)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.vendace.glm.formula, data = larval.vendace, method = "LRT")
rand(larval.yolk.vendace.glm.final)

## Calculate estimated marginal means - be very patient!
larval.yolk.vendace.glm.emm <- emmeans(larval.yolk.vendace.glm.final, ~ temperature)

## Pairwise
pairs(larval.yolk.vendace.glm.emm, simple = "temperature", adjust = "tukey", type = "response") 


# STATISTICAL ANALYSIS - YOLK-SAC VOLUME - WHITEFISH --------------------------------------------

## fit full model
larval.yolk.whitefish.glm.full <- lmer(y_vol_mm3 ~ 1 + temperature + (1|family) + (1|male) + (1|female) + (1|block), 
                                     data = larval.whitefish)

## backward elimination to select best model
larval.yolk.whitefish.glm <- step(larval.yolk.whitefish.glm.full)
( larval.yolk.whitefish.glm.formula <- get_model(larval.yolk.whitefish.glm)@call[["formula"]])

## fit best model
larval.yolk.whitefish.glm.final <- lmer(larval.yolk.whitefish.glm.formula, data = larval.whitefish)

## likelihood ratio test for fixed and random effects
mixed(larval.yolk.whitefish.glm.formula, data = larval.whitefish, method = "LRT")
rand(larval.yolk.whitefish.glm.final)

## Calculate estimated marginal means - be very patient!
larval.yolk.whitefish.glm.emm <- emmeans(larval.yolk.whitefish.glm.final, ~ temperature)

## Pairwise
pairs(larval.yolk.whitefish.glm.emm, simple = "temperature", adjust = "tukey", type = "response") 


#### CALCULATE MEAN AND SE FOR NA & FI POPULATIONS -----------------------------------------------

## Length-at-Hatch
larval.tl.summary <- larval %>% filter(!is.na(length_mm), length_mm != 0) %>% 
  group_by(population, temperature, group) %>% 
  summarize(mean.tl = mean(length_mm),
            se.tl = sd(length_mm)/sqrt(n()),
            n = n())

## Yolk-sac Volume
larval.yolk.summary <- larval %>% filter(!is.na(y_vol_mm3), y_vol_mm3 != 0) %>% 
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


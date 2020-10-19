#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)
library(lme4)
library(emmeans)
library(buildmer)

emm_options(pbkrtest.limit = 14000)


# LOAD LARVAL LENGTH DATA ---------------------------------------------------------------------

larval.yolk.lo.2 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-2")
larval.yolk.ls.2 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-2")
larval.yolk.lk.2 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-2")
larval.yolk.lo.4.5 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-4.5")
larval.yolk.ls.4.5 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-4.5")
larval.yolk.lk.4.5 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-4.5")
larval.yolk.lo.7.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-7")
larval.yolk.ls.7.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-7")
larval.yolk.lk.7.0 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-7")
larval.yolk.lo.9.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LO-9")
larval.yolk.ls.9.0 <- read_excel("data/2020-Artedi-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LS-9")
larval.yolk.lk.9.0 <- read_excel("data/2019-Finland-Temperature-Experiment-LarvaeMeasurements.xlsx", sheet = "LK-9")

# Combine each population, temperature, and species
larval.yolk <- bind_rows(larval.yolk.lo.2, larval.yolk.ls.2, larval.yolk.lk.2,
                       larval.yolk.lo.4.5, larval.yolk.ls.4.5, larval.yolk.lk.4.5,
                       larval.yolk.lo.7.0, larval.yolk.ls.7.0, larval.yolk.lk.7.0,
                       larval.yolk.lo.9.0, larval.yolk.ls.9.0, larval.yolk.lk.9.0) %>% 
  filter(!is.na(y_vol_mm3)) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 4.5, 7, 9),
                              labels = c("2.0°C", "4.5°C", "7.0°C", "9.0°C")),
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
rm(larval.yolk.lo.2, larval.yolk.ls.2, larval.yolk.lk.2, larval.yolk.lo.4.5, larval.yolk.ls.4.5, larval.yolk.lk.4.5,
   larval.yolk.lo.7.0, larval.yolk.ls.7.0, larval.yolk.lk.7.0, larval.yolk.lo.9.0, larval.yolk.ls.9.0, larval.yolk.lk.9.0)


#### STATISTICAL ANALYSIS - LENGTH-AT-HATCH - GLM ------------------------------------------------

## Find best fit model - be patient!
larval.yolk.model <- buildmer(y_vol_mm3 ~ temperature + group + temperature:group + (1|male) + (1|female) + (1|family) + (1|block), 
                            direction = 'forward', data = larval.yolk, REML = TRUE)
( larval.yolk.glm <- formula(larval.yolk.model@model))
## y_vol_mm3 ~ 1 + group + temperature + group:temperature + (1|female) + (1|block)

## Create generalized linear mixed models
# fit best (full) model
larval.yolk.glm.best <- lme4::lmer(larval.yolk.glm, data = larval.yolk)

## Create generalized linear mixed models with fixed-effects removed for LRTs
# fit best model with temperature main effect removed
larval.yolk.glm.best.temp <- lme4::lmer(y_vol_mm3 ~ 1 + group + (1|female) + (1|block),
                                      data = larval.yolk)
# fit best model with population main effect removed
larval.yolk.glm.best.pop <- lme4::lmer(y_vol_mm3 ~ 1 + temperature + (1|female) + (1|block),
                                     data = larval.yolk)

## Create generalized linear mixed models with random-effects removed for LRTs
# fit best model with female random-effect removed
larval.yolk.glm.best.female <- lme4::lmer(y_vol_mm3 ~ 1 + group + temperature + group:temperature + (1|block),
                                        data = larval.yolk)
# fit best model with plate random-effect removed
larval.yolk.glm.best.block <- lme4::lmer(y_vol_mm3 ~ 1 + group + temperature + group:temperature + (1|female),
                                       data = larval.yolk)

## Calculate LRT for both temperature and population fixed-effects and interaction
# temperature:population (significant; main effects irrelevant)
drop1(larval.yolk.glm.best, test = "Chisq")
# temperature
anova(larval.yolk.glm.best.temp, larval.yolk.glm.best, test = "Chisq")
# population
anova(larval.yolk.glm.best, larval.yolk.glm.best.pop, test = "Chisq")

## Calculate LRT for female and block random-effects
# female
anova(larval.yolk.glm.best.female, larval.yolk.glm.best, test = "Chisq")
# block
anova(larval.yolk.glm.best.block, larval.yolk.glm.best, test = "Chisq")

## Calculate estimated marginal means - be very patient!
larval.yolk.glm.emm <- emmeans(larval.yolk.glm.best, ~ temperature*group)

## Pairwise, cld, confidence intervals
pairs(larval.yolk.glm.emm, simple = "temperature", adjust = "fdr", type = "response") 
larval.yolk.glm.emm.confint <- multcomp::cld(larval.yolk.glm.emm, type = "response", adjust = "fdr", sort = F, alpha = 0.05, Letters = LETTERS) %>% 
  mutate(.group = gsub("[[:space:]]", "", .group),
         lower.CL = ifelse(lower.CL < 0, 0, lower.CL))

## Save output to prevent having to re-run time consuming models
write.csv(larval.yolk.glm.emm.confint, "data/emmeans/larval_yolk_glm_emm.csv", row.names = FALSE)


#### VISUALIZATIONS ----------------------------------------------------------

larval.yolk.glm.emm.confint <- larval.yolk.glm.emm.confint %>% 
  mutate(jit.y = ifelse(group == "LK-Vendace" & temperature == "7.0°C", upper.CL-0.025, upper.CL+0.04),
         hjust = ifelse(group == "LK-Vendace" & temperature == "7.0°C", 2, 0.5))

ggplot(larval.yolk.glm.emm.confint, aes(x = temperature, y = emmean, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 1.0, position = position_dodge(0.18)) +
  geom_point(size = 3.25, position = position_dodge(0.18)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(0.18),
                size = 0.8, width = 0.2, linetype = "solid", show.legend = FALSE) +
  geom_text(aes(label = .group, y = jit.y, hjust = hjust), size = 3, 
            position = position_dodge(0.18), show.legend = FALSE) +
  scale_x_discrete(expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(-0.05, 1.8), breaks = seq(0, 1.8, 0.2), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace   ", "LK-Whitefish   ", "LS-Cisco   ", "LO-Cisco")) +
  labs(y = expression("Yolk-sac Volume (mm"^3*")"), x = 'Incubation Temperature (°C)') +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 20, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 20, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/larvae/2020-YSV-BW-Confint.png", width = 8.5, height = 6, dpi = 300)

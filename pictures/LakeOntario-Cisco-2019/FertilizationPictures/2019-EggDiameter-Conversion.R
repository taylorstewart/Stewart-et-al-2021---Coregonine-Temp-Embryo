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


## ===========================================================
## Create list of file names
## ===========================================================
measurment.files <- list.files("/Users/Taylor/CloudStation/Cisco-Climate-Change/Development-Pictures/LakeOntario-Cisco-2019/FertilizationPictures/DiamaterMeasurments",
                               full.names = TRUE)[-1]


## ===========================================================
## 
## ===========================================================
data <- do.call(bind_rows, lapply(measurment.files, function(i) {
  data <- read.csv(i, header = TRUE) %>% 
    mutate(family = substr(i, 152, 157),
           ## 6.3x magnification
           diameter_mm = (pix/10.96)/10) %>% 
    dplyr::select(family, rep, diameter_pix = pix, diameter_mm)
}))


data.survival <- read_excel("/Users/Taylor/CloudStation/Cisco-Climate-Change/Development-Pictures/LakeOntario-Cisco-DevelopmentStatus.xlsx", sheet = "SurvivalData", range = "A1:S5159") %>% 
  dplyr::select(family, female, male, block, treatment, plate_group, survival, hatch, shrinkage, mutated, rearing, dpf) %>% 
  filter(treatment == 2)

data.survival.family <- data.survival %>% group_by(family) %>% 
  summarize(egg.n = n(),
            larvae.n = sum(survival),
            survival.perc = (larvae.n/egg.n)*100) %>% 
  dplyr::select(family, survival.perc)

data.all <- left_join(data, data.survival.family) %>% 
  mutate(Survival = factor(ifelse(is.na(survival.perc), "Fertilization Failure",
                           ifelse(survival.perc < 25, "< 25%", 
                           ifelse(survival.perc >= 25 & survival.perc < 50, "25% ≥ x < 50%",
                           ifelse(survival.perc >= 50 & survival.perc < 75, "50% ≥ x < 75%", "> 75%")))),
                           ordered = TRUE, levels = c("Fertilization Failure", "< 25%", "25% ≥ x < 50%", "50% ≥ x < 75%", "> 75%"))) %>% 
  mutate(female = substr(family, 1, 3),
         female = gsub("F", "Female #", female))


ggplot(data.all, aes(x = family, y = diameter_mm)) +
  geom_boxplot(aes(fill = Survival)) + 
  scale_y_continuous(limits = c(2.18, 2.93), breaks = seq(2.2, 2.8, 0.2), expand = c(0,0)) +
  scale_fill_manual(values = c("#bd0026", "#f03b20", "#fd8d3c", "#fecc5c", "#ffffb2")) +
  labs(y = "Embryo Diameter (mm)", x = "Family") +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(), axis.ticks.length = unit(1.5, 'mm'), 
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 20, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        legend.title = element_text(size = 25), legend.text = element_text(size = 20), 
        legend.key.size = unit(1.3, 'cm'),
        strip.text = element_text(size = 20), plot.margin = unit(c(5, 5, 5, 5), "mm")) +
  facet_wrap(~ female, scales = "free_x", ncol = 3)

ggsave("/Users/Taylor/CloudStation/Cisco-Climate-Change/Development-Pictures/LakeOntario-Cisco-2019/FertilizationPictures/DiamaterMeasurments/embryo-diameter.png", width = 20, height = 12, dpi = 300)

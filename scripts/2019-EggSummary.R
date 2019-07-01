library(dplyr)
library(readxl)

data <- read_excel("data/LakeOntario-Cisco-2019.xlsx", sheet = "HatchingData") %>% 
  dplyr::select(family, female, male, block, treatment, plate_group, survival, hatch, shrinkage, mutated, rearing, dpf) %>% 
  ## remove poor survival families
  filter(family != "F01M01", family != "F01M02", family != "F01M03", family != "F01M04",
         family != "F03M01", family != "F03M02", family != "F03M03", family != "F03M04") %>% 
  mutate(family = factor(family),
         female = substr(family, 1, 3),
         female = factor(gsub("F", "Female ", female)))

data.summary <- data %>% group_by(treatment) %>% 
  summarize(egg.n = n(),
            larvae.n = sum(survival),
            survival.perc = (larvae.n/egg.n)*100,
            mortality.perc = 100-survival.perc,
            hatch.perc = (sum(ifelse(survival == 1 & hatch == 1, 1, 0))/larvae.n)*100,
            mutations.perc = (sum(mutated)/larvae.n)*100)

data.family.summary <- data %>% group_by(treatment, family) %>% 
  summarize(egg.n = n(),
            larvae.n = sum(survival),
            survival.perc = (larvae.n/egg.n)*100,
            mortality.perc = 100-survival.perc,
            hatch.perc = (sum(ifelse(survival == 1 & hatch == 1, 1, 0))/larvae.n)*100,
            mutations.perc = (sum(mutated)/larvae.n)*100)


# CLEAR THE ENVIRONMENT FIRST ---------------------------------------------

rm(list = ls(all.names = TRUE))


# LOAD PACKAGES -----------------------------------------------------------

library(dplyr)
library(readxl)
library(magrittr)
library(ggplot2)
library(nortest)
library(car)
library(lme4)
library(MuMIn)
library(ade4)
library(emmeans)

options(na.action = "na.fail")
emm_options(pbkrtest.limit = 14000)


# LOAD INCUBATION TEMPERATURE DATA ----------------------------------------

ADD.2020 <- read.csv("data/2020-Artedi-ADD.csv", header = TRUE) %>% 
  dplyr::select(population, temperature, ADD) %>% 
  group_by(population, temperature) %>% 
  mutate(dpf = 1:n())


# LOAD HATCHING DATA ------------------------------------------------------

hatch.USA.2020 <- read_excel("data/2020-Artedi-Temperature-Experiment.xlsx", sheet = "2020HatchingData") %>% 
  mutate(year = 2020) %>% 
  filter(is.na(notes) | notes != "empty well") %>% 
  filter(block != "A" | population != "superior") %>% 
  mutate(eye = as.numeric(eye),
         hatch = as.numeric(hatch)) %>% 
  filter(!is.na(eye), !is.na(hatch)) %>% 
  left_join(ADD.2020) %>% 
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.albula <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi vendace") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

hatch.Finland.lavaretus <- read_excel("data/2019-Finland-Temperature-Experiment.xlsx", sheet = "L. Konnevesi whitefish") %>% 
  mutate(year = 2019,
         premature = 0) %>% 
  dplyr::select(year, population, species, family, male, female, block, plate, temperature, eye, premature, hatch, dpf, ADD)

## Combine all populations and years
hatch <- bind_rows(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus) %>% 
  mutate(population = factor(population, levels = c("konnevesi", "superior", "ontario"), ordered = TRUE),
         temperature = factor(temperature, ordered = TRUE, 
                              levels = c(2, 4.5, 7, 9),
                              labels = c(1, 2, 3, 4)),
         female = factor(female, levels = seq(1, 12, 1),
                         labels = c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12")),
         male = factor(male, levels = seq(1, 16, 1),
                       labels = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")),
         # Create a variable with population and species combined
         group = factor(interaction(population, species), ordered = TRUE,
                        levels = c("konnevesi.albula", "konnevesi.lavaretus", "superior.artedi", "ontario.artedi"),
                        labels = c("LK-Vendace", "LK-Whitefish", "LS-Cisco", "LO-Cisco")))

## Clean up environment
rm(hatch.USA.2020, hatch.Finland.albula, hatch.Finland.lavaretus, ADD.2020)


# STATISTICAL ANALYSIS - INCUBATION PERIOD (DPF) - GLM --------------------

# filter to only hatched embryos
hatch.dpf <- hatch %>% filter(!is.na(dpf), hatch == 1)

## Check assumptions
##  Normality
ad.test(hatch.dpf$dpf)
## p = 0.122; accept normality

## Homogeneity of Variance by source
bartlett.test(dpf ~ interaction(temperature, group), data = hatch.dpf)
## p = 0.033; unequal variance!!

# create linear mixed model
hatch.dpf.glm <- lmer(dpf ~ temperature + group + temperature * group + 
                      (1|male) + (1|female) + (1|family),
                      contrasts = list(temperature = 'contr.sum', group = 'contr.sum'),
                      data = hatch.dpf)

boot.anova <- function(data, formula, conf.int = conf.int, dec = dec, reps = reps, 
                       pw.comp = pw.comp, seed = seed) {
  # Set seed (for replication)
  set.seed(seed)
  
  # Fit non-bootstrapped ANOVA to obtain Sum Sq and df
  glmFit <- lmer(dpf ~ temperature + group + temperature * group + (1|male) + (1|female) + (1|family), 
                 data = hatch.dpf, 
                 contrasts = list(temperature = 'contr.sum'))
  # Get DV name
  dvName <- colnames(model.frame(glmFit))[1]                    
  # Get IV names
  ivNames <- as.character(c(glmFit@call[["formula"]][[3]][[2]][[2]][[2]][[2]][[2]], glmFit@call[["formula"]][[3]][[2]][[2]][[2]][[2]][[3]]))
  
  # Fit ANOVA model with Type III Sum of Squares
  anovaFit <- Anova(glmFit, type = 3, test.statistic = "F", icontrasts = c("contr.sum"))
  
  # Extract F
  FValues <- anovaFit[["F"]]       # F-values
  Fdf <- anovaFit[["Df"]]                # df
  
  # Additional information for use later
  FLength <- length(FValues) - 1         # Number of F-statistics - Residual
  
  # Pairwise comparisons (if requested)
  if(is.null(pw.comp) == FALSE) {
    # Initialise variables to prevent errors
    postHocNames <- NULL
    postHocCoeff <- NULL
    postHocComparisons <- NULL
    postHocDF <- NULL
    CompName <- NULL
    
    for(i in 1:length(pw.comp)) {
      # Parses expression to be passed to lsmeans()
      exprString <- paste("pairwise ~", pw.comp[i])
      exprEval <- eval(parse(text = exprString))
      suppressMessages(postHocRes <- summary(emmeans(glmFit, exprEval)))
      
      # If interaction term specified
      if(grepl("|", pw.comp[i], fixed = TRUE) == FALSE) {
        CompName <- as.character(postHocRes$contrasts[, 1])
      } else {
        # Comparison name: "A: B"
        CompName <- paste(as.character(postHocRes$contrasts[, 1]), 
                          as.character(postHocRes$contrasts[, 2]),
                          sep = ": ")  
      }
      
      # Comparison(s) name(s)
      postHocNames <- c(postHocNames, CompName)
      
      # Get comparison coeffecient(s)
      postHocCoeff <- c(postHocCoeff, as.numeric(postHocRes$contrasts[["estimate"]]))
      
      # Get number of comparisons (to adjust confidence interval for
      # multiple comparisons)
      postHocComparisons <- c(postHocComparisons, 
                              length(as.numeric(postHocRes$contrasts[["estimate"]])))
      
      # Get df
      postHocDF <- c(postHocDF, as.numeric(postHocRes$contrasts[["df"]]))
    }
  }




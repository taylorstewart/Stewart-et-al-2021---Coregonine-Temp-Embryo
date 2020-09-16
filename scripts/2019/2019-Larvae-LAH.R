##############################################################
##############################################################
##  Cisco Dissertation - Pilot Study 2017
##
##  LARVAL LENGTH-AT-HATCH SCRIPT
##
##############################################################
##############################################################
## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))

## ===========================================================
## Load Packages -- used here
## ===========================================================
library(dplyr)       # manipulating data
library(magrittr)    # for %<>%
library(readxl)      # reading data
library(ggplot2)     # visualizations
library(multcomp)    # for glht()
library(car)         # for Boot()
library(AICcmodavg)  # for aictab()
library(lsmeans)     # calculating least-squares means
library(FSA)         # for residPlot() and fitPlot()
library(rcompanion)  # for pairwisePermutationTest()

## ===========================================================
## Set the random seed for reproducibility (i.e., randomization
##   is used in the application of the CI's).
## ===========================================================
set.seed(84621684)

## ===========================================================
## Load and Initial Manipulations of the Larval Hatch 
##   Yolk Sac Area Data
## ===========================================================
larvae.tl <- read.csv("data/2019/LakeOntario-Cisco-2019-Lengths.csv", header = TRUE) %>% 
  mutate(trt = factor(ifelse(trt == 2, "2.0", "4.5")))


## -----------------------------------------------------------
## Calculate means and std. error for each treatment
## -----------------------------------------------------------
larvae.tl.summary <- larvae.tl %>% group_by(trt) %>% 
  summarize(mean.tl = mean(length.mm),
            sd.tl = sd(length.mm),
            n = n(),
            se.tl = sd.tl/sqrt(n)) %>% 
  arrange(trt)


##############################################################
## ANALYSIS
##############################################################
## -----------------------------------------------------------
## Fit model and check assumptions
## -----------------------------------------------------------
lm.tl <- lm(length.mm ~ trt, data = larvae.tl)
summary(lm.tl)
fitPlot(lm.tl)


##------------------------------------------------------------
## Check Assumptions
##------------------------------------------------------------
##  Normality
shapiro.test(larvae.tl$length.mm)     ## p < 0.001; fail normality
residPlot(lm.tl)               ## Visually left-skewed

## Homogeneity of Variance by source
## Flinger-Killeen Test is a non-parametric test which is very robust against departures from normality.
fligner.test(length.mm ~ trt, data = larvae.tl)  ## p = 0.0428; approximately equal variance

## Kruskal-Wallis
kruskal.test(larvae.tl$length.mm, larvae.tl$trt)

pairwise.wilcox.test(larvae.tl$length.mm, larvae.tl$trt, p.adjust.method = "BH")


## create data frame with estimated marginal means and significant letters
tl.emm <- emmeans(lm.tl, ~ trt, weights = "cells")
contrast(tl.emm, "pairwise", by = NULL)

perm.pair <- pairwisePermutationTest(length.mm ~ trt, data = larvae.tl, method = 'holm')
perm.pair.ltrs <- cldList(comparison = perm.pair$Comparison, p.value = perm.pair$p.adjust, threshold  = 0.05) %>% 
  dplyr::select(group = Letter)


tl.emm.df <- data.frame(print(tl.emm)) %>% bind_cols(perm.pair.ltrs) %>% 
  mutate(n.fish = c(188, 418),
         label = paste0(trt, "\n(", n.fish, ")"))
colnames(tl.emm.df) <- c('trt', 'emmean', 'SE', 'df', 'lower.CL', 'upper.CL', 'group', 'n.fish', 'label')



##############################################################
## VISUALIZATION
##############################################################
## Using permutation pairwise test with bootstrapped estimates and confidence intervals
ggplot(tl.emm.df, aes(x = trt, y = emmean, group = 1)) +
  geom_point(size = 3) + 
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, size = 0.85) +
  annotate("text", x = 1, y = 10.72, label = "n=188", size = 6, vjust = 0.4) + 
  annotate("text", x = 2, y = 10.72, label = "n=418", size = 6, vjust = 0.4) + 
  #geom_text(aes(label = group), hjust = -1.0, vjust = 0.35, size = 10, colour = "black") + 
  scale_y_continuous(limits = c(10.7, 11.5), breaks = seq(10.7, 11.5, 0.2), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.3, 0)) +
  labs(y = "Mean Size-at-Hatch (mm ± SE)", x = 'Temperature (°C)') +
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.line = element_line(), axis.ticks.length = unit(2.0, 'mm'),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 25, margin = margin(0, 15, 0, 0)),
        axis.title.x = element_text(size = 25, margin = margin(15, 0, 0, 0)),
        text = element_text(family = 'Helvetica'), plot.margin = unit(c(8, 5, 5, 5), "mm"))

ggsave("figures/2019-LakeOntario-LAH.png", dpi = 600, width = 7, height = 8)



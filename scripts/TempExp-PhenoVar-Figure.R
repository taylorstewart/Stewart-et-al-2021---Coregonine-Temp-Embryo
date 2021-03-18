#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(egg)


# RUN EACH PHENO SCRIPT -----------------------------------------------------------------------

source("scripts/TempExp-Embryo-PhenoVar.r")
source("scripts/TempExp-Larval-PhenoVar.r")

phenoVar.all <- bind_rows(phenoVar.embryo.all, phenoVar.larval.all) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD", "LAH", "YSV"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)", "Length-at-Hatch", "Yolk-sac Volume")),
         variance = ifelse(variance == 0, 0.001, variance))


#### VISUALIZATION -------------------------------------------------------------------------------

## Create base plot
plot.phenoVar <- ggplot(mapping = aes(x = group, y = variance, group = component, fill = component)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = ifelse(variance - error < 0, 0, variance - error), 
                    ymax = ifelse(variance + error > 100, 99.5, variance + error)), 
                position = position_dodge(0.9), size = 0.2, width = 0.4, color = "gray15", show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#7bccc4", "#f0f9e8", "#bae4bc", "#2b8cbe"),
                    labels = c("Female  ", "Male  ", "Female x Male  ", "Error")) +
  #annotation_custom(textGrob("Coldest-Cold-Warm-Warmest", gp = gpar(fontsize = 15, col = "grey30")), 
  #                  xmin = 2, xmax = 2, ymin = -7.5, ymax = -7.5) +
  coord_cartesian(clip = "off") +
  labs(y = "% of Total Phenotypic Variation", x = "Study Group") + #\nTemperature Treatment") +
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_line(color = "gray90", size = 0.27), 
        panel.grid.major = element_line(color = "gray90", size = 0.27),
        panel.border = element_rect(color = "black", fill = NA, size = 0.27),
        #axis.line = element_line(size = 0.27), 
        axis.title.x = element_text(color = "Black", size = 10, margin = margin(4.6, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 10, margin = margin(0, 4.6, 0, 0)),
        axis.text.x = element_text(size = 8.2, margin = margin(0.9, 0, 0.9, 0)),
        axis.text.y = element_text(size = 8.2),
        axis.ticks.length = unit(0.9, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 9.1),
        legend.key.size = unit(0.6, 'cm'),
        legend.position = "top",
        legend.margin = margin(3, 0, -2, 0, unit = 'mm'),
        strip.text = element_text(size = 8),
        strip.background = element_rect(color = "transparent", fill = "white"),
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), 'mm')) +
  facet_wrap(~trait)

## Plot each row separately 
plot.phenoVar1 <- plot.phenoVar %+% filter(phenoVar.all, trait %in% c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")) + 
  theme(axis.title.x = element_blank())
plot.phenoVar2 <- plot.phenoVar %+% filter(phenoVar.all, trait %in% c("Length-at-Hatch", "Yolk-sac Volume")) + 
  theme(legend.position = "none")

## Combine each row
plot.phenoVar.all <- grid.arrange(grobs = lapply(
  list(plot.phenoVar1, plot.phenoVar2),
  set_panel_size,
  width = unit(5.3, "cm"),
  height = unit(5.3, "cm")
))

## Save figure
ggsave("figures/Fig4.tiff", plot = plot.phenoVar.all, width = 6.9, height = 5.7, dpi = 600)


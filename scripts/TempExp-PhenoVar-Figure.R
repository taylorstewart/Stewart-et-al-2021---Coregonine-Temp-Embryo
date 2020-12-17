#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


# RUN EACH PHENO SCRIPT -----------------------------------------------------------------------

source("scripts/TempExp-Embryo-PhenoVar.r")
source("scripts/TempExp-Larval-PhenoVar.r")

phenoVar.all <- bind_rows(phenoVar.embryo.all, phenoVar.larval.all) %>% 
  mutate(trait = factor(trait, ordered = TRUE, levels = c("survival", "dpf", "ADD", "LAH", "YSV"),
                        labels = c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)", "Length-at-Hatch", "Yolk-sac Volume")))


#### VISUALIZATION -------------------------------------------------------------------------------

plot.phenoVar <- ggplot(mapping = aes(x = group, y = variance, group = component.trt, fill = component)) + 
  geom_bar(stat = "identity", size = 0.5, position = position_dodge(0.9), color = "black") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_fill_manual(values = c("#7bccc4", "#f0f9e8", "#bae4bc", "#2b8cbe"),
                    labels = c("Female  ", "Male  ", "Female x Male  ", "Error")) +
  annotation_custom(textGrob("Coldest-Cold-Warm-Warmest", gp = gpar(fontsize = 15, col = "grey30")), 
                    xmin = 2, xmax = 2, ymin = -7.5, ymax = -7.5) +
  coord_cartesian(clip = "off") +
  labs(y = "% of Total Phenotypic Variation", x = "Study Group\nTemperature Treatment") +
  theme_bw() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16, margin = margin(2, 0, 25, 0)),
        axis.text.y = element_text(size = 16),
        axis.ticks.length = unit(2, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.25, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 16),
        strip.background = element_rect(color = "white", fill = "white"),
        plot.margin = unit(c(1, 1, 1, 1), 'mm')) +
  facet_wrap(~trait)


plot.phenoVar1 <- plot.phenoVar %+% filter(phenoVar.all, trait %in% c("Embryo Survival", "Incubation Period (DPF)", "Incubation Period (ADD)")) + 
  theme(axis.title.x = element_blank())
plot.phenoVar2 <- plot.phenoVar %+% filter(phenoVar.all, trait %in% c("Length-at-Hatch", "Yolk-sac Volume")) + 
  theme(legend.position = "none")

plot.phenoVar.all <- grid.arrange(grobs = lapply(
  list(plot.phenoVar1, plot.phenoVar2),
  set_panel_size,
  width = unit(10.5, "cm"),
  height = unit(13, "cm")
))

ggsave("figures/2020-PhenoVar.png", plot = plot.phenoVar.all, width = 14, height = 14, dpi = 300)


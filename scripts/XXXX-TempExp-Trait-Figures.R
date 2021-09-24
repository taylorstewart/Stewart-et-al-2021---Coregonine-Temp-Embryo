#### CLEAR THE ENVIRONMENT FIRST -----------------------------------------------------------------

rm(list = ls(all.names = TRUE))


#### LOAD PACKAGES -------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(egg)


# RUN EACH PHENO SCRIPT -----------------------------------------------------------------------

source("scripts/TempExp-Embryo-LHT.R")
source("scripts/TempExp-Larvae-MT.R")

traitsOverall.all <- bind_rows(hatch.survival.summary, hatch.dpf.summary, hatch.ADD.summary, larval.tl.summary, larval.yolk.summary)
traitsStand.all <- bind_rows(hatch.survival.summary.stand, hatch.dpf.summary.stand, hatch.ADD.summary.stand, larval.tl.summary.stand, larval.yolk.summary.stand)


#### VISUALIZATIONS - MEANS ----------------------------------------------------------------------

## Embryo Survival
plot.survival <- ggplot(filter(traitsOverall.all, trait == "survival"), aes(x = temperature, y = (mean.trait*100), group = group, shape = group, color = group, linetype = group)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = (mean.trait - se.trait) * 100, ymax = (mean.trait + se.trait) * 100), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean Embryo Survival (%)", x = "Temperature (°C)") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27, color = "gray95"), 
        panel.grid.major = element_line(size = 0.27, color = "gray95"),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Days Post Fertilization
plot.dpf <- ggplot(filter(traitsOverall.all, trait == "dpf"), aes(x = temperature, y = mean.trait, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(30, 225), breaks = seq(50, 225, 25), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean Days Post-fertiliztion") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27, color = "gray95"), 
        panel.grid.major = element_line(size = 0.27, color = "gray95"),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Accumulated Degree-Days
plot.add <- ggplot(filter(traitsOverall.all, trait == "ADD"), aes(x = temperature, y = mean.trait, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.9, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.13),
                size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(250, 850), breaks = seq(250, 850, 100), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean Accumulated\nDegree-days (°C)") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27, color = "gray95"), 
        panel.grid.major = element_line(size = 0.27, color = "gray95"),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Length-at-Hatch
plot.tl <- ggplot(filter(traitsOverall.all, trait == "LAH"), aes(x = temperature, y = mean.trait, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.2, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.1), size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(6.5, 12), breaks = seq(7, 12, 1), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = "Mean Length-at-Hatch (mm)", color = "Populations") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27, color = "gray95"), 
        panel.grid.major = element_line(size = 0.27, color = "gray95"),
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Yolk-sac Volume
plot.ysv <- ggplot(filter(traitsOverall.all, trait == "YSV"), aes(x = temperature, y = mean.trait, group = group, color = group, shape = group, linetype = group)) + 
  geom_line(size = 0.4, position = position_dodge(0.13)) +
  geom_point(size = 1.2, position = position_dodge(0.13), stroke = 0.6) +
  geom_errorbar(aes(ymin = mean.trait - se.trait, ymax = mean.trait + se.trait), 
                position = position_dodge(0.13), size = 0.4, linetype = "solid", show.legend = FALSE) +
  scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 1.45), breaks = seq(0, 1.4, 0.2), expand = c(0, 0)) +
  scale_color_grey("combine", start = 0.0, end = 0.8,
                   labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_shape_manual("combine", values = c(2, 5, 1, 0), 
                     labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  scale_linetype_manual("combine", values = c("solid", "dashed", "dotted", "solid"), 
                        labels = c("LK-Vendace ", "LK-Whitefish ", "LS-Cisco ", "LO-Cisco")) +
  labs(y = expression("Mean Yolk-sac Volume (mm"^3*")")) +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27, color = "gray95"), 
        panel.grid.major = element_line(size = 0.27, color = "gray95"),
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Combine all figures
plot.overall.all <- grid.arrange(
  arrangeGrob(get_legend(plot.survival)),
  arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.tl + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ysv + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 1),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob(""),
              plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()), 
              plot.add + theme(legend.position = "none", axis.title.x = element_blank()),
              textGrob(""),
              nrow = 1,
              widths = c(0.5, 1, 1, 0.5)),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob("Temperature (°C)", x = 0.5, just = "bottom", gp = gpar(cex = 0.7, fontfamily = "Arial"))),
  heights = c(0.15, 1.0, 0.05, 1.0, 0.05, 0.04)
)

ggsave("figures/summaryForDefense.tiff", plot = plot.overall.all, width = 8, height = 5, dpi = 250)


#### VISUALIZATIONS - STANDARDIZED ---------------------------------------------------------------

## Plot Standardized Survival
plot.survival.stand <- ggplot(filter(traitsStand.all, trait == "survival"), aes(x = group, y = percent.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_hline(yintercept = 0, color = "gray30", linetype = "solid", size = 0.3) +
  geom_errorbar(aes(ymin = (percent.diff - se.trait.stand), ymax = (percent.diff + se.trait.stand)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-85, 5), breaks = seq(-80.0, 0, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized Survival (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 


## Plot Standardized DPF
plot.dpf.stand <- ggplot(filter(traitsStand.all, trait == "dpf"), aes(x = group, y = percent.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_hline(yintercept = 0, color = "gray30", linetype = "solid", size = 0.3) +
  geom_errorbar(aes(ymin = (percent.diff - se.trait.stand), ymax = (percent.diff + se.trait.stand)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-75, 5), breaks = seq(-75, 0, 25), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized DPF (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 

## Plot Standardized ADD
plot.add.stand <- ggplot(filter(traitsStand.all, trait == "ADD"), aes(x = group, y = percent.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_hline(yintercept = 0, color = "gray30", linetype = "solid", size = 0.3) +
  geom_errorbar(aes(ymin = (percent.diff - se.trait.stand), ymax = (percent.diff + se.trait.stand)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-5, 105), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized ADD (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm')) 

## Plot Standardized LAH
plot.tl.stand <- ggplot(filter(traitsStand.all, trait == "LAH"), aes(x = group, y = percent.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_hline(yintercept = 0, color = "gray30", linetype = "solid", size = 0.3) +
  geom_errorbar(aes(ymin = (percent.diff - se.trait.stand), ymax = (percent.diff + se.trait.stand)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-20, 5), breaks = seq(10, -20, -5), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized LAH (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Plot Standardized YSV
plot.ysv.stand <- ggplot(filter(traitsStand.all, trait == "YSV"), aes(x = group, y = percent.diff, group = temp.treatment, fill = temp.treatment)) + 
  geom_bar(stat = "identity", size = 0.2, position = position_dodge(0.9), color = "black") +
  geom_hline(yintercept = 0, color = "gray30", linetype = "solid", size = 0.3) +
  geom_errorbar(aes(ymin = (percent.diff - se.trait.stand), ymax = (percent.diff + se.trait.stand)), 
                position = position_dodge(0.9), size = 0.3, width = 0.4, show.legend = FALSE) +
  #scale_x_continuous(limits = c(1.75, 9.15), breaks = c(2, 4, 4.4, 6.9, 8, 8.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-10, 550), breaks = seq(0, 500, 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("#0571b0", "#92c5de", "#f4a582", "#ca0020")) +
  labs(y = "Standardized YSV (%)", x = "Study Group") +
  theme_bw() +
  theme(panel.grid.minor = element_line(size = 0.27), 
        panel.grid.major = element_line(size = 0.27),
        axis.line = element_line(size = 0.1), 
        axis.title.x = element_text(color = "Black", size = 8.4, margin = margin(3.8, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 8.4, margin = margin(0, 3.8, 0, 0)),
        axis.text.x = element_text(size = 6.5, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 6.9),
        axis.ticks.length = unit(0.8, 'mm'),
        axis.ticks = element_line(size = 0.27), 
        legend.title = element_blank(),
        legend.text = element_text(size = 7.6),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "top",
        plot.margin = unit(c(1.9, 1.9, 1.9, 1.9), 'mm'))

## Combine all figures
plot.stand.all <- grid.arrange(
  arrangeGrob(get_legend(plot.survival.stand)),
  arrangeGrob(plot.survival.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              plot.tl.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              plot.ysv.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              nrow = 1),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob(""),
              plot.dpf.stand + theme(legend.position = "none", axis.title.x = element_blank()), 
              plot.add.stand + theme(legend.position = "none", axis.title.x = element_blank()),
              textGrob(""),
              nrow = 1,
              widths = c(0.5, 1, 1, 0.5)),
  arrangeGrob(textGrob("")),
  arrangeGrob(textGrob("Population", x = 0.5, just = "bottom", gp = gpar(cex = 0.7, fontfamily = "Arial"))),
  heights = c(0.15, 1.0, 0.05, 1.0, 0.05, 0.04)
)

ggsave("figures/summaryForDefense-Stand.tiff", plot = plot.stand.all, width = 8, height = 5, dpi = 250)


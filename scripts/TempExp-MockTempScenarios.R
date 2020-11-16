## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))

## ===========================================================
## Load packages
## ===========================================================
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(grid)
library(gridExtra)
library(gtable)
library(geosphere)
library(zoo)

## ===========================================================
## Load data and prepare
## ===========================================================
temp.n <- read.csv("data/Mock_Temperature/HOBO-Lake-Superior-SandIsland-ModelingMock-Normal.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"),
         group = "Normal")
temp.h <- read.csv("data/Mock_Temperature/HOBO-Lake-Superior-SandIsland-ModelingMock-Hot.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"),
         group = "Hot")

temp <- bind_rows(temp.n, temp.h) %>% 
  filter(date < as.POSIXct("2018-06-02")) %>% 
  mutate(group = factor(group, ordered = TRUE, levels = c("Normal", "Hot"))) 

temp.spawn.start <- temp %>% group_by(group) %>% 
  filter(date < "2018-01-15") %>% 
  filter(mean.sst <= 4.0) %>% 
  filter(row_number() == 1) %>% 
  dplyr::select(group, spawn.start = date)

temp.spawn.end <- temp %>% group_by(group) %>% 
  filter(date < "2018-01-30") %>% 
  filter(mean.sst <= 3.0) %>% 
  filter(row_number() == 1) %>% 
  dplyr::select(group, spawn.end = date)

temp.hatch.start <- temp %>% group_by(group) %>% 
  filter(date > "2018-03-01") %>% 
  filter(mean.sst >= 3.0) %>% 
  filter(row_number() == 1) %>% 
  dplyr::select(group, hatch.start = date)

temp.hatch.end <- temp %>% group_by(group) %>% 
  filter(date > "2018-03-01") %>% 
  filter(mean.sst >= 4.0) %>% 
  filter(row_number() == 1) %>% 
  dplyr::select(group, hatch.end = date)


temp.spawn <- left_join(temp.spawn.start, temp.spawn.end) %>% 
  group_by(group) %>% 
  mutate(mean.date = mean(as.POSIXct(c(spawn.end, spawn.start), format = c("%Y-%m-%d"))))
temp.spawn.n <- temp.spawn %>% filter(group == "Normal")
temp.spawn.h <- temp.spawn %>% filter(group == "Hot")

temp.hatch <- left_join(temp.hatch.start, temp.hatch.end) %>% 
  group_by(group) %>% 
  mutate(mean.date = mean(as.POSIXct(c(hatch.end, hatch.start), format = c("%Y-%m-%d"))))
temp.hatch.n <- temp.hatch %>% filter(group == "Normal")
temp.hatch.h <- temp.hatch %>% filter(group == "Hot")

temp.summary <- temp %>% 
  left_join(temp.spawn.end) %>% 
  left_join(temp.hatch.start) %>% 
  group_by(group) %>% 
  #filter(date > spawn.end, date < hatch.start) %>% 
  summarize(mean.temp = mean(mean.sst),
            sd.temp = sd(mean.sst))

temp.incubation <- temp %>% 
  left_join(temp.spawn.end) %>% 
  left_join(temp.hatch.start) %>% 
  group_by(group) %>% 
  filter(date > spawn.end, date < hatch.start) %>% 
  summarize(n = n())

plot <- ggplot(data = temp, aes(x = date, y = mean.sst, group = group, color = group)) +
  annotate("rect", xmin = temp.spawn.n$spawn.start, xmax = temp.spawn.n$spawn.end, ymin = 3, ymax = 4, color = NA, fill = "#67a9cf", alpha = 0.4) +
  annotate("rect", xmin = temp.spawn.h$spawn.start, xmax = temp.spawn.h$spawn.end, ymin = 3, ymax = 4, color = NA, fill = "#ef8a62", alpha = 0.4) +
  annotate("rect", xmin = temp.hatch.n$hatch.start, xmax = temp.hatch.n$hatch.end, ymin = 3, ymax = 4, color = NA, fill = "#67a9cf", alpha = 0.4) +
  annotate("rect", xmin = temp.hatch.h$hatch.start, xmax = temp.hatch.h$hatch.end, ymin = 3, ymax = 4, color = NA, fill = "#ef8a62", alpha = 0.4) +
  ## Blue Label
  geom_segment(y = 6.5, x = temp.spawn.n$mean.date, yend = 4.1, xend = temp.spawn.n$mean.date, size = 1.5, color = "#67a9cf", arrow = arrow()) +
  geom_segment(y = 6.5, x = temp.hatch.n$mean.date, yend = 4.1, xend = temp.hatch.n$mean.date, size = 1.5, color = "#67a9cf", arrow = arrow()) +
  geom_segment(y = 6.5, x = temp.spawn.n$mean.date-(15*60*60), yend = 6.5, xend = temp.spawn.n$mean.date+(24*60*60*50), size = 1.5, color = "#67a9cf") +
  geom_segment(y = 6.5, x = temp.hatch.n$mean.date+(15*60*60), yend = 6.5, xend = temp.hatch.n$mean.date-(24*60*60*50), size = 1.5, color = "#67a9cf") +
  annotate("text", x = temp.spawn.n$spawn.end+(24*60*60*71.5), y = 6.5, label = "143 Days Incubation", color = "#67a9cf", size = 5) +
  ## Orange Label
  geom_segment(y = 5.5, x = temp.spawn.h$mean.date, yend = 4.1, xend = temp.spawn.h$mean.date, size = 1.5, color = "#ef8a62", arrow = arrow()) +
  geom_segment(y = 5.5, x = temp.hatch.h$mean.date, yend = 4.1, xend = temp.hatch.h$mean.date, size = 1.5, color = "#ef8a62", arrow = arrow()) +
  geom_segment(y = 5.5, x = temp.spawn.h$mean.date-(15*60*60), yend = 5.5, xend = temp.spawn.h$mean.date+(24*60*60*10), size = 1.5, color = "#ef8a62") +
  geom_segment(y = 5.5, x = temp.hatch.h$mean.date+(15*60*60), yend = 5.5, xend = temp.hatch.h$mean.date-(24*60*60*10), size = 1.5, color = "#ef8a62") +
  annotate("text", x = temp.spawn.h$spawn.end+(24*60*60*31), y = 5.5, label = "62 Days Incubation", color = "#ef8a62", size = 5) +
  geom_line(size = 1.25) +
  geom_hline(yintercept = 4, alpha = 0.2) +
  geom_hline(yintercept = 3, alpha = 0.2) +
  scale_x_datetime(date_labels = "%m", date_breaks = "1 month", expand = c(0.0001,0.0)) +
  scale_y_continuous(limits = c(-0.25, 10), breaks = seq(0, 10, 2), expand = c(0.0, 0.0)) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62"), labels = c("2.0", "5.0")) +
  labs(y = "Water Temperature (°C)", x = "Month", color = "Mean Water Temperature (°C)") +
  theme_cowplot() +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(), 
        axis.title.y = element_text(size = 20, margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(size = 20, margin = margin(10, 0, 0, 0)),
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 17),
        axis.ticks.length = unit(2, 'mm'),
        legend.key = element_rect(fill = "white"), 
        legend.key.width = unit(3, 'lines'), 
        legend.key.height = unit(2, 'lines'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15), 
        legend.position = c(0.12, 0.85),
        plot.margin = unit(c(5, 10, 2, 2), "mm"))

ggdraw() +
  draw_image("figures/misc/egg_blue.png", x = -0.31, y = -0.16, scale = 0.1) +
  draw_image("figures/misc/egg_orange.png", x = -0.18, y = -0.059, scale = 0.07) +
  draw_image("figures/misc/larvae_blue.png", x = 0.4, y = -0.13, scale = 0.1) +
  draw_image("figures/misc/larvae_orange.png", x = 0.26, y = -0.061, scale = 0.08) +
  draw_plot(plot)

ggsave("figures/Mock_IncubationTemperature.png", height = 6, width = 9, dpi = 300)


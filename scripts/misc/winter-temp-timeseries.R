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


## ===========================================================
## Read in Superior data
## ===========================================================
temp.ls <- read_excel("/Users/Taylor/CloudStation/Cisco-Climate-Change/Modeling/coregonidR/data/LakeSuperior_SandIsland_temperature_logger_2016-18.xlsx",
                   sheet = "2016-2018") %>% 
  mutate(date = as.POSIXct(timestamp, tz = "America/New_York")) %>% 
  filter(date >= "2016-11-01" & date < "2017-06-01") %>% 
  mutate(two.weeks = round_date(date, "1 week"),
         lake = "superior") %>% 
  group_by(lake, two.weeks) %>% 
  summarize(temp.c = mean(temp.c)) %>% 
  mutate(year = year(two.weeks),
         month = month(two.weeks),
         day = day(two.weeks),
         year = ifelse(year == 2016, 2018, 2019),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d")) %>% 
  select(lake, date, temp.c)


## ===========================================================
## Read in Konnevesi data
## ===========================================================
temp.lk <- read_excel("/Users/Taylor/CloudStation/Cisco-Climate-Change/Modeling/coregonidR/data/LakeKonnevesi_temperature_2009-12.xlsx",
                      sheet = "Sheet1") %>% 
  filter(date >= "2009-10-31", date <= "2010-06-02") %>% 
  mutate(two.weeks = round_date(date, "1 week"),
         lake = "konnevesi") %>% 
  group_by(lake, two.weeks) %>% 
  summarize(temp.c = mean(temp.c)) %>% 
  mutate(year = year(two.weeks),
         month = month(two.weeks),
         day = day(two.weeks),
         year = ifelse(year == 2009, 2018, 2019),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d")) %>% 
  select(lake, date, temp.c)


## ===========================================================
## Read in Geneva data
## ===========================================================
temp.lg <- read.csv("/Users/Taylor/CloudStation/Cisco-Climate-Change/Modeling/coregonidR/data/LakeGeneva_temperature_sonde_2015_15m.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y"),
         date = sub('..', '', date),
         date = as.POSIXct(paste0("20", date), format = "%Y-%m-%d"),
         year = year(date),
         month = month(date),
         day = day(date),
         year = ifelse(year == 2014, 2018, 2019),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))
  

## ===========================================================
## Read in Bourget data
## ===========================================================
temp.lb <- read.csv("/Users/Taylor/CloudStation/Cisco-Climate-Change/Modeling/coregonidR/data/LakeBourget_temperature_sonde_2015_15m.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y"),
         date = sub('..', '', date),
         date = as.POSIXct(paste0("20", date), format = "%Y-%m-%d"),
         year = year(date),
         month = month(date),
         day = day(date),
         year = ifelse(year == 2013, 2018, 2019),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))


## ===========================================================
## Read in Constance data
## ===========================================================
temp.lc <- data.frame(lake = "constance",
                      date = as.POSIXct(c("2018-11-1", "2018-12-1", "2019-1-1", "2019-2-1", "2019-3-1", "2019-4-1", "2019-5-1", "2019-6-1"), format = "%Y-%m-%d"),
                      temp.c = c(10.42, 7.42, 6.16, 5.42, 5.6, 8.11, 10.57, 14.91))


## ===========================================================
## Spawning periods
## ===========================================================
temp.spawn.ls <- data.frame(lake = "superior",
                            date = as.POSIXct("2018-12-11"),
                            temp.c = 4.0)

temp.spawn.lk <- data.frame(lake = "konnevesi",
                            date = as.POSIXct("2018-11-10"),
                            temp.c = 3.89)

temp.spawn.lc <- data.frame(lake = "constance",
                            date = as.POSIXct("2018-12-20"),
                            temp.c = 6.65)

temp.spawn.lb <- data.frame(lake = "bourget",
                            date = as.POSIXct("2018-12-28"),
                            temp.c = 7.32)

temp.spawn.lg <- data.frame(lake = "geneva",
                            date = as.POSIXct("2019-01-10"),
                            temp.c = 8.7)

temp.spawn <- bind_rows(temp.spawn.ls, temp.spawn.lk, temp.spawn.lc, temp.spawn.lb, temp.spawn.lg)


## ===========================================================
## Hatching periods
## ===========================================================

temp.hatch.ls <- data.frame(lake = "superior",
                            date = as.POSIXct("2019-05-05"),
                            temp.c = 4.4)

temp.hatch.lk <- data.frame(lake = "konnevesi",
                            date = as.POSIXct("2019-05-06"),
                            temp.c = 4.18)

temp.hatch.lc <- data.frame(lake = "constance",
                            date = as.POSIXct("2019-03-18"),
                            temp.c = 6.7)

temp.hatch.lb <- data.frame(lake = "bourget",
                            date = as.POSIXct("2019-03-14"),
                            temp.c = 6.6)

temp.hatch.lg <- data.frame(lake = "geneva",
                            date = as.POSIXct("2019-03-10"),
                            temp.c = 6.45)

temp.hatch <- bind_rows(temp.hatch.ls, temp.hatch.lk, temp.hatch.lc, temp.hatch.lb, temp.hatch.lg)


## ===========================================================
## Combine lake temps
## ===========================================================
temp <- bind_rows(temp.ls, temp.lk, temp.lb, temp.lg, temp.lc) %>% 
  mutate(lake = factor(lake, levels = c("konnevesi", "superior", "constance", "geneva", "bourget"), ordered = TRUE))


## ===========================================================
## Plot time-series
## ===========================================================
ggplot(temp, aes(x = date, y = temp.c, group = lake)) +
  #geom_rect(aes(xmin = min(date), xmax = max(date), ymin = 2, ymax = 9), fill = "grey90") +
  geom_path(size = 1.75, aes(color = lake)) + #, linetype = lake)) + 
  geom_point(data = temp.spawn, aes(x = date, y = temp.c, group = lake), size = 5, shape = 4, stroke = 2) +
  geom_point(data = temp.hatch, aes(x = date, y = temp.c, group = lake), size = 5, shape = 1, stroke = 2) +
  scale_x_datetime(date_labels = "%m", date_breaks = "1 month", expand = c(0.0, 0.0)) +
  scale_y_continuous(limits = c(-0.25, 16), breaks = seq(0, 16, 4), expand = c(0.0, 0.0)) +
  scale_color_manual(labels = c("L. Konnevesi (FI; 62°N)  ", "L. Superior (US; 47°N)  ", "L. Constance (DE; 47°N)  ", "L. Geneva (FR; 46°N)  ", "L. Bourget (FR; 45°N)  "),
                     values = c("#a6cee3", "#1f78b4", "gray50", "#b2df8a", "#33a02c")) +
  #scale_linetype_manual(labels = c("L. Konnevesi (FI; 62°N)  ", "L. Superior (US; 47°N)  ", "L. Geneva (FR; 46°N)  ", "L. Bourget (FR; 45°N)  "),
  #                      values = c("twodash", "dashed", "solid", "dotted")) +
  labs(y = "Water Temperature (°C)", x = "Month") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(), 
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin = margin(10, 0, 0, 0)), 
        axis.title.y = element_text(size = 20, margin = margin(0, 10, 0, 0)), 
        axis.ticks.length = unit(2.5, 'mm'),
        legend.position = "top", 
        legend.text = element_text(size = 15),
        legend.title = element_blank(), 
        legend.key.size = unit(1, 'lines'),
        legend.key.width = unit(3.8, 'lines'), 
        legend.key.height = unit(2, 'lines'),
        legend.key = element_rect(fill = "white"), 
        strip.background = element_rect(fill = "white"), 
        plot.margin = unit(c(5, 5, 5, 5), "mm"))

ggsave("figures/temp-profiles.png", width = 12, height = 7, dpi = 300)


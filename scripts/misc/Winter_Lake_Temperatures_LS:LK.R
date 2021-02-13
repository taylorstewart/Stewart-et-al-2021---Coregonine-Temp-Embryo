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
temp.ls <- read_excel("data/Winter_Lake_Temperatures/LakeSuperior_SandIsland_temperature_logger_2016-18.xlsx",
                   sheet = "2016-2018") %>% 
  mutate(date = as.POSIXct(timestamp, tz = "America/New_York")) %>% 
  filter(date >= "2017-10-15" & date < "2018-06-01") %>% 
  mutate(two.weeks = round_date(date, "7 days"),
         lake = "superior") %>% 
  group_by(lake, two.weeks) %>% 
  summarize(temp.c = mean(temp.c)) %>% 
  mutate(year = year(two.weeks),
         month = month(two.weeks),
         day = day(two.weeks),
         year = ifelse(year == 2017, 2019, 2020),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d")) %>% 
  select(lake, date, temp.c)


## ===========================================================
## Read in Superior data
## ===========================================================
temp.lo <- read_excel("data/Winter_Lake_Temperatures/LakeOntario_ChaumontBay.xlsx",
                      sheet = "LakeOntario_ChaumontBay") %>% 
  mutate(date = as.POSIXct(date, tz = "America/New_York")) %>% 
  filter(date >= "2018-10-15" & date < "2019-06-01") %>% 
  mutate(two.weeks = round_date(date, "7 days"),
         lake = "ontario") %>% 
  group_by(lake, two.weeks) %>% 
  summarize(temp.c = mean(temp.c)) %>% 
  mutate(year = year(two.weeks),
         month = month(two.weeks),
         day = day(two.weeks),
         year = ifelse(year == 2018, 2019, 2020),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d")) %>% 
  select(lake, date, temp.c)


## ===========================================================
## Read in Konnevesi data
## ===========================================================
temp.lk <- read_excel("data/Winter_Lake_Temperatures/LakeKonnevesi_temperature.xlsx",
                      sheet = "Sheet1") %>% 
  filter(date >= "2017-10-15", date <= "2018-06-01") %>% 
  mutate(two.weeks = round_date(date, "7 days"),
         lake = "konnevesi") %>% 
  group_by(lake, two.weeks) %>% 
  summarize(temp.c = mean(temp.c)) %>% 
  mutate(year = year(two.weeks),
         month = month(two.weeks),
         day = day(two.weeks),
         year = ifelse(year == 2017, 2019, 2020),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d")) %>% 
  select(lake, date, temp.c)


## ===========================================================
## Read in Geneva data
## ===========================================================
temp.lg <- read.csv("/Users/taylor/SynologyDrive/Cisco-Climate-Change/Coregonine-DegreeDay-Modeling/data/LakeGeneva_temperature_sonde_2015_15m.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y"),
         date = sub('..', '', date),
         date = as.POSIXct(paste0("20", date), format = "%Y-%m-%d"),
         year = year(date),
         month = month(date),
         day = day(date),
         year = ifelse(year == 2014, 2019, 2020),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))


## ===========================================================
## Read in Bourget data
## ===========================================================
temp.lb <- read.csv("/Users/taylor/SynologyDrive/Cisco-Climate-Change/Coregonine-DegreeDay-Modeling/data/LakeBourget_temperature_sonde_2015_15m.csv", header = TRUE) %>% 
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y"),
         date = sub('..', '', date),
         date = as.POSIXct(paste0("20", date), format = "%Y-%m-%d"),
         year = year(date),
         month = month(date),
         day = day(date),
         year = ifelse(year == 2013, 2019, 2020),
         date = as.POSIXct(paste0(year, "-", month, "-", day), format = "%Y-%m-%d"))


## ===========================================================
## Read in Constance data
## ===========================================================
temp.lc <- data.frame(lake = "constance",
                      date = as.POSIXct(c("2019-11-1", "2019-12-1", "2020-1-1", "2020-2-1", "2020-3-1", "2020-4-1", "2020-5-1", "2020-6-1"), format = "%Y-%m-%d"),
                      temp.c = c(10.42, 7.42, 6.16, 5.42, 5.6, 8.11, 10.57, 14.91))


## ===========================================================
## Spawning periods
## ===========================================================
temp.spawn.ls <- data.frame(lake = "superior",
                            date = as.POSIXct("2019-12-01"),
                            temp.c = 4.0)

temp.spawn.lo <- data.frame(lake = "ontario",
                            date = as.POSIXct("2019-12-07"),
                            temp.c = 5.15)

temp.spawn.lk.a <- data.frame(lake = "konnevesi",
                            date = as.POSIXct("2019-10-27"),
                            temp.c = 5.9)

temp.spawn.lk.l <- data.frame(lake = "konnevesi",
                              date = as.POSIXct("2019-11-12"),
                              temp.c = 4.0)

temp.spawn.lc <- data.frame(lake = "constance",
                            date = as.POSIXct("2019-12-20"),
                            temp.c = 6.65)

temp.spawn.lb <- data.frame(lake = "bourget",
                            date = as.POSIXct("2019-12-28"),
                            temp.c = 7.32)

temp.spawn.lg <- data.frame(lake = "geneva",
                            date = as.POSIXct("2020-01-20"),
                            temp.c = 8)

temp.spawn <- bind_rows(temp.spawn.ls, temp.spawn.lo, temp.spawn.lk.a, temp.spawn.lk.l,
                        temp.spawn.lc, temp.spawn.lb, temp.spawn.lg)


## ===========================================================
## Hatching periods
## ===========================================================

temp.hatch.ls <- data.frame(lake = "superior",
                            date = as.POSIXct("2020-05-08"),
                            temp.c = 4.275)

temp.hatch.lo <- data.frame(lake = "ontario",
                            date = as.POSIXct("2020-05-04"),
                            temp.c = 4.22)

temp.hatch.lk <- data.frame(lake = "konnevesi",
                            date = as.POSIXct("2020-05-01"),
                            temp.c = 4.19)

temp.hatch.lc <- data.frame(lake = "constance",
                            date = as.POSIXct("2020-03-18"),
                            temp.c = 6.7)

temp.hatch.lb <- data.frame(lake = "bourget",
                            date = as.POSIXct("2020-03-14"),
                            temp.c = 6.6)

temp.hatch.lg <- data.frame(lake = "geneva",
                            date = as.POSIXct("2020-03-10"),
                            temp.c = 6.45)

temp.hatch <- bind_rows(temp.hatch.ls, temp.hatch.lo, temp.hatch.lk, 
                        temp.hatch.lc, temp.hatch.lb, temp.hatch.lg)


## ===========================================================
## Combine lake temps
## ===========================================================
temp <- bind_rows(temp.ls, temp.lo, temp.lk, temp.lc, temp.lb, temp.lg) %>% 
  mutate(lake = factor(lake, levels = c("konnevesi", "superior", "constance", "geneva", "bourget", "ontario"), ordered = TRUE))


## ===========================================================
## Plot time-series
## ===========================================================
ggplot(temp, aes(x = date, y = temp.c, group = lake)) +
  geom_path(size = 1.75, aes(color = lake)) +
  geom_point(data = temp.spawn, aes(x = date, y = temp.c, group = lake), size = 5, shape = 4, stroke = 2, show.legend = FALSE) +
  geom_point(data = temp.hatch, aes(x = date, y = temp.c, group = lake), size = 5, shape = 1, stroke = 2, show.legend = FALSE) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 4.5, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 7, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 9, linetype = "dashed", color = "gray50") +
  scale_x_datetime(date_labels = "%m", date_breaks = "1 month", expand = c(0.0, 0.0)) +
  scale_y_continuous(limits = c(-0.25, 16), breaks = seq(0, 16, 2), expand = c(0.0, 0.0)) +
  scale_color_manual(labels = c("Konnevesi (FI; 62°N)  ", "Superior (US; 47°N)  ", "Constance (DE; 47°N)", 
                                "Geneva (FR; 46°N)", "Bourget (FR; 45°N)", "Ontario (US; 44°N)"),
                     values = c("gray50", "#253494", "#2c7fb8", "#41b6c4", "#7fcdbb", "#c7e9b4")) +
  labs(y = "Water Temperature (°C)", x = "Month") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(), 
        axis.text = element_text(size = 17),
        axis.title.x = element_text(size = 20, margin = margin(10, 0, 0, 0)), 
        axis.title.y = element_text(size = 20, margin = margin(0, 10, 0, 0)), 
        axis.ticks.length = unit(2.5, 'mm'),
        legend.position = "top", 
        legend.text = element_text(size = 17),
        legend.title = element_blank(), 
        legend.key.size = unit(2, 'lines'),
        legend.key = element_rect(fill = "white"), 
        strip.background = element_rect(fill = "white"), 
        plot.margin = unit(c(2, 5, 5, 5), "mm"))

ggsave("figures/Temp-Profiles-All-AllTemps.png", width = 12, height = 7, dpi = 300)


# Clear the environment first ---------------------------------------------
rm(list = ls(all.names = TRUE))

# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(lubridate)
library(magrittr)

# Create list of data files -----------------------------------------------
sdr.files <- list.files('data/SDR_Oxygen/Hatch', recursive = TRUE, full.names = TRUE, pattern = "xlsx")

# Combine data files ------------------------------------------------------
data.oxygen <- do.call(rbind, lapply(sdr.files, function(i) {
  print(i)
  num <- which(sdr.files == i)
    
  tmp <- read_excel(i, sheet = "Oxygen Orginal", skip = 12)[c(1, 3:26)] %>% 
    mutate(num = num)
  
  n <- ceiling(nrow(tmp) * 0.10)
  
  tmp %<>% 
    ## remove first 60 minutes when temperature may be unstable
    slice(-1:-n) %>% 
    ## remove last 60 minutes when temperature may be unstable
    slice(1:(n()-n)) %>%
    ## format datetime
    rename(datetime = "Date/Time") %>% 
    mutate(datetime = gsub("[.]", "-", datetime),
           datetime = as.POSIXct(datetime, format = "%d-%m-%y %H:%M:%S", tz = "America/New_York"),
           datetime = round_date(datetime, unit = "minute"),
           ## add population from file name
           population = substr(i, 35, 36),
           ## change numerical population to long name
           population = gsub("LS", "Superior", population),
           population = gsub("LO", "Ontario", population),
           ## add temperature from file name
           temperature = substr(i, 41, 41),
           ## change numerical treatment to temperature
           temperature = as.numeric(gsub("4", "9.0", temperature)),
           temperature = as.numeric(gsub("3", "7.0", temperature)),
           temperature = as.numeric(gsub("2", "4.5", temperature)),
           temperature = as.numeric(gsub("1", "2.0", temperature)))
}))

## load survival data
data.survival <- read_excel("data/2020-Artedi-SDR-Oxygen.xlsx", sheet = "OxoDishes") %>% 
  mutate(survival = factor(survival))


## ===========================================================
## Data Manipulation - 30 Minutes
## ===========================================================
## Create new datetime to summarize by every 30 minutes
data.oxygen.30 <- data.oxygen %>% mutate(date = as.Date(datetime, format = "%Y/%m/%d", tz = "America/New_York"),
                                         hour = hour(datetime),
                                         minute = minute(datetime),
                                         minute.30 = ifelse(minute <= 30, 0, 30),
                                         datetime = as.POSIXct(paste0(date, " ", hour, ":", minute.30), format = "%Y-%m-%d %H:%M")) %>% 
  select(datetime, population, temperature, num, 2:25)

## gather well columns into tidy format
data.oxygen.30 %<>% pivot_longer(5:28, "well")

## Summarize by every 30 minutes
data.oxygen.30 %<>% group_by(datetime, population, temperature, num, well) %>% 
  summarize(value = mean(value))


## ===========================================================
## Data Manipulation - 1 Minute
## ===========================================================
## gather well columns into tidy format
data.oxygen.1 <- data.oxygen %>% pivot_longer(2:25, "well")


## ===========================================================
## 
## ===========================================================
data <- left_join(data.oxygen.30, data.survival)

data.hatched <- data %>% filter(datetime <= hatch.date) %>% 
  group_by(population, temperature, well) %>% 
  summarize(n = n(),
            mean.DO = mean(value)) %>% 
  group_by(population, temperature) %>% 
  summarize(weight.mean.DO = weighted.mean(mean.DO, n),
            sd.DO = sd(mean.DO),
            se.DO = sd.DO/sqrt(n())) %>% 
  mutate(temperature = factor(temperature))


## ===========================================================
## Visualization
## ===========================================================
ggplot(data.hatched, aes(x = temperature, y = weight.mean.DO, group = population, fill = population)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), size = 0.5, color = "black") +
  geom_errorbar(aes(ymin = weight.mean.DO-se.DO, ymax = weight.mean.DO+se.DO),
                width = 0.3, position = position_dodge(0.9), size = 0.7) +
  scale_y_continuous(limits = c(0, 12.4), breaks = seq(8.2, 12.4, 0.6), expand = c(0, 0)) +
  scale_fill_manual(labels = c("Ontario    ", "Superior"), 
                    values = c("#a6cee3", "#1f78b4")) +
  coord_cartesian(ylim = c(8.0, 12.4)) +
  labs(x = "Incubation Temperature (°C)", y = "Dissolved Oxygen (mg/L ± SE)") +
  theme_classic() +
  theme(axis.title.x = element_text(color = "Black", size = 22, margin = margin(10, 0, 0, 0)),
        axis.title.y = element_text(color = "Black", size = 22, margin = margin(0, 10, 0, 0)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18),
        axis.ticks.length = unit(1.5, "mm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1.0, 'cm'),
        legend.position = "top",
        strip.text = element_text(size = 15),
        plot.margin = unit(c(5, 5, 5, 5), 'mm'))

ggsave("figures/embryo/2020-SDR-Oxygen.png", width = 12, height = 7, dpi = 300)



hatch.survival <- hatch %>% filter(eye != 0)
hatch.filt <- filter(hatch.survival, group == "LO-Cisco", temperature == "9.0°C")

s <- length(levels(hatch.filt$female))
d <- length(levels(hatch.filt$male))
k <- hatch.filt %>% group_by(family) %>% 
  summarize(n = n()) %>% ungroup() %>% 
  summarize(mean.n = mean(n)) %>% pull(mean.n)

herit.surv.model <- mmer(hatch ~ 1,
             random = ~ male + female + male:female + block + plate,
             rcov = ~ units,
             data = hatch.filt)

herit.surv.model <- lmer(hatch ~ 1 + (1|male) + (1|female) + (1|male:female) + (1|block) + (1|plate),
                         data = data.temp.group, REML = F)


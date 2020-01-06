# # Repeatability Analyses # # 

# load packages and dataset
library(tidyverse)
library(rptR)
library(lme4)

# # using rptR to analyze data # #

# checking out what effects are relevent first
l1  <- lmer(startNum ~ year + (1|cgPlaId), data = phen_19967) 
l1b <- lmer(startNum ~ year + expNm + (1|cgPlaId), data = phen_19967)
summary(l1b)
anova(l1, l1b) # Hmmm seems like we should adjust for year planted....

l2  <- lmer(dur ~ year + (1|cgPlaId), data = phen_19967) # year matters with duration - why is this? AIC is lower
l2b <- lmer(dur ~ year + expNm + (1|cgPlaId), data = phen_19967)
anova(l2, l2b)


l3 <- lmer(dur ~ 1 +    (1|cgPlaId), data = phen_19967)
summary(l2)
summary(l3)
anova(l2, l3)
# conclusion: use year to analyze both...

# repeatability of start time - including year as fixed effect

# are the data normal? Pretty much....
phen_19967 %>% ggplot()+
  geom_histogram(aes(startNum))+
  facet_wrap(~year)
# Models
r1 <- rpt(startNum ~ (1|cgPlaId), grname = "cgPlaId", data = phen_19967, datatype = "Gaussian",
          nboot = 0)

# Now adjusting for year to year variation as well as expNm
r1 <- rpt(startNum ~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 0, npermut = 0)
summary(r1)

# adding in row and position to see if it changes anything
r1rp <- rpt(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
            data = phen_19967, datatype = "Gaussian",
            nboot = 0, npermut = 0)
summary(r1rp)

# testing differences without using rptR
l4  <- lmer(startNum ~ year + yrPlanted + (1|cgPlaId), data = phen_19967) 
l4r <- lmer(startNum ~ year + yrPlanted + row + (1|cgPlaId), data = phen_19967)
l4p <- lmer(startNum ~ year + yrPlanted + pos + (1|cgPlaId), data = phen_19967)
l4b <- lmer(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), data = phen_19967)
anova(l4b)

anova(l4, l4b) # p = 0.014
anova(l4, l4r) # row doesn't matter
anova(l4, l4p) # position does matter (should include position....)

# try adding site as a random effect instead of cgPlaId - phenology not as repeatable for site
phen_1996 %>% 
  group_by(site, year) %>%
  summarize(count = n()) %>%
  ggplot()+
  geom_col(aes(site, count), color = "black", fill = "white") +
  labs(x = NULL)+
  facet_wrap(~year)+
  theme(axis.text.x = element_text(size = rel(1.1), angle = 45))

# fitting with site only
r2 <- rpt(startNum ~ year + (1|site), grname = "site", data = phen_1996, datatype  = "Gaussian",
          nboot = 100) # why is this a singular fit?
summary(r2)

# fitting with cgPlaID AND site
r3 <- rpt(startNum ~ year + (1|cgPlaId) + (1|site), grname = c("cgPlaId", "site"), data = phen_1996, datatype = "Gaussian", nboot = 100)
summary(r3)

# what about repeatability of headct - Poisson distribution?

r3 <- rpt(headCt ~ (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Poisson", 
          nboot = 500)
summary(r3)
print(r3)

# repeatability of duration - sort of normally distributed, although left skewed
phen_1996 %>%
  ggplot(aes(year, dur))+
  geom_count(alpha = 0.2)+
  stat_summary(fun.y = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)+
  labs(x = NULL, 
       y = "Flowering duration (days)",
       size = "Count")+
  theme(axis.text.x = element_text(angle = 45, size = rel(1.1), hjust = 1),
        axis.title.y = element_text(size = rel(1.15)),
        legend.background = element_rect(color = "black"))

r4 <- rpt(dur ~ (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Gaussian", nboot = 1000)

r5 <- rpt(dur ~ year + (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Gaussian", nboot = 1000)
summary(r5)

r6 <- rpt(dur ~ (1|site), grname = "site", data = phen_1996, datatype = "Gaussian", nboot = 100)
summary(r6)

phen_1996 %>%
  ggplot(aes(headCt, dur))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), colour = "red")+
  theme_bw()
cor(phen_1996$headCt, phen_1996$dur)

# # repeatability of DAM - should be the same as adjusted repeatability # #

rr1 <- rpt(dam ~ 1 + (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Gaussian",
           nboot = 500)
summary(rr1) # pretty close










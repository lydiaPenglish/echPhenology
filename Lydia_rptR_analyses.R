# Script to conduct repeatability analysis uses rptR package on FFD and flowering duration

# load packages and dataset
library(tidyverse)
library(rptR)
library(lme4)
data("phen_dataset")
my_cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
             "#A6761D","#666666")

# # using rptR to analyze data # #

# 1. FFD

# checking out what factors should be included as fixed effects (aka what should we be adjusting for)
l1  <- lmer(startNum ~ 1 + (1|cgPlaId), data = phen_19967) 
l1a <- lmer(startNum ~ year + (1|cgPlaId), data = phen_19967)
l1b <- lmer(startNum ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l1,  l1a) # Model fits WAY better with year
anova(l1a, l1b) # surprisingly model fits much better with yrPlanted

# adding in row and position
l1c <- lmer(startNum ~ year + yrPlanted + row + (1|cgPlaId), data = phen_19967)
l1d <- lmer(startNum ~ year + yrPlanted + pos + (1|cgPlaId), data = phen_19967)
l1e <- lmer(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), data = phen_19967)
anova(l1c, l1b) # row doesn't matter
anova(l1d, l1b) # position does matter
# same results via backwards elimination (vs addition)

# are the FFD data normal? Pretty much....
phen_19967 %>% ggplot()+
  geom_histogram(aes(startNum))+
  facet_wrap(~year)

# repeatbility model with everything adjusted
r1 <- rpt(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE)
summary(r1) 
# repeatability model without adjusted for row/pos
r1b <- rpt(startNum ~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r1b)

# 2. Duration

# checking out what factors should be included as fixed effects (aka what should we be adjusting for)
l2  <- lmer(dur ~ 1 + (1|cgPlaId), data = phen_19967) 
l2a <- lmer(dur ~ year + (1|cgPlaId), data = phen_19967)
l2b <- lmer(dur ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l2,  l2a) # Model fits better with year
anova(l2a, l2b) # Model also fits better with yrPlanted

# adding in row and position
l2c <- lmer(dur ~ year + yrPlanted + row + (1|cgPlaId), data = phen_19967)
l2d <- lmer(dur ~ year + yrPlanted + pos + (1|cgPlaId), data = phen_19967)
l2e <- lmer(dur ~ year + yrPlanted + row + pos + (1|cgPlaId), data = phen_19967)
anova(l2c, l2b) # row doesn't matter
anova(l2d, l2b) # position also doesn't matter
# same results via backwards elimination (vs addition)

# are the duration data normal? Meh, yeah looks ok...
phen_19967 %>% ggplot()+
  geom_histogram(aes(dur))+
  facet_wrap(~year)

# repeatbility model with everything adjusted
r2 <- rpt(dur ~ year + yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000)
summary(r2) 
# repeatability model without adjusted for row/pos
r2b <- rpt(dur~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r2b)



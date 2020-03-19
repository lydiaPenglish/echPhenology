# Script to conduct repeatability analysis uses rptR package on FFD and flowering duration

# load packages and dataset
library(tidyverse)
library(rptR)
library(lme4)
data("phen_dataset")

# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# 1. FFD

# # examining fix factors to adjust our repeatability estimates # # 

# base model 
l1  <- lmer(startNum ~ 1 + (1|cgPlaId), data = phen_19967) 

# add in year
l1a <- lmer(startNum ~ year + (1|cgPlaId), data = phen_19967)
anova(l1a, l1)    # keep!

# add in cohort
l1b <- lmer(startNum ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l1a, l1b)   # keep!

# add in headCt 
l1c <- lmer(startNum ~ year + yrPlanted + headCt + (1|cgPlaId), data = phen_19967)
anova(l1b, l1c)   # keep

# in locs
l1d <- lmer(startNum ~ year + yrPlanted + headCt + row + (1|cgPlaId), data = phen_19967)
l1e <- lmer(startNum ~ year + yrPlanted + headCt + pos + (1|cgPlaId), data = phen_19967)
l1f <- lmer(startNum ~ year + yrPlanted + headCt + row + pos + (1|cgPlaId), data = phen_19967)

anova(l1d, l1c) # row does't matter
anova(l1e, l1c) # position matters
anova(l1c, l1f)

# use l1f as the final model (keep row even though doesn't sig effect model fit)
summary(l1f)

# check model
lmerTest::rand(l1f) # random effect matters
performance::r2(l1f)
performance::check_model(l1f)    # looks fine, a little bit of a weird tail

# are the FFD data normal? Pretty much....
phen_19967 %>% ggplot()+
  geom_histogram(aes(startNum))+
  facet_wrap(~year)

# repeatbility model with everything adjusted
r1 <- rpt(startNum ~ year + yrPlanted + headCt + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE)
summary(r1) 
print(r1)
# repeatability model without adjusted for row/pos
r1b <- rpt(startNum ~ year + yrPlanted + headCt + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r1b)

# same model but using dam instead of start num
r1c <- rpt(dam ~ yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r1c)
 
 
# 2. Duration

# checking out what factors should be included as fixed effects (aka what should we be adjusting for)
# base model
l2  <- lmer(dur ~ 1 + (1|cgPlaId), data = phen_19967) 

# adding year
l2a <- lmer(dur ~ year + (1|cgPlaId), data = phen_19967)
anova(l2, l2a)      # keep!

# adding cohort
l2b <- lmer(dur ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l2a, l2b)     # keep!

# adding headCt
l2c <- lmer(dur ~ year + yrPlanted + headCt + (1|cgPlaId), data = phen_19967)
anova(l2b, l2c)

# adding in row and position
l2d <- lmer(dur ~ year + yrPlanted + headCt + row + (1|cgPlaId), data = phen_19967)
l2e <- lmer(dur ~ year + yrPlanted + headCt + pos + (1|cgPlaId), data = phen_19967)
l2f <- lmer(dur ~ year + yrPlanted + headCt + row + pos + (1|cgPlaId), data = phen_19967)

anova(l2c, l2d)  # ns
anova(l2c, l2e)  # ns
anova(l2c, l2f)


# are the duration data normal? Meh, yeah looks ok...
phen_19967 %>% ggplot()+
  geom_histogram(aes(dur))+
  facet_wrap(~year)

# repeatbility model with everything adjusted
r2 <- rpt(dur ~ year + yrPlanted + headCt + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000)
summary(r2) 

# repeatability model without adjusted for row/pos
r2b <- rpt(dur~ year + yrPlanted + headCt + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r2b)


# Script to conduct repeatability analysis uses rptR package on FFD and flowering duration

# load packages and dataset
library(dplyr)
library(rptR)
library(lme4)
library(ggplot2)
data("phen_dataset")

# FFD: fitting mixed models to get correct fixed effects ------------------------------

# graph of distribution
phen_19967 %>%
  filter(year != "2016" & year !="2017") %>%
  ggplot(aes(startNum)) +
  geom_histogram(binwidth = 2)+
  facet_wrap(~year)

# # examining fix factors to adjust our repeatability estimates # # 

# base model 
l1  <- lmer(startNum ~ 1 + (1|cgPlaId), data = phen_19967) 
summary(l1)

# add in year
l1a <- lmer(startNum ~ year + (1|cgPlaId), data = phen_19967)
summary(l1a)
anova(l1a, l1)    # keep!
ggResidpanel::resid_panel(l1a)
performance::check_model(l1a)


# add in cohort, by itself and with year
l1b <-  lmer(startNum ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
l1b2 <- lmer(startNum ~ yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l1a, l1b)   # keep!
anova(l1, l1b2)   # keep!

# add in headCt - don't do this bc phenotype! 
# l1c <-  lmer(startNum ~ year + yrPlanted + headCt + (1|cgPlaId), data = phen_19967)
# l1c2 <- lmer(startNum ~ headCt + year + (1|cgPlaId), data = phen_19967)
# anova(l1b, l1c)   # keep
# summary(l1c2)
# performance::compare_performance(l1, l1c2)

# in locs
l1d <- lmer(startNum ~ year + yrPlanted + row + (1|cgPlaId), data = phen_19967)
l1e <- lmer(startNum ~ year + yrPlanted + pos + (1|cgPlaId), data = phen_19967)
l1f <- lmer(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), data = phen_19967)

anova(l1d, l1b) # row does't matter
anova(l1e, l1b) # position matters
anova(l1b, l1f)

# use l1f as the final model (keep row even though doesn't sig effect model fit)
summary(l1f)

# check model
lmerTest::rand(l1f) # random effect matters
performance::r2(l1f)
performance::check_model(l1f)    # looks fine, a little bit of a weird tail
ggResidpanel::resid_panel(l1f)

# FFD: repeatability model --------------------------------------------------------------

# repeatbility model with everything adjusted
r1 <- rpt(startNum ~ year + yrPlanted +  row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE)
summary(r1) 
print(r1)

# repeatability model with fixed effects in denominator of repeatbility estimator
r2 <- rpt(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE,
          adjusted = FALSE)
summary(r2)

r3 <- rpt(startNum ~ year + yrPlanted + row + pos + (1|cgPlaId), 
          grname = c("cgPlaId", "Fixed", "Residual"),
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE,
          ratio = FALSE)
summary(r3)


# repeatability model without adjusted for row/pos
r1b <- rpt(startNum ~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r1b)

# same model but using dam instead of start num
r1c <- rpt(dam ~ yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r1c)
 

# FFD: adding site ------------------------------------------------

# site as another random effect?
s1 <- rpt(startNum ~ year + yrPlanted +  row + pos + (1|site), grname = "site",
          data = phen_19967, datatype = "Gaussian",
          nboot = 100, npermut = 100, parallel = TRUE)
summary(s1)

# site AND cgPlaID as random effects?
s2 <- rpt(startNum ~ year + yrPlanted +  row + pos + (1|cgPlaId) + (1|site), 
          grname = c("cgPlaId","site"),
          data = phen_19967, datatype = "Gaussian",
          nboot = 100, npermut = 100, parallel = TRUE)
summary(s2)


# Duration: fitting models to get fixed effect adjustments -------------------------

# checking out what factors should be included as fixed effects (aka what should we be adjusting for)
# base model
l2  <- lmer(dur ~ 1 + (1|cgPlaId), data = phen_19967) 

# adding year
l2a <- lmer(dur ~ year + (1|cgPlaId), data = phen_19967)
anova(l2, l2a)      # keep!

# adding cohort
l2b  <- lmer(dur ~ year + yrPlanted + (1|cgPlaId), data = phen_19967)
l2b2 <- lmer(dur ~ yrPlanted + (1|cgPlaId), data = phen_19967)
anova(l2b, l2a)     # keep!
anova(l2, l2b2)     # keep!

# adding in row and position
l2d <- lmer(dur ~ year + yrPlanted +  row + (1|cgPlaId), data = phen_19967)
l2e <- lmer(dur ~ year + yrPlanted +  pos + (1|cgPlaId), data = phen_19967)
l2f <- lmer(dur ~ year + yrPlanted +  row + pos + (1|cgPlaId), data = phen_19967)

anova(l2b, l2d)  # ns
anova(l2b, l2e)  # ns
anova(l2b, l2f)  # ns - get rid of...?

# are the duration data normal? Meh, yeah looks ok...
phen_19967 %>%
  filter(year != "2016" & year !="2017") %>%
  ggplot()+
  geom_histogram(aes(dur))+
  facet_wrap(~year)

# checking model:
ggResidpanel::resid_panel(l2f)
performance::check_model(l2f)

g2f <- glmer(dur ~ year + yrPlanted + (1|cgPlaId), data = phen_19967,
             family = poisson(link = "log"))   # doesn't converge row/pos included
performance::check_overdispersion(g2f)
performance::check_convergence(g2f)
performance::check_model(g2f)
ggResidpanel::resid_panel(g2f)

# Duration: Repeatability models -----------------------------------------------------

# repeatbility model with everything adjusted
r2 <- rpt(dur ~ year + yrPlanted + row + pos + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 1000, npermut = 1000, parallel = TRUE)
summary(r2) 

# repeatability model without adjusted for row/pos
r2b <- rpt(dur~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
           data = phen_19967, datatype = "Gaussian",
           nboot = 1000, npermut = 1000)
summary(r2b)

# try poisson repeatability without row/pos - takes a really long time and says it won't converge
r3 <- rpt(dur ~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Poisson", link = ("log"),
          nboot = 0, npermut = 0, parallel = TRUE)
summary(r3)

# Head count: fitting models ---------------------------------------

# Does head count vary by year? Kind of....
phen_19967 %>%
  ggplot(aes(year, headCt))+
  geom_count(alpha = 0.5)+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)

h0 <- glmer(headCt ~ 1 + (1|cgPlaId), 
            data = phen_19967, family = poisson(link = "log"))

h1 <- glmer(headCt ~ yrPlanted + (1|cgPlaId), 
            data = phen_19967, family = poisson(link = "log"))
anova(h0, h1)  # ok keep year planted
performance::check_overdispersion(h1)

# adding year, putting in optimizer from here (https://stats.stackexchange.com/questions/164457/r-glmer-warnings-model-fails-to-converge-model-is-nearly-unidentifiable)
h2 <- glmer(headCt ~ yrPlanted + year + (1|cgPlaId), 
            data = phen_19967, family = poisson(link = "log"),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(h2)
anova(h2, h1)

r4 <- rpt(headCt ~ yrPlanted + year + (1|cgPlaId), 
          grname = "cgPlaId",
          data = phen_19967, datatype = "Poisson", 
          nboot = 10, npermut = 10, parallel = TRUE)
summary(r4)

# Is head count related to start date? ---------------------------------------

# look first:
phen_19967 %>%
  ggplot(aes(headCt, startNum))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~year)

mod_fun <- function(df){
  lm(startNum ~ headCt, data = df, family = poisson(link = "log"))
}

mods <- phen_19967 %>%
  group_by(year) %>%
  tidyr::nest() %>%
  mutate(model  = purrr::map(data, mod_fun),
         tidy   = purrr::map(model, broom::tidy),
         glance = purrr::map(model, broom::glance)) %>%
  tidyr::unnest(tidy) %>%
  tidyr::unnest(glance, names_repair = "universal") %>%
  filter(term == "headCt") %>%
  mutate_if(is.numeric, ~round(., 3))


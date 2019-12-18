# # Repeatability Analyses # # 

# load packages and dataset
library(tidyverse)
library(rptR)
library(lme4)

# # cor between start date and end date # # 

cor(phen_1996$startNum, phen_1996$endNum) # 0.87
phen_1996 %>%
  ggplot(aes(startNum, endNum))+
  geom_point()

# # correlatios between (standardized) FFDs of different years

cor_tab <- phen_19967 %>%
  select(cgPlaId, year, dam) %>%
  arrange(cgPlaId) %>%
  pivot_wider(names_from = year, values_from = dam) %>%
  column_to_rownames("cgPlaId")
cormat <- round(cor(cor_tab, use = "pairwise.complete.obs", method = "spearman"), 2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
up_cormat <- get_upper_tri(cormat)
up_melt <- reshape2::melt(up_cormat, na.rm = TRUE)

corP <- round(cor.test(cor_tab, use = "pairwise.complete.obs", method = "spearman"), 2)

# # using rptR to analyze data # #

# checking out what effects are relevent first
l1  <- lmer(startNum ~ year + (1|cgPlaId), data = phen_19967) 
l1b <- lmer(startNum ~ year + expNm + (1|cgPlaId), data = phen_19967)
summary(l1b)
anova(l1, l1b) # Hmmm seems like we should adjust for year planted....

l2  <- lmer(dur ~ year + (1|cgPlaId), data = phen_19967) # year matters with duration - why is this? AIC is lower
l2b <- lmer(dur ~ year + expNm + (1|cgPlaId), data = phen_19967)
anova(l2, l2b)


l3 <- lmer(dur ~ 1 +    (1|cgPlaId), data = phen_1996)
summary(l2)
summary(l3)
anova(l1, l1b)
# conclusion: use year to analyze both...

# repeatability of start time - including year as fixed effect

# are the data normal? Pretty much....
phen_19967 %>% ggplot()+
  geom_histogram(aes(startNum))+
  facet_wrap(~year)
# Models
r1 <- rpt(startNum ~ (1|cgPlaId), grname = "cgPlaId", data = phen_19967, datatype = "Gaussian",
          nboot = 0)
# Another version of the same model where you get raw variances for the fixed effects and residuals
r1 <- rpt(startNum ~ year + yrPlanted + (1|cgPlaId), grname = "cgPlaId",
          data = phen_19967, datatype = "Gaussian",
          nboot = 500)
summary(r1)

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


# ---- adding columns for consistency analysis -- previous research methods ----

# DAM = day after median. Subtracts the median first flowering date of all plants in a year from the 
# first flowering date of an individual. Examples: a dam of -5 means a plant flowered 5 days before
# the median first day, and a dam of 0 means the plant started flowering on the first day. 

# DAF = day after first. 

damrank <- p1996 %>%
  group_by(year)%>%
  summarize(medianDate = median(startDtEarly)) %>%
  left_join(p1996, ., by = "year")%>%
  dplyr::mutate(dam = as.numeric(startDtEarly - medianDate))%>%
  dplyr::mutate(dur = endNum - startNum)%>%
  # rank: ranks FFDs (doesn't care if days occur in between)
  group_by(year) %>%
  dplyr::mutate(rank = dense_rank(startDtEarly))%>%
  dplyr::mutate(daf = startDtEarly - min(startDtEarly)) %>%
  # year as factor
  ungroup(year)%>%
  dplyr::mutate(year = forcats::as_factor(as.character(year)))%>%
  dplyr::mutate(site = forcats::as_factor(site))










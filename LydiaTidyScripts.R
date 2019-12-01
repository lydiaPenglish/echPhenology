# Script to tidy Echinacea phenology data and save as an .rda file or use in a markdown doc
#       created by:    Lydia English
#       last modified: 11/24/2019 by LE

# # Load packages # # 

library(tidyverse) # using tidy syntax
library(lubridate) # to help with dates
library(lme4)      # mixed models
library(rptR)      # repeatability package
theme_set(theme_bw())

# # Load phenology datasets # #

p.05.17x <- read_csv("data-raw/exPt1Phenology.csv")                 # data from 2005 - 2017
p.10x <- read_delim("data-raw/2010.CG1.Phenology.txt", delim = ",") # 2010 data
ped <- read_csv("data-raw/96979899qGenPedigreeLE.csv")              # pedigree data for site origins
rowpos <- read_csv("data-raw/1996rowPosData.csv")                   # row/position data

# # Data wrangling # # 

# 2010 data 
p.10 <- p.10x %>%
  dplyr::rename("year" = "Year")%>%
  dplyr::rename("startDtEarly"= "startDateEarly")%>%
  dplyr::rename("startDtLate" = "startDateLate")%>%
  dplyr::rename("endDtEarly"  = "endDateEarly")%>%
  dplyr::rename("endDtLate"   = "endDateLate")%>%
  # converting dates from characters to dates 
  mutate_at(vars(startDtEarly:endDtLate), lubridate::ymd_hms)%>%
  mutate_at(vars(startDtEarly:endDtLate), lubridate::date)%>%
  select(-note)

# 2005 - 2017 data
p.05.17 <- p.05.17x %>%
  dplyr::rename("year" = "phenYear")%>%
  group_by(year, cgPlaId) %>%
  # getting number of heads, deleting cgHdId
  dplyr::mutate(headCt = n()) %>%
  ungroup() %>%
  select(-c(row, pos, ttColor, phenNote, cgHdId))

# need one row per plant - consolidating flowering phenology to earliest and latest dates
#    if a plant has multiple heads

# earliest starts
starts <- p.05.17 %>%
  group_by(year, cgPlaId)%>%
  filter(startDtEarly == min(startDtEarly))%>%
  select(year, cgPlaId, startDtEarly, headCt)
# latest ends
ends <-   p.05.17 %>%
  group_by(year, cgPlaId)%>%
  filter(endDtLate    == max(endDtLate)) %>%
  select(year, cgPlaId, endDtLate, headCt)

# join starts and ends back together
p.05.17f <- left_join(starts, ends, by = c("year", "cgPlaId", "headCt")) %>%
  arrange(year, cgPlaId) %>%
  dplyr::distinct()%>%
  # get rid of bad dates
  filter(startDtEarly > 1940-01-01)%>%
  filter(endDtLate    > 1940-01-01)

# adding back in 2010 data
p.all <- bind_rows(p.05.17f, p.10)%>%
  arrange(year, cgPlaId) %>%
  select(-c(endDtEarly, startDtLate))

# reference pop - 1996, 1997, 1998, 1999
ALL1990s <- p.all %>% filter(cgPlaId < 2500) 

# 1996 dataset
x1996 <-p.all %>% filter(cgPlaId < 647) %>%
            # add columns with julian date
          dplyr::mutate(startNum = yday(startDtEarly))%>%
          dplyr::mutate(endNum   = yday(endDtLate))%>%
            # add column with the number of times something flowered
          group_by(cgPlaId) %>%
          dplyr::mutate(phenCt = n()) %>%
          ungroup() %>%
          filter(phenCt > 1)

# pedigree data 
ped96 <- ped %>%
  filter(cgplaid < 647) %>%
  # rename columns
  dplyr::rename("site" = "siteOfOriginPedigree")%>%
  dplyr::rename("cgPlaId" = "cgplaid")%>%
  # recode levels
  mutate(site = recode(site, 
                       "AA" = "aa",
                       "Eriley" = "eri",
                       "Lf" = "lf",
                       "NWLF" = "nwlf",
                       "SPP" = "spp",
                       "Stevens" = "sap",
                       "Nessman" = "ness"
                       ))

# merging data
phen_1996 <- left_join(x1996, ped96, by = "cgPlaId")%>%
  # get rid of empty levels
  filter(site != "",
         site != "Unknown",
         site != "Tower",
         cgPlaId != 0) %>%
  # add in row position data
         left_join(., rowpos, by = "cgPlaId")%>%
  ungroup(year)%>%
  mutate(dur = endNum - startNum)%>%
  mutate(year = forcats::as_factor(year)) %>%
  group_by(year)%>%
  mutate(medianDate = median(startDtEarly)) %>%
  mutate(dam = as.numeric(startDtEarly - medianDate)) %>%
  dplyr::mutate(rank = dense_rank(startDtEarly))%>%
  dplyr::mutate(daf = startDtEarly - min(startDtEarly)) %>%
  ungroup(year) %>%
  dplyr::mutate(site = forcats::as_factor(site))

# # saving data for easy access in the future # #  
save(phen_1996, file = "phen_dataset.rda")

# # trying out a plot where I look at plants that flowered more than 7 times and their spread of flowering times
phen_1996 %>%
  filter(phenCt > 7) %>%
  group_by(cgPlaId) %>%
  summarize(meanDAM = mean(dam),
            sdDAM   = sd(dam),
            n       = n(),
            seDAM   = sdDAM/sqrt(n))%>%
  mutate(newID = rank(meanDAM)) %>%
  ggplot(aes(newID, meanDAM))+
  geom_point()+
  geom_errorbar(aes(ymin=meanDAM-seDAM, ymax=meanDAM+seDAM), colour="black", width=.1)

phen_7 <- phen_1996 %>%
  filter(phenCt > 7) %>%
  group_by(cgPlaId) %>%
  mutate(medDAM = median(dam)) %>%
  ungroup() %>%
  distinct(cgPlaId, .keep_all = TRUE) %>%
  arrange(desc(medDAM)) %>%
  mutate(newID = row_number()) %>%
  select(cgPlaId, medDAM, newID) %>%
  left_join(., phen_1996, by = "cgPlaId")

phen_7 %>% 
  ggplot(aes(x = newID, y = dam)) +
  geom_boxplot(aes(group = newID))+
  #geom_point(alpha = 0.2)+
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip()

# # looking at correlations

# cor between start date and end date

cor(phen_1996$startNum, phen_1996$endNum)
phen_1996 %>%
  ggplot(aes(startNum, endNum))+
  geom_point()

# # using rptR to analyze data # #

# checking out what effects are relevent first
l1 <- lmer(startNum ~ year + (1|cgPlaId), data = p1996) 
l2 <- lmer(dur ~ year + (1|cgPlaId), data = phen_1996) # year matters with duration - why is this? AIC is lower
l3 <- lmer(dur ~ 1 +    (1|cgPlaId), data = phen_1996)
summary(l2)
summary(l3)
anova(l2, l3)
# conclusion: use year to analyze both...

# repeatability of start time - including year as fixed effect
r1 <- rpt(startNum ~ (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Gaussian",
                  nboot = 1000)
# Another version of the same model where you get raw variances for the fixed effect and residuals
r1 <- rpt(startNum ~ year + (1|cgPlaId), grname = c("cgPlaId", "Fixed", "Residual", "Overdispersion"),
          data = phen_1996, datatype = "Gaussian",
                  nboot = 500, npermut = 500, 
          parallel = TRUE, ratio = FALSE)
summary(r1)

# try adding site as a random effect instead of cgPlaId - phenology not as repeatable for site
r2 <- rpt(startNum ~ 1 + (1|site), grname = ("site"), data = phen_1996, datatype  = "Gaussian",
          nboot = 0) # why is this a singular fit?
summary(r2)

# what about repeatability of headct - Poisson distribution

r3 <- rpt(headCt ~ (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Poisson", 
          nboot = 500)
summary(r3)
print(r3)

# repeatability of duration
phen_1996 %>%
  ggplot(aes(dur))+
  geom_bar(color = "black", fill = "white")+
  xlab("Flowering duration (days)")+
  theme_bw()

r4 <- rpt(dur ~ (1|cgPlaId), grname = "cgPlaId", data = p1996, datatype = "Gaussian", nboot = 1000)
r5 <- rpt(dur ~ year + (1|cgPlaId), grname = "cgPlad", data = p1996, datatype = "Gaussian", nboot = 1000)
summary(r4)
plot(r4)

g1 <- glm(headCt ~ site, family = "poisson", data = damrank)
summary(g1)

p1996 %>%
  ggplot(aes(site, headCt))+
  geom_boxplot(aes(group = site))

phen_1996 %>%
  ggplot(aes(headCt, dur))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), colour = "red")+
  theme_bw()
cor(phen_1996$headCt, phen_1996$dur)

g2 <- glm(headCt ~ dur, family = "poisson", data = p1996)
summary(g2)

g3 <- lm(headCt ~ poly(dur, 4), data = p1996)
summary(g3)

# - repeatability of DAM - should be the same as adjusted repeatability

rr1 <- rpt(dam ~ 1 + (1|cgPlaId), grname = "cgPlaId", data = phen_1996, datatype = "Gaussian",
           nboot = 500)
summary(rr1)

# # Correlation between DAMS in years # # 

cor_tab <- phen_1996 %>%
  select(cgPlaId, year, dam) %>%
  pivot_wider(names_from = year, values_from = dam) %>%
  column_to_rownames("cgPlaId")
cor(cor_tab, use = "pairwise.complete.obs", method = "spearman")
cor.test(cor_tab$`2009`, cor_tab$`2005`)

# -- adding columns for consistency analysis -- previous research methods-- #

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

# average DAM
# exclude plants that just flowered once
DR_twoPlus <- damrank %>%
  filter(phenCt >=2)
DR_twoPlus %>%
  group_by(cgPlaId)%>%
  summarize(meanDAM  = mean(dam),
            varDAM   = var(dam),
            phenCt   = n())

# -- correlations -- #

# 1. Between start and end time

DAMS %>%
  ggplot(aes(startNum, endNum))+
  geom_point()+
  geom_smooth(method = lm)
cor.test(~ endNum + startNum, data = damrank) 

# 2. Between duration and head count
damrank %>%
  ggplot(aes(headCt, dur))+
  geom_boxplot(aes(group = headCt))
# with headCt as continuous
p1 <- lm(dur ~ poly(headCt, 2), data = damrank) # fits better with polynomial
summary(p1)

# -- helpful plots -- #

# 1. Differences in flowering time between years

p1996 %>%
  ggplot(aes(year, startNum))+
  geom_point()+
  ylab("First Flowering Date (julian)")+
  xlab("Year")+
  scale_x_continuous(breaks = c(2005,2007,2009, 2011, 2013, 2015))

# 2. Day after medians

damrank %>%
  ggplot(aes(year, dam))+
  geom_boxplot(aes(group = year))

# -- using rptR to analyze data -- #

# checking out what effects are relevent first
l1 <- lmer(startNum ~ year + (1|cgPlaId), data = p1996)
summary(l1)


r1 <- rpt(startNum ~ (1|cgPlaId), grname = "cgPlaId", data = damrank, datatype = "Gaussian",
                  nboot = 1000)

r1 <- rpt(startNum ~ year + (1|cgPlaId), grname = "cgPlaId", data = damrank, datatype = "Gaussian",
                  nboot = 1000)
r1
summary(r1)
# try and get Katherine's package to work
# ggResidpanel::resid_panel(r1)

# try adding site as a random effect - phenology not repeatable for site
r2 <- rpt(startNum ~ year + ( 1+ year|site), grname = "site", data = damrank, datatype  = "Gaussian")
print(r2)

# what about repeatability of headct - Poisson distribution

r3 <- rpt(headCt ~ (1|cgPlaId), grname = "cgPlaId", data = damrank, datatype = "Poisson", nboot = 0)
print(r3)

# repeatability of duration
hist(damrank$dur)

r4 <- rpt(dur ~ (1|cgPlaId), grname = "cgPlaId", data = damrank, datatype = "Gaussian", nboot = 1000)
summary(r4)
plot(r4)

g1 <- glm(headCt ~ site, family = "poisson", data = damrank)
summary(g1)
print(g1)
damrank %>%
  ggplot(aes(site, headCt))+
  geom_boxplot(aes(group = site))


### site comparisons - tidied

library(tidyverse)
#library(echPhenology)
library(lubridate)
library(lme4)
library(rptR)

# -- read in datasets -- #

# using dataset of 2005 - 2015, just need to add 2010


p.05.17x <- read_csv("exPt1Phenology.csv")
p.10x <- read_delim("2010.CG1.Phenology.txt", delim = ",") # 2010 data
ped <- read_csv("96979899qGenPedigreeLE.csv")
rowpos <- read_csv("1996rowPosData.csv")


# -- data wrangling -- # 

p.10 <- p.10x %>%
  dplyr::rename("year" = "Year")%>%
  dplyr::rename("startDtEarly"= "startDateEarly")%>%
  dplyr::rename("startDtLate" = "startDateLate")%>%
  dplyr::rename("endDtEarly"  = "endDateEarly")%>%
  dplyr::rename("endDtLate"   = "endDateLate")%>%
  # converting characters into dates
  mutate_at(vars(startDtEarly:endDtLate), lubridate::ymd_hms)%>%
  mutate_at(vars(startDtEarly:endDtLate), lubridate::date)%>%
  select(-note)
p.05.17 <- p.05.17x %>%
  dplyr::rename("year" = "phenYear")%>%
  group_by(year, cgPlaId) %>%
  dplyr::mutate(headCt = n()) %>%
  select(-c(row, pos, ttColor, phenNote, cgHdId))

# -- one row per plant:taking full range of flowering time if plant has multiple heads --# 

starts <- p.05.17 %>%
  group_by(year, cgPlaId)%>%
  filter(startDtEarly == min(startDtEarly))%>%
  select(year, cgPlaId, startDtEarly, headCt)
ends <-   p.05.17 %>%
  group_by(year, cgPlaId)%>%
  filter(endDtLate    == max(endDtLate)) %>%
  select(year, cgPlaId, endDtLate, headCt)

# join back togeter
p.05.17f <- left_join(starts, ends) %>%
  dplyr::distinct()%>%
  arrange(year, cgPlaId) %>%
  # get rid of bad dates
  filter(startDtEarly > 1940-01-01)%>%
  filter(endDtLate    > 1940-01-01)

# adding back in 2010 data
p.all <- bind_rows(p.05.17f, p.10)%>%
  arrange(year, cgPlaId) %>%
  select(-c(endDtEarly, startDtLate))

# reference pop
ALL1990s <- p.all %>% filter(cgPlaId < 2500) 

# 1996 dataset
x1996 <-    p.all %>% filter(cgPlaId < 647) %>%
            # add columns with julian date
            dplyr::mutate(startNum = yday(startDtEarly))%>%
            dplyr::mutate(endNum   = yday(endDtLate))%>%
            # add column with the number of times something flowered
            group_by(cgPlaId) %>%
            dplyr::mutate(phenCt = n())

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
p1996 <- left_join(x1996, ped96, by = "cgPlaId")%>%
  # get rid of empty levels
  filter(site != "",
         site != "Unknown",
         site != "Tower",
         cgPlaId != 0) %>%
  # add in row position data
         left_join(., rowpos, by = "cgPlaId")%>%
  ungroup(year)%>%
  mutate(dur = endNum - startNum)%>%
  mutate(year = forcats::as_factor(year))

# -- saving data for easy access in the future -- # 

save(p1996, file = "phen_dataset.rda")


# -- using rptR to analyze data -- #

# checking out what effects are relevent first
l1 <- lmer(startNum ~ year + (1|cgPlaId), data = p1996) 
l2 <- lmer(dur ~ year + (1|cgPlaId), data = p1996)
summary(l2)
# year matters

# repeatability of start time - including year as fixed effect
r1 <- rpt(startNum ~ (1|cgPlaId), grname = "cgPlaId", data = damrank, datatype = "Gaussian",
                  nboot = 1000)

r1 <- rpt(startNum ~ year + (1|cgPlaId), grname = "cgPlaId", data = p1996, datatype = "Gaussian",
                  nboot = 1000)
summary(r1)

# try adding site as a random effect instead of cgPlaId - phenology not as repeatable for site
r2 <- rpt(startNum ~ year + (1|site), grname = ("site"), data = p1996, datatype  = "Gaussian",
          nboot = 1000)
summary(r2)

# what about repeatability of headct - Poisson distribution

r3 <- rpt(headCt ~ (1|cgPlaId), grname = "cgPlaId", data = p1996, datatype = "Poisson", nboot = 0)
print(r3)

# repeatability of duration
p1996 %>%
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

p1996 %>%
  ggplot(aes(headCt, dur))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ poly(x, 3, raw = TRUE), colour = "red")+
  theme_bw()

g2 <- glm(headCt ~ dur, family = "poisson", data = p1996)
summary(g2)

g3 <- lm(headCt ~ poly(dur, 4), data = p1996)
summary(g3)

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
  group_by(year) %>%
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


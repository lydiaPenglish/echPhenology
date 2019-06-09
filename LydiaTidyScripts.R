### site comparisons - tidied

library(tidyverse)
library(echPhenology)
library(lubridate)

# -- read in datasets -- #

# using this dataset of 2005 - 2015, just need to add 2010

p.05.15 <- read_csv("exPt1Phenology20151130.csv")
ped <- read_csv("96979899qGenPedigreeLE.csv")
rowpos <- read_csv("1996rowPosData.csv")


p13x <- read_csv("2013.cg1.phenology.csv") # 2013 data phenology data
p.05.12x <- read_delim("allYearsNot2010.CG1.Phenology.txt", delim = ",") # 2005-2012 phenology data, no 2010
p.10x <- read_delim("2010.CG1.Phenology.txt", delim = ",") # 2010 data


# all data?
allP <- read_csv("exPt1Phenology20151130.csv")

# -- data wrangling -- # 

p.05.12 <- p.05.12x %>%
  dplyr::rename("cgHdId" = "cgHeadId") %>%
  # making column for number of heads
  group_by(year, cgPlaId) %>%
  mutate(headCt = n()) %>%
  # converting dates to POSIXt
  mutate_at(vars(startDateEarly:endDateLate), as.character) %>%
  mutate_at(vars(startDateEarly:endDateLate), lubridate::ymd)

p.10 <- p.10x %>%
  dplyr::rename("year" = "Year")%>%
  # converting characters into dates
  mutate_at(vars(startDateEarly:endDateLate), lubridate::ymd_hms)%>%
  mutate_at(vars(startDateEarly:endDateLate), lubridate::date)

p13 <- p13x %>%
  dplyr::rename("year" = "phenYear") %>%
  dplyr::rename("startDateEarly"= "startDtEarly") %>%
  dplyr::rename("startDateLate" = "startDtLate") %>%
  dplyr::rename("endDateEarly"  = "endDtEarly") %>%
  dplyr::rename("endDateLate"   = "endDtLate") %>%
  # making column for number of flowering heads
  group_by(year, cgPlaId) %>%
  mutate(headCt = n())%>%
  # converting date-time to just date
  mutate_at(vars(startDateEarly:endDateLate), lubridate::date)

# merging three datasets together:

p.temp <- bind_rows(p.05.12, p.10) %>%
          bind_rows(., p13)

# -- one row per plant:taking full range of flowering time if plant has multiple heads --# 

# filter by early and late dates
starts <- p.temp %>%
  group_by(year, cgPlaId)%>%
  filter(startDateEarly == min(startDateEarly))%>%
  select(year, cgPlaId, startDateEarly, headCt)
ends <- p.temp %>%
  group_by(year, cgPlaId)%>%
  filter(endDateLate    == max(endDateLate)) %>%
  select(year, cgPlaId, endDateLate, headCt)

# join back togeter
p.05.13 <- left_join(starts, ends) %>%
  dplyr::distinct()%>%
  arrange(year, cgPlaId) %>%
  # get rid of bad dates
  filter(startDateEarly > 1940-01-01)%>%
  filter(endDateLate    > 1940-01-01)

# reference pop
ALL1990s <- p.05.13 %>% filter(cgPlaId < 2500) 

# 1996
x1996 <-    p.05.13 %>% filter(cgPlaId < 647)





# checking to make sure this dataset is accurate
tester <- allP %>%
  filter(phenYear < 2014)%>%
  rename("year" = "phenYear")%>%
  rename("startDateEarly" = "startDtEarly")%>%
  rename("endDateLate" = "endDtLate")%>%
  group_by(year, cgPlaId) %>%
  mutate(headCt = n())%>%
  select(year, cgPlaId, startDateEarly, endDateLate, headCt)%>%
  arrange(year, cgPlaId)
testStart <- tester %>%
  group_by(year, cgPlaId)%>%
  filter(startDateEarly == min(startDateEarly))%>%
  select(year, cgPlaId, startDateEarly, headCt)
testEnd <- tester %>%
  group_by(year, cgPlaId)%>%
  filter(endDateLate    == max(endDateLate)) %>%
  select(year, cgPlaId, endDateLate, headCt)
tester2 <- left_join(testStart, testEnd)%>%
  dplyr::distinct()%>%
  arrange(year, cgPlaId) %>%
  # get rid of bad dates
  filter(startDateEarly > 1940-01-01)%>%
  filter(endDateLate    > 1940-01-01)

baddies <- anti_join(p.05.13, tester2) # looks good - just need 2010 data

# Formating dates to fit with the phenology package
p.plant <- updateDT(p.plant, "startDtEarly", "sd.E")
p.plant <- updateDT(p.plant, "startDtLate", "sd.L")
p.plant <- updateDT(p.plant, "endDtEarly", "ed.E")
p.plant <- updateDT(p.plant, "endDtLate", "ed.L")


# reference pop
ALL1990s <- p.plant[p.plant$cgPlaId < 2500, ] 

# 1996
x1996 <- p.plant[p.plant$cgPlaId < 647, ]

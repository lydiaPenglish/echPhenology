### site comparisons - tidied

library(tidyverse)
#library(echPhenology)
library(lubridate)

# -- read in datasets -- #

# using dataset of 2005 - 2015, just need to add 2010

p.05.15x <- read_csv("exPt1Phenology20151130.csv")
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
p.05.15 <- p.05.15x %>%
  dplyr::rename("year" = "phenYear")%>%
  group_by(year, cgPlaId) %>%
  mutate(headCt = n()) %>%
  select(-c(row, pos, ttColor, phenNote, cgHdId))

# -- one row per plant:taking full range of flowering time if plant has multiple heads --# 

starts <- p.05.15 %>%
  group_by(year, cgPlaId)%>%
  filter(startDtEarly == min(startDtEarly))%>%
  select(year, cgPlaId, startDtEarly, headCt)
ends <-   p.05.15 %>%
  group_by(year, cgPlaId)%>%
  filter(endDtLate    == max(endDtLate)) %>%
  select(year, cgPlaId, endDtLate, headCt)

# join back togeter
p.05.15f <- left_join(starts, ends) %>%
  dplyr::distinct()%>%
  arrange(year, cgPlaId) %>%
  # get rid of bad dates
  filter(startDtEarly > 1940-01-01)%>%
  filter(endDtLate    > 1940-01-01)

# adding back in 2010 data
p.all <- bind_rows(p.05.15f, p.10)%>%
  arrange(year, cgPlaId) %>%
  select(-c(endDtEarly, startDtLate))

# reference pop
ALL1990s <- p.all %>% filter(cgPlaId < 2500) 

# 1996 dataset
x1996 <-    p.all %>% filter(cgPlaId < 647) %>%
            # add columns with julian date
            mutate(startNum = yday(startDtEarly))%>%
            mutate(endNum   = yday(endDtLate))%>%
            # add column with the number of times something flowered
            group_by(cgPlaId) %>%
            mutate(phenCt = n())


# lydia start here!
# pedigree data 
ped96 <- ped[ped$cgplaid < 647, ]
names(ped96)[names(ped96)=="siteOfOriginPedigree"] <- "site"
names(ped96)[names(ped96)=="cgplaid"] <- "cgPlaId"
# fixing levels in pedigree data
levels(ped96$site)[levels(ped96$site) == "AA"] <- "aa"
levels(ped96$site)[levels(ped96$site) == "Eriley"] <- "eri"
levels(ped96$site)[levels(ped96$site) == "Lf"] <- "lf"
levels(ped96$site)[levels(ped96$site) == "NWLF"] <- "nwlf"
levels(ped96$site)[levels(ped96$site) == "SPP"] <- "spp"
levels(ped96$site)[levels(ped96$site) == "Stevens"] <- "sap"
levels(ped96$site)[levels(ped96$site) == "Nessman"] <- "ness"

# merging data
x1996 <- merge(x1996, ped96, by.x = "cgPlaId", by.y = "cgPlaId")

### get rid of empty levels 
x1996 <- x1996[x1996$site != "", ]
x1996 <- x1996[x1996$site != "Unknown", ]
x1996 <- x1996[x1996$site != "Tower", ] # Only 1 plant flowered
x1996$site <- x1996$site[ , drop=TRUE]
levels(x1996$site)

# adding in row and position data 
x1996 <- merge(x1996, rowpos, by.x = "cgPlaId", by.y = "cgPlaId")

x2013 <- x1996[x1996$year=="2013", ]
x2012 <- x1996[x1996$year=="2012", ]
x2011 <- x1996[x1996$year=="2011", ]
x2010 <- x1996[x1996$year=="2010", ]
x2009 <- x1996[x1996$year=="2009", ]
x2008 <- x1996[x1996$year=="2008", ]
x2007 <- x1996[x1996$year=="2007", ]
x2006 <- x1996[x1996$year=="2006", ]
x2005 <- x1996[x1996$year=="2005", ]









p.05.12 <- p.05.12x %>%
  dplyr::rename("cgHdId" = "cgHeadId") %>%
  # making column for number of heads
  group_by(year, cgPlaId) %>%
  mutate(headCt = n()) %>%
  # converting dates to POSIXt
  mutate_at(vars(startDateEarly:endDateLate), as.character) %>%
  mutate_at(vars(startDateEarly:endDateLate), lubridate::ymd)



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

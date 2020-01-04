# Script to tidy Echinacea phenology data and save as an .rda file or use in a markdown doc
#       created by:    Lydia English
#       last modified: 12/15/2019 by LE

# # Load packages # # 

library(tidyverse) # using tidy syntax
library(lubridate) # to help with dates

# # Load phenology datasets # #

p.05.17x <- read_csv("data-raw/exPt1Phenology.csv")                 # data from 2005 - 2017
p.10x <- read_delim("data-raw/2010.CG1.Phenology.txt", delim = ",") # 2010 data
ped <- read_csv("data-raw/96979899qGenPedigreeLE.csv")              # pedigree data for site origins
rowpos <- read_csv("data-raw/cg1CoreData.csv") %>%                  # row/position data
  select(cgPlaId:yrPlanted) %>%                                     # using all of CG1 data since I only had row/pos for 1996
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  mutate(yrPlanted = as.factor(yrPlanted))
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

# 1996 and 1997 datasets
x19967 <- p.all %>% filter(cgPlaId < 1237) %>%
            # add columns with julian date
          dplyr::mutate(startNum = yday(startDtEarly))%>%
          dplyr::mutate(endNum   = yday(endDtLate))%>%
            # add column with the number of times something flowered
          group_by(cgPlaId) %>%
          dplyr::mutate(phenCt = n()) %>%
          ungroup() 
# how many plants just flowered once
x19967 %>%
  filter(phenCt == 1) %>% distinct(cgPlaId) %>% tally()

x19967 <- x19967 %>%
          filter(phenCt > 1)

# pedigree data 
ped967 <- ped %>%
  filter(cgplaid < 1237) %>%
  # rename columns
  dplyr::rename("site" = "siteOfOriginPedigree")%>%
  dplyr::rename("cgPlaId" = "cgplaid")%>%
  # recode levels
  mutate(site = recode(site, 
                       "AA"      = "aa",
                       "Aanenson"= "aa",
                       "Eriley"  = "eri",
                       "Riley"   = "ri",
                       "Lf"      = "lf",
                       "NWLF"    = "nwlf",
                       "NWlf"    = "nwlf",
                       "SPP"     = "spp",
                       "Stevens" = "sap",
                       "Staff"   = "spp",
                       "Nessman" = "ness"
                       ))

# merging data
phen_19967 <- left_join(x19967, ped967, by = "cgPlaId")%>%
  # get rid of empty levels - I'm going to ignore this for now, since I actually don't care about site origin. 
   filter(cgPlaId != 0
    #    site != "",
    #    site != "Unknown",
    #    site != "Tower",
         ) %>%
  # add in row position data
  left_join(., rowpos, by = "cgPlaId")%>%
  ungroup(year)%>%
  mutate(dur = endNum - startNum,
         year = forcats::as_factor(year)) %>%
  filter(dur < 100) %>%
  group_by(year)%>%
  mutate(medianDate = median(startDtEarly)) %>%
  mutate(dam = as.numeric(startDtEarly - medianDate)) %>%
  dplyr::mutate(rank = dense_rank(startDtEarly))%>%
  dplyr::mutate(daf = startDtEarly - min(startDtEarly)) %>%
  ungroup(year) %>%
  dplyr::mutate(site = forcats::as_factor(site))

# # saving data for easy access in the future # #  
# save(phen_19967, file = "phen_dataset.rda")

# Also going to save separate dataframes for 1996 and 1997
phen_1996 <- phen_19967 %>% filter(cgPlaId < 647)
# save(phen_1996, file = "phen96_dataset.rda")
phen_1997 <- phen_19967 %>% filter(646 < cgPlaId & cgPlaId < 1237)
# save(phen_1997, file = "phen97_dataset.rda")

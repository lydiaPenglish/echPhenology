# Script to make Echinacea consistency figures # 
#     created by: Lydia English
#     last updates: 12/12/2019

# load packages and data
library(tidyverse)
library(patchwork) # for plottig multiple graphs
library(lmerTest)
library(viridis)
load("data/phen_dataset.rda") #both 1996 and 1997 dataset
theme_set(theme_bw())

# # Summary figs # # 

# Distribution of head count in each year
phen_19967 %>%
  ggplot()+
  geom_histogram(aes(headCt))+
  facet_wrap(~year)

# comparing with duration
p1 <- phen_19967 %>%
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

p2 <- phen_19967 %>%
  ggplot(aes(year, headCt))+
  geom_count(alpha = 0.2)+
  stat_summary(fun.y = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)+
  labs(x = NULL, 
       y = "Head Count",
       size = "Count")+
  theme(axis.text.x = element_text(angle = 45, size = rel(1.1), hjust = 1),
        axis.title.y = element_text(size = rel(1.15)),
        legend.background = element_rect(color = "black"))

p1 + p2
phen_19967 <- phen_19967 %>%
  mutate(yearN = as.numeric(year))
# should year be a factor or continuous?
m1 <- lmerTest::lmer(dur ~ headCt + year + expNm + (1|cgPlaId), data = phen_19967)
summary(m1)
anova(m1)

phen11 <- phen_19967 %>%
  filter(year == 2011)

l1 <- lm(dur ~ headCt + I(headCt^2), data = phen11)
summary(l1)

ggplot(phen_19967, aes(headCt, dur))+
  geom_point(alpha = 0.2)+
  facet_wrap(~year)

# # Trying to visualize with plants are flowering each year and when they start # #

# all row/pos

rowpos <- read_csv("data-raw/cg1CoreData.csv") %>%
  select(cgPlaId:yrPlanted) %>%
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  mutate(yrPlanted = as.factor(yrPlanted))

phen_19967 %>%
  filter(year != 2017) %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.5)+
  geom_point(aes(row, pos, color = dam))+
  scale_color_viridis(option="magma")+
  facet_wrap(~year)

phen_19967 %>%
  ggplot()+
  #geom_point(data = rowpos, aes(row, pos), size = 0.5)+
  geom_point(aes(row, pos, size = phenCt), color = "black", fill = "white")



# # trying out a plot where I look at plants that flowered more than 7 times and their spread of flowering times
# this one has standarized dates (mean +/- error bar)
phen_19967 %>%
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

# this one has boxplot - more full spread
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

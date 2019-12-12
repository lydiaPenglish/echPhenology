# Script to make Echinacea consistency figures # 
#     created by: Lydia English
#     last updates: 12/12/2019

# load packages and data
library(tidyverse)
library(patchwork) # for plottig multiple graphs
library(lmerTest)
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

l1 <- lmerTest::lmer(dur ~ headCt + (1|year), data = phen_19967)
summary(l1)

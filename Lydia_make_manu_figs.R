# Separate script to make some of figures for the manuscript #
# If figures aren't located in this script then the header will tell you what document to find them in
# All figures are saved in echPhenology/figs 

library(tidyverse)
library(patchwork)
library(RColorBrewer)
data("phen_dataset")                                     # phenology data
census_dat <- read_csv("data-raw/cg1CoreData.csv") %>%   # getting census data 
  # filter for columns with just flowering info
  select(cgPlaId, yrPlanted, starts_with("fl")) %>%
  # filter for plants in 1996 and 1997 cohort
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  # go from wide to long
  pivot_longer(cols = starts_with("fl"), names_to = "fl_year") %>%
  filter(value != 0)  %>%
  # get rid of the "fl" in front of all year values, make numeric
  mutate(fl_year = as.numeric(unlist(str_extract_all(fl_year, "(?<=fl)\\d{4}"))))

# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# ---- Figure 1. Census data summary, 4 graphs in one -----
#      i. Age when plants first begin flowering
#      ii. Number of plants flowering in each year  
#      iii. Number of times plants flowered within the study period
#      iv. Interval between flowering periods

# i. Age
firstYrs <- census_dat %>%
  group_by(cgPlaId) %>%
  mutate(firstYr = min(fl_year)) %>%
  ungroup() %>%
  distinct(cgPlaId, yrPlanted, firstYr) %>%
  mutate(ageAtFl = firstYr - yrPlanted)

# average 

firstYrs %>%
  group_by(yrPlanted) %>%
  summarize(avgAge = mean(ageAtFl),
            minAge = min(ageAtFl),
            maxAge = max(ageAtFl),
            medAge = median(ageAtFl))

p1 <- firstYrs %>%
  ggplot(aes(ageAtFl))+
  geom_histogram(stat = "count", color = "black", size = 1, fill = my_cols[8]) +
  labs(x = "Age at First Flower", y = "Number of Plants")+
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)) +
  facet_wrap(~yrPlanted, ncol = 1)+
  theme(axis.text        = element_text(size = rel(1.25)),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text       = element_text(size = rel(1.4)),
        axis.title       = element_text(size = rel(1.5)))
p1

# ii. Num of plants histogram
burn_years <- c(2006, 2008, 2011, 2013, 2015, 2018)

p2 <- census_dat %>%
  mutate(fl_year = as.factor(fl_year),
         burn = if_else(fl_year %in% burn_years, "burned", "not_burned")) %>%
  ggplot(aes(fl_year))+
  geom_histogram(aes(fill = burn), stat = "count", color = "black", size = 1)+
  scale_fill_manual(values = c(my_cols[2], my_cols[12]))+
  scale_x_discrete(breaks = c("1999", "2002", "2005", "2008", "2011", "2014",
                              "2017"))+
  guides(fill = FALSE)+
  labs(y = NULL,
       x = NULL)+
  theme(axis.text = element_text(size = rel(1.25)))
p2 
# iii. Phen Count histogram

census_phenCt <- census_dat %>%
  group_by(cgPlaId)%>%
  summarize(phenCt = sum(value)) %>%
  filter(phenCt != 1)
# Woah cgPlaID 372 has flowered 16 times from 1996 to 2018

census_phenCt %>%
  summarize(meanPC = mean(phenCt),
            maxPC  = max(phenCt),
            medPC  = median(phenCt))

p3 <- census_phenCt %>%
  ggplot()+
  geom_bar(aes(phenCt), size = 1, fill = my_cols[4], color = "black")+
  labs(x = "Phenology count",
       y = NULL)+
  theme(axis.text    = element_text(size = rel(1.25)),
        axis.title.x = element_text(size = rel(1.3)))
p3
#ggsave("phenCt_hist.png", plot = plot_pc) 

# iv. interval histogram 

fl_intv <- census_dat %>%
  select(cgPlaId, fl_year) %>%
  group_by(cgPlaId) %>%
  mutate(intv = fl_year - lag(fl_year, default = fl_year[1]))
fl_intv_sum <- fl_intv %>%
  select(-fl_year) %>%
  filter(intv != 0) %>%
  summarize(avg_intv = mean(intv),
            phenCt   = n() + 1) # No. of times a plant flowerd, adding one bc I took out a year when filtering for non-zero values
fl_intv_sum %>%
  summarize(mean_tot = mean(avg_intv),
            med_tot  = median(avg_intv),
            max_tot  = max(avg_intv),
            min_tot  = min(avg_intv))  
p4 <- fl_intv_sum %>%
  ggplot(aes(avg_intv))+
  geom_histogram(binwidth = 1, fill = my_cols[4], color = "black", size = 1) +
  #geom_vline(xintercept = 2.27, size = 1.5, lty = 2)+
  labs(x = "Average flowering interval",
       y = NULL)+
  scale_x_continuous(breaks = c(1:12))+
  theme(axis.text  = element_text(size = rel(1.25)),
        axis.title = element_text(size = rel(1.3)))
p4

# **** putting all the plots together ****

p1 + p2 / (p3 + p4)  + plot_annotation(tag_levels = 'A')

# ---- Figure 2. Phenology data summary, 2 graphs in one  ----
      # i. Timing of flowering phenology by year
      # ii. Duration of flowering phenology by year

# i. Timing 

burn_years <- c(2006, 2008, 2011, 2013, 2015)

phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) %>%
  group_by(year) %>%
  mutate(medianFFD = median(startNum))

b1 <- phen_all %>%
  #hmmm I'm not sure this is the best way to to this....
  ggplot(aes(year, as.Date(startNum, origin = "2018-01-01")))+
  geom_count(aes(color = burn), alpha = 0.5)+
  geom_errorbar(aes(ymin = as.Date(medianFFD, origin = "2018-01-01"), 
                    ymax = as.Date(medianFFD, origin = "2018-01-01")),
                width = 0.5, size = 1.5)+
  labs(y = "First flowering \nday",
       x = NULL, 
       size = "Count")+
  guides(color = FALSE)+
  scale_size_area()+
  scale_color_manual(values = c(my_cols[2], my_cols[12]))+
  #scale_size_area(max_size = 8)+
  theme(axis.text.x = element_blank(),
        axis.title  = element_text(size = rel(1.3)),
        legend.background = element_rect(color = "black"))
b1

# ii. Duration

b2 <- phen_all %>%
  group_by(year) %>%
  mutate(medianDur = median(dur))%>%
  ggplot(aes(year, dur))+
  geom_count(aes(color = burn), alpha = 0.5)+
  geom_errorbar(aes(ymin = medianDur, 
                    ymax = medianDur),
                width = 0.5, size = 1.5)+
  labs(x = NULL, 
       y = "Flowering duration \n(days)",
       size = "Count")+
  guides(color = FALSE)+
  scale_size_area()+
  scale_color_manual(values = c(my_cols[2], my_cols[12]))+
  theme(axis.text.x = element_text(size = rel(1.25)),
        axis.title.y = element_text(size = rel(1.3)),
        legend.background = element_rect(color = "black"))
b2
# **** together

b1/b2 + plot_annotation(tag_levels = 'A')

# ---- Figure 3. Correlation analysis ----

# This figure can be found in the script "Lydia_correlation_analysis.R" at the bottom where is says # visualization

# ---- Figure S1. CG1 ----

# this figure is generated in "Lydia_summary_data.R" under the tab "Plot of Common Garden Experiment"

# ---- Figure S2. Peak flowering day ----

# Plot located in "Lydia_summary_data.R" in the section:
# "How much does "peak" flowering day vary by year" 

# ---- Figure S3. Relationship between head count and duration of flowering time ----

h1 <- phen_19967 %>%
  ggplot(aes(headCt, dur))+
  geom_point(alpha = 0.2, size = 2)+
  geom_smooth(color = my_cols[10], size = 1.5, se = FALSE)+
  #geom_smooth(method = "lm", se = FALSE, color = my_cols[1], size = 1.25)+
  #geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = my_cols[2], size = 1.25)+
  #geom_smooth(method = "lm", formula = y ~ x + I(x^3), se = FALSE, color = my_cols[3], size = 1.25)+
  facet_wrap(~year) +
  labs(y = "Flowering duration (days)", x = "Head Count")+
  theme(strip.background = element_rect(fill = "white"),
        strip.text       = element_text(size = rel(1.25)),
        axis.title       = element_text(size = rel(1.25)),
        axis.text        = element_text(size = rel(1.2)))
h1

# ---- Figure S4. Z-scores for Moran's I analysis of FFD ----

# Plot located in "Lydia_SAC_analysis.R" under 
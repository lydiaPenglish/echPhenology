# Separate script to make some of figures for the manuscript #
# If figures aren't located in this script then the header will tell you what document to find them in
# All figures are saved in echPhenology/figs 

library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
data("phen_dataset")                                     # phenology data
census_dat <- readr::read_csv("data-raw/cg1CoreData.csv") %>%   # getting census data 
  # filter for columns with just flowering info
  select(cgPlaId, yrPlanted, starts_with("fl")) %>%
  # filter for plants in 1996 and 1997 cohort
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  # go from wide to long
  tidyr::pivot_longer(cols = starts_with("fl"), names_to = "fl_year") %>%
  filter(value != 0)  %>%
  # get rid of the "fl" in front of all year values, make numeric
  mutate(fl_year = as.numeric(unlist(stringr::str_extract_all(fl_year, "(?<=fl)\\d{4}"))))

# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# ---- Figure 1. Census data summary, 4 graphs in one -----
#      i. Age when plants first begin flowering
#      (ii. Number of plants flowering in each year) MOVED TO FIGURE 2  
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
            sdAge  = sd(ageAtFl),
            minAge = min(ageAtFl),
            maxAge = max(ageAtFl),
            medAge = median(ageAtFl))

p1 <- firstYrs %>%
  ggplot(aes(ageAtFl))+
  geom_histogram(stat = "count", color = "black", size = 1, fill = "#C669D5") +
  labs(x = "Age at First Flower", y = NULL)+
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)) +
  facet_grid(rows = vars(yrPlanted))+
  theme(axis.text        = element_text(size = rel(1.25)),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text       = element_text(size = rel(1.4)),
        axis.title       = element_text(size = rel(1.5)))
p1

# iii. Phen Count histogram

census_phenCt <- census_dat %>%
  filter(fl_year != 2018) %>%
  group_by(yrPlanted, cgPlaId)%>%
  summarize(phenCt = sum(value)) 
#%>%
 # filter(phenCt != 1)
# Woah cgPlaID 372 has flowered 16 times from 1996 to 2018

census_phenCt %>%
  summarize(meanPC = mean(phenCt),
            sdPC   = sd(phenCt),
            maxPC  = max(phenCt),
            medPC  = median(phenCt))

p3 <- census_phenCt %>%
  ggplot()+
  geom_bar(aes(phenCt), size = 1, fill = "#B3B3B3", color = "black")+
  facet_grid(rows = vars(yrPlanted))+
  scale_x_continuous(breaks = c(1, 5, 10, 15))+
  labs(x = "Flowering count",
       y = "Number of flowering plants")+
  theme(axis.text    = element_text(size = 13),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text       = element_text(size = 16),
        axis.title = element_text(size = 14))
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
            sd_tot   = sd(avg_intv),
            med_tot  = median(avg_intv),
            max_tot  = max(avg_intv),
            min_tot  = min(avg_intv))  
p4 <- fl_intv_sum %>%
  ggplot(aes(avg_intv))+
  geom_histogram(binwidth = 1, fill = "#B3B3B3", color = "black", size = 1) +
  #geom_vline(xintercept = 2.27, size = 1.5, lty = 2)+
  labs(x = "Flowering Interval",
       y = NULL)+
  scale_x_continuous(breaks = c(1:12))+
  theme(axis.text  = element_text(size = rel(1.25)),
        axis.title = element_text(size = rel(1.5)))
p4

# **** putting all the plots together ****

p1 + p3 / p4  + plot_annotation(tag_levels = 'A')

# ---- Figure 2. Phenology data summary, 2 graphs in one  ----
      # i. Number of flowering plants (moved from FIG 1)
      # ii. Timing of flowering phenology by year
      # iii. Duration of flowering phenology by year

burn_years <- c(2006, 2008, 2011, 2013, 2015)
phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) %>%
  group_by(year) %>%
  mutate(medianFFD = median(startNum))


# ii. Num of plants histogram 
p1 <- phen_all %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) %>%
  ggplot(aes(year))+
  geom_histogram(aes(fill = burn), stat = "count", color = "black", size = 1)+
  scale_fill_manual(values = c(my_cols[2], "#B3B3B3"))+
  scale_x_discrete(breaks = c("2005", "2007", "2009", "2011", "2013", "2015", 
                              "2017"))+
  guides(fill = FALSE)+
  labs(y = "Number of Flowering\n Plants",
       x = NULL)+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))
p1 

# i. Timing 

p2 <- phen_all %>%
  # hmmm I'm not sure this is the best way to to this....
  ggplot(aes(year, as.Date(startNum, origin = "2018-01-01"))) +
  geom_count(aes(fill = burn), color = "black", pch = 21) +
  geom_errorbar(aes(
    ymin = as.Date(medianFFD, origin = "2018-01-01"),
    ymax = as.Date(medianFFD, origin = "2018-01-01")
  ),
  width = 0.5, size = 1.5
  ) +
  labs(
    y = "First flowering Day",
    x = NULL,
    size = "Count"
  ) +
  guides(fill = FALSE) +
  scale_size_area(breaks = c(20, 40, 60, 80), max_size = 8) +
  scale_fill_manual(values = c(my_cols[2], "#B3B3B3")) +
  scale_y_date(date_breaks = "2 weeks", date_labels = "%b-%d") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.background = element_rect(color = "black")
  )
p2

p2/p1

# ii. Duration

p3 <- phen_all %>%
  group_by(year) %>%
  mutate(medianDur = median(dur)) %>%
  ggplot(aes(year, dur)) +
  geom_count(aes(fill = burn), color = "black", pch = 21) +
  geom_errorbar(aes(
    ymin = medianDur,
    ymax = medianDur
  ),
  width = 0.5, size = 1.5
  ) +
  labs(
    x = NULL,
    y = "Flowering duration",
    size = "Count"
  ) +
  guides(fill = FALSE) +
  scale_size_area(breaks = c(20, 40, 60, 80), max_size = 8) +
  scale_fill_manual(values = c(my_cols[2], "#B3B3B3")) +
  scale_x_discrete(breaks = c("2005", "2007", "2009", "2011", "2013", "2015", 
                              "2017"))+
  theme(
    axis.text.x = element_text(size = rel(1.25)),
    axis.title.y = element_text(size = rel(1.5)),
    legend.background = element_rect(color = "black")
  )
p3
# **** together

p1 + p2/p3 + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.5, 2))

# number of times plants have flowered

phen_all %>%
  ungroup() %>%
  select(cgPlaId, phenCt, yrPlanted) %>%
  distinct() %>%
  ggplot(aes(phenCt))+
  geom_histogram(stat = "count", color = "black", size = 1, fill = "#C669D5")+
  labs(x = "Phenology Count", y = NULL)+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))

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
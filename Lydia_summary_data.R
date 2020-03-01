# # Script for getting summary data and making summary figures about Echinacea phenology # #

# loading libraries and dataframes
library(tidyverse)
library(lme4)
library(lmerTest)
library(lubridate)
library(patchwork)
data("phen_dataset")
data("phen96_dataset")
data("phen97_dataset")
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

# graphing parameters
theme_set(theme_bw())
library(RColorBrewer)
brewer.pal(n = 8, name = "Dark2")
my_cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
             "#A6761D","#666666")
display.brewer.pal(n = 8, name = "Dark2")
display.brewer.pal(n = 7, name = 'PuBuGn')
brewer.pal(n = 7, name = "PuBuGn")

# ---- Plot of the common garden experiment ----
rowpos <- read_csv("data-raw/cg1CoreData.csv") %>%
  select(cgPlaId:yrPlanted) %>%
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  mutate(yrPlanted = as.factor(yrPlanted))

cg1 <- rowpos %>%
  ggplot(aes(row, pos))+
  geom_point(shape = 4)+
  geom_hline(yintercept = 959.5, lty = 2, size = 1)+
  coord_fixed()+
  labs(x = NULL, y = NULL)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
cg1
ggsave("common_garden.png", plot = cg1, path = "./figs")

# ---- Some basic information - short code snippets ----

# 1. Number of sites that plants came from 
ped <- read_csv("data-raw/96979899qGenPedigreeLE.csv") 
ped_sites <- ped %>%
  group_by(expNm, siteOfOriginPedigree) %>%
  summarize(n = n())

# 2. Number of plants that didn't flower 

phen_19967 %>%
  anti_join(rowpos, ., by = "cgPlaId") # n = 770

phen_19967 %>%
  distinct(cgPlaId) %>% tally()

# 3. How many times did plants flower on avg
census_phenCt <- census_dat %>%             
  group_by(cgPlaId)%>%
  summarize(phenCt = sum(value)) %>%
  filter(phenCt != 1)

census_phenCt %>%                           # mean = 6.63 within census period
  summarize(mPC   = mean(phenCt),
            maxPC = max(phenCt))

phen_19967 %>%                              # mean = 5.63 within phenology period
  summarize(mPC   = mean(phenCt),
            maxPC = max(phenCt))

# 4. Number of plants that flowered more than 10 times? Answer = 2
phen_19967 %>%
  filter(phenCt > 10) %>%
  distinct(cgPlaId)

# 5. Range of number of flowering plants per year

phen_19967 %>%
  group_by(year) %>%
  summarize(n = n())

# ---- Avg age of a plant when it first flowered  and flowering interval ----

# get earliest FFD for each cgPlaId
minYrsToFl <- phen_19967 %>%
  group_by(cgPlaId, expNm) %>%
  summarize(earliest = min(startDtEarly)) %>%         # 3 plants only flowered 2014, 2015 - 362, 552, 1076
  arrange(desc(earliest)) %>%
  mutate(minYrFl = year(earliest),
         YrsToFl = minYrFl - 2005 + 1)

minYrsToFl %>%                                        # majority of plants flowered within first 2 years of experiment, but others waited
  ggplot(aes(YrsToFl))+
  geom_bar()    

# cross referencing CG1 core dataset

wide_info <- read_csv("data-raw/cg1CoreData.csv") %>%
  select(cgPlaId, yrPlanted, starts_with("fl")) %>%
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  select(-yrPlanted) %>%
  pivot_longer(cols = starts_with("fl"), names_to = "year") %>%
  pivot_wider(names_from = cgPlaId, values_from = value) %>%
  mutate(year = str_extract_all(year, "(?<=fl)\\d{4}"))

long_info <- read_csv("data-raw/cg1CoreData.csv") %>%
  select(cgPlaId, yrPlanted, starts_with("fl")) %>%
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  pivot_longer(cols = starts_with("fl"), names_to = "year") %>%
  filter(value != 0)  %>%
  mutate(year = as.numeric(unlist(str_extract_all(year, "(?<=fl)\\d{4}"))))
  
firstYrs <- long_info %>%
  group_by(cgPlaId) %>%
  mutate(firstYr = min(year)) %>%
  ungroup() %>%
  distinct(cgPlaId, yrPlanted, firstYr) %>%
  mutate(ageAtFl = firstYr - yrPlanted)

FYs <- firstYrs %>%
  ggplot(aes(ageAtFl))+
  geom_histogram(stat = "count", color = "black", size = 1, fill = "#D0D1E6") +
  labs(x = "Plant Age", y = NULL)+
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21)) +
  facet_grid(rows = vars(yrPlanted))+
  theme(axis.text        = element_text(size = rel(1.25)),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text       = element_text(size = rel(1.5)),
        axis.title       = element_text(size = rel(1.3)))

FYs
#ggsave("AgePlantsBeginFlowering.png", FYs, path = "./figs")

# looking at flowering interval 
fl_intv <- long_info %>%
  select(cgPlaId, year) %>%
  group_by(cgPlaId) %>%
  mutate(intv = year - lag(year, default = year[1]))
fl_intv_sum <- fl_intv %>%
  select(-year) %>%
  filter(intv != 0) %>%
  summarize(avg_intv = mean(intv),
            phenCt   = n() + 1) # No. of times a plant flowerd, adding one bc I took out a year when filtering for non-zero values
fl_intv_sum %>%
  summarize(mean_tot = mean(avg_intv))       # overall average interval is 2.27yrs

# graphics
# 1. histogram
fl_intv_plot <- fl_intv_sum %>%
  ggplot(aes(avg_intv))+
  geom_histogram(binwidth = 1, fill = "white", color = "black", size = 1) +
  geom_vline(xintercept = 2.27, size = 1.5, lty = 2)+
  labs(x = "Average number of years \nbetween flowering years",
       y = NULL)+
  scale_x_continuous(breaks = c(1:12))+
  theme(axis.text  = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.25)))
fl_intv_plot
ggsave("dist_flowering_interval.png", plot = fl_intv_plot, path = "./figs")

# somehow want to visualize interval with phenCt

fl_intv_phenCt <- fl_intv_sum %>%
  ggplot(aes(phenCt, avg_intv))+
  geom_boxplot(aes(group = phenCt), size = 1)+
  geom_jitter(width = 0.25, alpha = 0.5, color = "#A6BDDB")+
  scale_size_continuous(breaks = c(5,10,15,20))+
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14, 16))+
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13))+
  labs(x = "Number of times plants flowered",
       y = "Average number of years between \nflowering periods")+
  theme(axis.text  = element_text(size = rel(1.15)),
        axis.title = element_text(size = rel(1.25)))
fl_intv_phenCt
ggsave("fl_interval_phenCt.png", plot = fl_intv_phenCt, path = "./figs")

# calculating for just phenology dataset
# fl_intv <- phen_19967 %>%
#   select(cgPlaId, year) %>%
#   arrange(cgPlaId) %>%
#   group_by(cgPlaId) %>%
#   mutate(year = as.numeric(year), 
#          intv = year - lag(year, default = year[1]))


# ---- Distribution of and relationship between head count and duration ----
burn_years <- c(2006, 2008, 2011, 2013, 2015)

phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned"))

# # -- Distribution of head count and duration in each year -- # # 
phen_19967 %>%
  ggplot()+
  geom_histogram(aes(headCt))+
  facet_wrap(~year)

phen_19967 %>%
  group_by(year) %>%
  summarize(meanD = mean(dur),
            sdD   = sd(dur),
            cvD   = sdD/meanD*100,
            nD    = n())

# comparing with duration
p1 <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) %>%
  ggplot(aes(year, dur))+
  geom_count(aes(color = burn), alpha = 0.5)+
  stat_summary(fun.y = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)+
  labs(x = NULL, 
       y = "Flowering duration (days)",
       size = "Count")+
  guides(color = FALSE)+
  scale_color_manual(values = c(my_cols[2], my_cols[8]))+
  theme(axis.text.x = element_text(angle = 45, size = rel(1.2), hjust = 1),
        axis.title.y = element_text(size = rel(1.3)),
        legend.background = element_rect(color = "black"),
        legend.position = "bottom")

p2 <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) %>%
  ggplot(aes(year, headCt))+
  geom_count(aes(color = burn), alpha = 0.5)+
  stat_summary(fun.y = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)+
  labs(x = NULL, 
       y = "Head Count",
       size = "Count")+
  guides(color = FALSE)+
  scale_color_manual(values = c(my_cols[2], my_cols[8]))+
  theme(axis.text.x = element_text(angle = 45, size = rel(1.2), hjust = 1),
        axis.title.y = element_text(size = rel(1.3)),
        legend.background = element_rect(color = "black"),
        legend.position = "bottom")
hdCt_dur <- p1 + p2
hdCt_dur
ggsave("HeadCt_and_duration.png", plot = hdCt_dur)

# distribution of head count in each year (in 1996 vs 1997 garden)
# this isn't the best way to visualize this, but hard with a poison distribution (hist wasn't great, boxplots were garbage)
hdCt_96_97 <- phen_19967 %>%
  mutate(expNm = as_factor(as.character(expNm)))%>%
  ggplot(aes(expNm, headCt))+
  geom_count(alpha = 0.5)+
  facet_wrap(~year)+
  labs(y = "Head Count", x = NULL, n = "Number of plants")+
  coord_flip()+
  theme(strip.background = element_rect(fill = "white"))
hdCt_96_97
ggsave("hdCt_96_vs_97.png", plot = hdCt_96_97)

# relationship between head Ct and duration - this isn't good. Fix this later....

# distribution of duration
phen_19967 %>%
  ggplot(aes(dur))+
  geom_histogram()+
  facet_wrap(~year)

m1 <- lmerTest::lmer(dur ~ headCt + year + expNm + (1|cgPlaId), data = phen_19967)
# ^ mixed model where year is a fixed effect...says that dur and headct are related
ggResidpanel::resid_panel(m1)
summary(m1)
anova(m1)

m2 <- lm(dur ~ headCt + burn, data = phen_all)
summary(m2)
car::Anova(m2, type = "II")

# OR.... here I am doing an individual model for each year
nested_df <- phen_19967 %>%
  group_by(year) %>%
  nest() %>%
  mutate(model = map(data, ~lm(dur ~ headCt + expNm, data = .)),
         tidy = map(model, broom::tidy))
models <- as_tibble(nested_df %>%
  unnest(tidy) %>%
  select(year, term:p.value))
models %>%
  filter(p.value < 0.05) # always significant (althought 2017 is weird and in opposite direction)

l1 <- lm(dur ~ headCt + I(headCt^2), data = phen11)
summary(l1)
# plot of headCt vs duration 
dur_hd_yr <- phen_19967 %>%
  filter(year != 2017) %>%
  ggplot(aes(headCt, dur))+
  geom_point(alpha = 0.2, size = 2)+
  geom_smooth(method = "lm", se = FALSE, color = my_cols[1], size = 1.25)+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = FALSE, color = my_cols[2], size = 1.25)+
  geom_smooth(method = "lm", formula = y ~ x + I(x^3), se = FALSE, color = my_cols[3], size = 1.25)+
  facet_wrap(~year) +
  labs(y = "Duration", x = "Head Count")+
  theme(strip.background = element_rect(fill = "white"),
        axis.title       = element_text(size = rel(1.25)),
        axis.text        = element_text(size = rel(1.2)))
dur_hd_yr
ggsave("dur_vs_hd_by_yr.png", plot = dur_hd_yr)

# # -- plot of average head count as a function of row/pos -- # # 
library(viridis)
avg_hd_ct <- phen_19967 %>%
  group_by(cgPlaId) %>%
  mutate(meanHdCt = mean(headCt)) %>%
  ungroup() %>%
  distinct(cgPlaId, row, pos, meanHdCt) %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.25)+
  geom_point(aes(row, pos, color = meanHdCt), size = 1.5)+
  geom_hline(yintercept = 959.5, lty = 2)+
  labs(x = NULL, y = NULL, color = "Avg Head Count")+
  coord_fixed()+
  scale_color_viridis(direction = -1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(color = "black"))
avg_hd_ct
ggsave("avgHdCt_row_pos.png", plot = avg_hd_ct)

# ---- Average duration per flowering head; how much does that vary per year ----

# need raw data again
p.05.17x <- read_csv("data-raw/exPt1Phenology.csv") %>%
  dplyr::rename("year" = "phenYear")%>%
  group_by(year, cgPlaId) %>%
  # getting number of heads, deleting cgHdId
  dplyr::mutate(headCt = n()) %>%
  ungroup() %>%
  select(-c(row, pos, ttColor, phenNote)) %>%
  filter(startDtEarly > 1940-01-01 & endDtLate > 1940-01-01) %>%
  mutate(dur = as.numeric(endDtLate - startDtEarly),
         year_f = as.factor(as.character(year))) %>%
  # this gets out one bad ones (wrong dates or end dates before start dates)
  filter(dur < 40 & dur > 0)
# average across all years
p.05.17x %>%
  group_by(year) %>%
  summarize(meanD = mean(dur),
            sdD   = sd(dur),
            cvD   = sdD/meanD*100,
            nD    = n())

p2 <- p.05.17x %>%
  ggplot(aes(year_f, dur))+
  geom_count()+
  labs(x = NULL,
       y = "Flowering duration (days)",
       size = "Number of flowering heads")+
  stat_summary(fun.y = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.3,
               size = 2,
               color = "#7570B3")+
  scale_size_continuous(breaks = c(50,100,200,400))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color = "black"),
        axis.title.y = element_text(size = rel(1.1)))

ggsave("duration_heads.png", plot = p2)

# ---- Do burn years have a longer flowering duration than non-burn years? ----

burn_years <- c(2006, 2008, 2011, 2013, 2015)

phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned"))

# mixed model, just looking at burn/not burn
m1 <- lmerTest::lmer(dur ~ burn + (1|cgPlaId), data = phen_all)
summary(m1)
anova(m1)
ggResidpanel::resid_panel(m1) # this doesn't fit super well 

m2 <- lmer(dur ~ burn*year + (1|cgPlaId), data = phen_all)
summary(m2)
anova(m2)

hist(phen_all$dur); mean(phen_all$dur) ; sd(phen_all$dur)

library(MASS)

g1 <- glm.nb(dur ~ burn, data = phen_all)
summary(g1)

g2 <- glmer.nb(dur ~ burn + (1|cgPlaId), data = phen_all)
summary(g2)
ggResidpanel::resid_panel(g2)

# wilcoxon rank sum test

w1 <- wilcox.test(dur ~ burn, phen_all)
w1

# burn years have more flowering plants? 
totals <- phen_19967 %>%
  group_by(year) %>%
  summarize(total_phen = n()) %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned"))

m4 <- lm(total_phen ~ burn, data = totals)  
summary(m4)
anova(m4)
ggResidpanel::resid_panel(m4)
confint(m4)

# ---- Is the FFD different between burn and non-burn years? Yes ---- 

burn_years <- c(2006, 2008, 2011, 2013, 2015)

phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned"))

l1 <- lm(startNum ~ burn, data = phen_all)
summary(l1)
l2 <- lmer(startNum ~ burn + (1|cgPlaId), data = phen_all)
summary(l2)
anova(l2)
ggResidpanel::resid_panel(l2)

# summary stat 
phen_summary <- phen_all %>%
  group_by(year) %>%
  summarize(avg_ffd = mean(startNum),
            n_ffd   = n(),
            sd_ffd  = sd(startNum)) %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned")) 

l2 <- lm(avg_ffd ~ burn, data = phen_summary)
summary(l2)

# ---- Do plants that flower earlier have more heads? FINISH later ----
# first try with one year 
phen_2006 <- phen_19967 %>% filter(year == "2006")
hist(phen_2006$headCt); mean(phen_2006$headCt) ; sd(phen_2006$headCt) # pretty close, poisson should be ok

g1 <- glm.nb(headCt ~ startNum, data = phen_2006)
summary(g1)
ggResidpanel::resid_panel(g1)
par(mfrow = c(2,2))
plot(g1)
par(mfrow = c(1,1)) # i think negative binomial looks better?

g2 <- glm(headCt ~ startNum, data = phen_2006)
summary(g2)
ggResidpanel::resid_panel(g2)

phen_19967 %>%
  ggplot(aes(startNum, headCt))+
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~year)

# try as a mixed model - there seems to be evidence for this. 

g2 <- glmer.nb(headCt ~ startNum + (1|cgPlaId), 
                  data = phen_19967)
summary(g2)
plot(g2)

g3 <- glmer(headCt ~ startNum + (1|cgPlaId), family = poisson(link = "log"), data = phen_19967)
summary(g3)
ggResidpanel::resid_panel(g3)

manyHds <- phen_19967 %>% filter(headCt > 6) %>%
  distinct(cgPlaId) %>% unlist()

lookies <- phen_19967 %>%
  filter(cgPlaId %in% manyHds) 

# are floweirng time and duration correlated?

phen_19967 %>%
  ggplot(aes(startNum, dur))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~year)

f1 <- lmer(dur ~ startNum*year + (1|cgPlaId), data = phen_19967)
summary(f1)
anova(f1)

# ---- How much does median and "peak" flowering day vary by year ----

# median FFD

phen_19967 %>%
  group_by(year) %>%
  summarize(medStartNum = median (startNum)) %>%
  mutate(medFFD = as_date(medStartNum))

# this is it's own section because it was a pain in the butt to do....

# what is peak flowering day? 
# note to self - tibbles don't work with EP packages...(grrrrrrr)
library(echPhenology)

# getting flowering schedules for each year using a loop
years <- phen_19967 %>% select(year) %>% distinct() %>% unlist() %>% unname() %>% as.character(.)
output <- list()
for(i in years){
  x <- phen_19967 %>% filter(year == i) %>% as.data.frame(.)
  FS <- makeFS(x, startCol = "startDtEarly", endCol = "endDtLate")
  output[[i]] <- summaryFS(FS)
}

#extracting peak dates for each year from loop output
peak_dates <- output %>%
  map_df("peakDate") %>%
  pivot_longer(everything(), names_to = "year", values_to = "peak_date") 

# summary data for peak dates
peak_dates %>%
  mutate(day_num = yday(peak_date)) %>%
  summarize(mean_peak   = mean(day_num),
            median_peak = median(day_num),
            sd_peak     = sd(day_num),
            se_peak     = sd_peak/sqrt(13),
            min_peak    = min(day_num),
            max_peak    = max(day_num)) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value")

peak_dates %>%
  mutate(day_num = yday(peak_date)) %>%
  summary(day_num)

#extracting number of flowering plants on each date for using in graphics
count_dates <- output %>%
  map_df("fl.density") %>%
  mutate(year = year(day))

# making into dataframe
FS_df <- phen_19967 %>% 
  left_join(., peak_dates, by = "year") %>%
  mutate(stDayMonth  = format(as.Date(startDtEarly), format = "%m-%d"),
         endDayMonth = format(as.Date(endDtLate),    format = "%m-%d")) %>%
  group_by(year) %>%
  mutate(ord = row_number(startDtEarly)) %>%
  ungroup() %>%
  dplyr::select(year, cgPlaId, startDtEarly, endDtLate, yrPlanted, peak_date, stDayMonth:ord)

# Tried to graph everything and just facet by year but the doesn't work with dates (WHYYYYYYYYY)
# Then I tried to just have month labels and add in peak flowering dates, but you can't add a line that's
# not a number :( 

# this is ok I guess, but it needs better x-axis tick marks and it doesn't have peak day highlighted
FS_df %>%
  ggplot()+
  geom_segment(aes(x = stDayMonth, xend = endDayMonth,
                   y = ord, yend = ord)) +
  labs(y = NULL, x = NULL)+
  facet_wrap(~year)

# New solution: make plot for every year and then patchwork the whole thing together (not ideal but satisfactory)
# example of one
p2005 <- FS_df %>%
  filter(year == 2005) %>%
  ggplot()+
  geom_segment(aes(x = startDtEarly, xend = endDtLate,
                   y = ord, yend = ord))+
  labs(x = NULL, y = NULL)+
  geom_vline(xintercept = as.Date("2005-07-11"))+
  geom_point(data = count_dates %>% filter(year == 2005),
             aes(day, count))+
  scale_x_date(limits = c(as.Date("2005-06-19"), as.Date("2005-09-01")))+
  scale_y_continuous(limits = c(0,300))
p2005

# Loop to get all the plots
plot_list = list()
for(i in years){
  peak_dt <- peak_dates %>% filter(year == i) %>% pull(peak_date)
  plot_list[[i]] <-  FS_df %>%
    filter(year == i) %>%
    mutate(mytitle = i) %>%
    ggplot()+
    geom_segment(aes(x = startDtEarly, xend = endDtLate,
                     y = ord, yend = ord), alpha = 0.25)+
    labs(x = NULL, y = NULL)+
    geom_vline(xintercept = peak_dt)+
    geom_point(data = count_dates %>% filter(year == i),
               aes(day, count), color = "#7570B3")+
    scale_x_date(limits = c(as.Date(paste(i,"-06-19", sep = "")), as.Date(paste(i,"-09-01", sep = ""))))+
    scale_y_continuous(limits = c(0,300))+
    facet_grid(. ~ mytitle)+
    theme(strip.background = element_rect(color = "black",
                                          fill = "white"),
          axis.text.x      = element_text(angle = 45, hjust = 1))
}

library(patchwork)
# figure out a more efficient way to do this....
all_FS <- wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]],
           plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]])
ggsave("flowering_schedules.png", plot = all_FS)

# ---- Trying to visualize which plants are flowering each year and when they start ----

# all row/pos
rowpos <- read_csv("data-raw/cg1CoreData.csv") %>%
  select(cgPlaId:yrPlanted) %>%
  filter(yrPlanted == 1996 | yrPlanted == 1997) %>%
  mutate(yrPlanted = as.factor(yrPlanted))

p05 <- phen_19967 %>%
  filter(year == 2005) %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.25)+
  geom_point(aes(row, pos, color = startNum), size = 1.5)+
  scale_color_viridis(limits = c(170,210), option="magma")+
  coord_fixed()+
  guides(color = FALSE)+
  theme_bw()
p06 <- phen_19967 %>%
  filter(year == 2006) %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.25)+
  geom_point(aes(row, pos, color = startNum), size = 1.5)+
  scale_color_viridis(limits = c(170,210), option="magma")+
  coord_fixed()+
  guides(color = FALSE)+
  theme_bw()
p07 <- phen_19967 %>%
  filter(year == 2007) %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.25)+
  geom_point(aes(row, pos, color = startNum), size = 1.5)+
  scale_color_viridis(limits = c(170,210), option="magma")+
  coord_fixed()+
  theme_bw()
library(patchwork)
composite <- (p05+p06+p07 + plot_layout(guides = 'collect'))
ggsave("ex_startNum050607.png", plot = composite)

all_yrs <- phen_19967 %>%
  filter(year != "2017") %>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.1)+
  geom_point(aes(row, pos, color = startNum), size = 1.5)+
  scale_color_viridis()+
  #scale_color_distiller(palette = 'PuBuGn')+
  #coord_fixed()+
  theme_bw()+
  facet_wrap(~year)+
  theme(strip.background = element_rect(fill = "white"))
all_yrs
ggsave("allyrs.png", plot = all_yrs)

phen_19967 %>%
  filter(year == 2014) %>%
  ggplot()+
  geom_text(aes(row, pos, label = cgPlaId))+
  coord_fixed()+
  theme_bw()

phen_ct <- phen_19967 %>%
  mutate(phenCt_f = as.character(phenCt),
         phenCt_f = ifelse(phenCt > 7, "8+", phenCt_f))%>%
  ggplot()+
  geom_point(data = rowpos, aes(row, pos), size = 0.5)+
  geom_point(aes(row, pos, color = phenCt_f), size = 2.5)+
  geom_hline(yintercept = 959.5, lty = 2)+
  labs(x = NULL, y = NULL, color = "Phenology Count")+
  coord_fixed()+
  scale_color_viridis(discrete = TRUE, direction = -1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(color = "black"))
ggsave("cg_phenCt.png", plot = phen_ct)

# ---- How different is flowering time between 1996 adn 1997 garden ---- 
phen_19967 %>%
  mutate(expNm = as_factor(as.character(expNm)))%>%
  ggplot(aes(expNm, startNum))+
  geom_count(alpha = 0.5)+
  stat_summary(fun.y = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin = ..y..),
               width = 0.4,
               size = 2)+
  facet_wrap(~year)+
  labs(y = "Start Num", x = NULL, n = "Number of plants")+
  coord_flip()+
  theme(strip.background = element_rect(fill = "white"))



# ---- Misc other figures, not very useful anymore ----

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
phen_7 <- phen_19967 %>%
  filter(phenCt > 7) %>%
  group_by(cgPlaId) %>%
  mutate(medDAM = median(dam)) %>%
  ungroup() %>%
  distinct(cgPlaId, .keep_all = TRUE) %>%
  arrange(desc(medDAM)) %>%
  mutate(newID = row_number()) %>%
  select(cgPlaId, medDAM, newID) %>%
  left_join(., phen_19967, by = "cgPlaId")
phen_7 %>% 
  ggplot(aes(x = newID, y = dam)) +
  geom_boxplot(aes(group = newID))+
  #geom_point(alpha = 0.2)+
  geom_hline(yintercept = 0, lty = 2) +
  coord_flip() 


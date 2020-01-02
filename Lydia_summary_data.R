# # Script for some summary data about Echinacea flowering phenology # #

library(tidyverse)
library(lme4)
library(lmerTest)
library(lubridate)
data("phen_dataset")
data("phen96_dataset")
data("phen97_dataset")
theme_set(theme_bw())

# # -- How much does "peak" flowering day vary by year
# what is peak flowering day? 
# note to self - tibbles don't work with EP packages...(grrrrrrr)
library(echPhenology)

# getting phenology summary data for all years
years <- phen_19967 %>% select(year) %>% distinct() %>% unlist() %>% unname() %>% as.character(.)
output <- list()
for(i in years){
  x <- phen_19967 %>% filter(year == i) %>% as.data.frame(.)
  FS <- makeFS(x, startCol = "startDtEarly", endCol = "endDtLate")
  output[[i]] <- summaryFS(FS)
}

#extracting peak dates for each year
peak_dates <- output %>%
  map_df("peakDate") %>%
  pivot_longer(everything(), names_to = "year", values_to = "peak_date") 
# summary data for peak date
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

# Bleh, I'm just going to have to do everything one by one...
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

plot_list = list()
for(i in years){
  peak_dt <- peak_dates %>% filter(year == i) %>% pull(peak_date)
  plot_list[[i]] <-  FS_df %>%
    filter(year == i) %>%
    mutate(mytitle = i) %>%
    ggplot()+
    geom_segment(aes(x = startDtEarly, xend = endDtLate,
                     y = ord, yend = ord), alpha = 0.5)+
    labs(x = NULL, y = NULL)+
    geom_vline(xintercept = peak_dt)+
    geom_point(data = count_dates %>% filter(year == i),
               aes(day, count))+
    scale_x_date(limits = c(as.Date(paste(i,"-06-19", sep = "")), as.Date(paste(i,"-09-01", sep = ""))))+
    scale_y_continuous(limits = c(0,300))+
    facet_grid(. ~ mytitle)+
    theme(strip.text       = element_text(size = rel(1.25)),
          strip.background = element_rect(color = "black",
                                          fill = "white"))
}

library(patchwork)
# figure out a more efficient way to do this
all_FS <- wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]],
           plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]])
ggsave("flowering_schedules.png", plot = all_FS)

iris %>% 
  rowid_to_column(var = "specimen") %>% 
  #  gather(contains("."), key = "Measurement", value = "cm") %>% 
  group_by(Species) %>% 
  nest %>% 
  mutate(plot = map2(Species, data,  function (.x,.y){
    ggplot(data = .y, aes(x =  Sepal.Length, y = Sepal.Width)) +
      geom_smooth() +
      ggtitle(label = .x)
    
  })) -> iris.with.plots

my_patchwork = function(data, plot_col="plot") {
  
  reduce(data[[plot_col]], `%+%`)
  
}

iris.with.plots %>% 
  my_patchwork()

      
# Average duration per flowering head; how much does that vary per year



# # -- distribution of head count in each year (in 1996 vs 1997 garden). Evenly spread? 

# histogram - ok, but different numbers of plants so hard to see differences in distributions
phen_19967 %>%
  mutate(expNm = as_factor(as.character(expNm)))%>%
  ggplot(aes(headCt, color = expNm))+
  geom_histogram(position = "identity", fill = "white", alpha = 0.5)+
  facet_wrap(~year)
#boxplot/count plot - it does look like the headCount seems to be similarly distributed between each 
#garden and each year. 1996 has a higher overall n each year (I believe), so it usually has 1 to 2 plants
# that have more heads than in 1997 garden. Sampling effect?
phen_19967 %>%
  mutate(expNm = as_factor(as.character(expNm)))%>%
  ggplot(aes(expNm, headCt))+
  geom_boxplot()+
  facet_wrap(~year)
phen_19967 %>%
  mutate(expNm = as_factor(as.character(expNm)))%>%
  ggplot(aes(expNm, headCt))+
  geom_count(alpha = 0.2)+
  facet_wrap(~year)

# # -- Do burn years have a longer flowering duration than non-burn years? Answer = no but
# head count is associated with duration, which is expected. 

burn_years <- c(2006, 2008, 2011, 2013, 2015)

phen_all <- phen_19967 %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "not_burned"))

m1 <- lm(dur ~ burn:headCt, data = phen_all)
summary(m1)
m2 <- lmerTest::lmer(dur ~ burn + headCt + burn*headCt + (1|cgPlaId), data = phen_all)
summary(m2)
anova(m2)
m3 <- lmerTest::lmer(dur ~ burn + headCt + (1|cgPlaId), data = phen_all)
anova(m3)
anova(m3, m2) # no difference, shoudl use m3 without interaction term. 

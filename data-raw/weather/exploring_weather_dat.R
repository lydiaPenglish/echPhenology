### Looking at relationship temperature/GDD and flowering time 

library(tidyverse)
data("phen_dataset")
theme_set(theme_bw())
burn_years <- c(2006, 2008, 2011, 2013, 2015)

# Data from NOAA - max temp, min temp, precip
weather <- read_csv("data-raw/weather/weather_dat.csv") %>%
  dplyr::rename_all(tolower) %>%
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  select(-c(station, name))

# trying to look at spring temperatures

weather_avg <- 
  weather %>%
  mutate(month = lubridate::month(date, label = TRUE),
         year  = lubridate::year(date),
         tavg  = (tmax + tmin)/2,
         year = as.factor(year)) %>%
  group_by(month, year) %>%
  summarize(tavg = mean(tavg, na.rm = TRUE),
            psum = sum(prcp, na.rm = TRUE)) %>%
  filter(month %in% c("Mar", "Apr", "May", "Jun")) %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "nb")) 

# join with phenology data - min FFD and median FFD
phen_weather <- 
  phen_19967 %>%
  group_by(year) %>%
  summarize(mean_FFD = mean(startNum),
            min_FFD  = min(startNum)) %>%
  left_join(weather_avg, ., by = "year") %>%
  # will add in later
  filter(year != 2019 & year != 2018) %>%
  arrange(year)

pw_wider <- phen_weather %>%
  pivot_wider(names_from = month, values_from = c(tavg, psum))

# quick modeling...

# i. Looking at when plant's first start flowering

# March temps
mar1 <- lm(mean_FFD ~ tavg_Mar, pw_wider)
summary(mar1)  # March temp seems to matter & higher March temps = earlier flowering

pw_wider %>%
  ggplot(aes(tavg_Mar, mean_FFD))+
  geom_point(aes(color = burn))

mar2 <- lm(tavg_Mar ~ burn, pw_wider) 
summary(mar2)

# April temps
apr1 <- lm(mean_FFD ~ tavg_Apr, data = pw_wider)
summary(apr1)   # April matters a little less

pw_wider %>%
  ggplot(aes(tavg_Apr, mean_FFD))+
  geom_point(aes(color = burn))

apr2 <- lm(tavg_Apr ~ burn, data = pw_wider)
summary(apr2)   # no sig difference in apr mean temps

# May temps 
may1 <- lm(mean_FFD ~ tavg_May, data = pw_wider)
summary(may1)

pw_wider %>%
  ggplot(aes(tavg_May, mean_FFD))+
  geom_point(aes(color = burn))

may2 <- lm(tavg_May ~ burn, data = pw_wider)
summary(may2) # no sig difference in May temps between burned and unburned years

# look at precip

mar_p <- lm(mean_FFD ~ psum_Jun, pw_wider)
summary(mar_p)

pw_wider %>%
  ggplot(aes(tot_prcp_May, mean_FFD))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

# Should I put them all together...?

all_try <- lm(mean_FFD ~ avg_avg_May, pw_wider)
summary(all_try)

# mixed model 

long_weather <- left_join(phen_19967, pw_wider, by = "year") %>%
  mutate(year_n = as.numeric(year))

mm1 <- lmerTest::lmer(startNum ~ psum_Mar + psum_Apr + psum_May + 
                        tavg_Mar + tavg_Apr + 
                        tavg_Jun +
              headCt +
              (1|cgPlaId), long_weather)
summary(mm1)
anova(mm1)
performance::r2(mm1)
lmerTest::rand(mm1)
performance::check_model(mm1)
performance::check_collinearity(mm1)

mm2 <- lmerTest::lmer(startNum ~ headCt + psum_Mar + psum_Jun + year +
                        (1|cgPlaId), long_weather)
summary(mm2)
anova(mm2)

library(nlme)
n1 <- lme(startNum ~ psum_Mar + psum_Apr + psum_May + 
            tavg_Mar + tavg_Apr + 
            tavg_Jun +
            headCt, random = ~1|cgPlaId, 
          data = long_weather)
n2 <- lme(startNum ~ psum_Mar + psum_Apr + psum_May + 
                  tavg_Mar + tavg_Apr + 
                  tavg_Jun +
                  headCt, random = ~1|cgPlaId, 
                  data = long_weather, correlation = corAR1())
anova(n2, n3)
performance::compare_performance(n1, n2)
ggResidpanel::resid_panel(n2)

# Morris extension weather data...growing degree day 

readxl::read_xlsx("data-raw/weather/morris_extension/2005.xlsx", sheet = 1, skip = 2)

read_weather_xlsx = function(f, into) {
  readxl::read_xlsx(f, sheet = 1, skip = 2) %>%
    dplyr::mutate(file=f) %>%
    tidyr::separate(file, into) 
}

read_weather_dir = function(path, pattern, into) {
  files = list.files(path = path,
                     pattern = pattern,
                     recursive = TRUE,
                     full.names = TRUE)
  plyr::ldply(files, read_weather_xlsx, into = into)
}

gdds <- 
  read_weather_dir(path = "data-raw/weather/morris_extension",
                   pattern = "*xlsx",
                   into = c("weather", "year", "extension")) %>%
  filter(`50` != "--" & `40` != "--") %>%
  mutate(GDD_50 = as.numeric(`50`),
         GDD_40 = as.numeric(`40`),
         date   = as.Date(Date)) %>%
  select(date, GDD_50, GDD_40) %>%
  mutate(month = lubridate::month(date, label = TRUE),
         year  = as.factor(lubridate::year(date))) %>%
  filter(month %in% c("Apr", "May", "Jun")) 

gd_nums <- gdds %>%
  filter(GDD_50 != 0) %>%
  group_by(year) %>%
  mutate(min_gdd = min(GDD_50),
         min_date = min(date)) %>%
  select(year, min_gdd, min_date) %>%
  distinct() %>%
  mutate(min_num = lubridate::yday(min_date))

phen_gdd <- 
  phen_19967 %>%
  group_by(year) %>%
  summarize(mean_FFD = mean(startNum),
            min_FFD  = min(startNum)) %>%
  left_join(gd_nums, by = "year")


phen_gdd %>%
  ggplot(aes(min_num, mean_FFD))+
  geom_point()


tt <-
  gdd_plant_dat %>%
  select(date, GDD_50, GDD_40, month, year, medianDate) %>%
  filter(date == medianDate) %>%
  distinct() %>%
  mutate(burn  = if_else(year %in% burn_years, "burn", "not_burned"))

t_mod <- lm(GDD_50 ~ burn, tt)  
anova(t_mod)
summary(t_mod)



# need to make flowering schedule plot again.. (I think?)

# median FFD

phen_weather %>%
  group_by(year) %>%
  summarize(medStartNum = median(startNum)) %>%
  mutate(medFFD = lubridate::as_date(medStartNum))

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

#extracting number of flowering plants on each date for using in graphics
count_dates <- output %>%
  map_df("fl.density") %>%
  mutate(year = lubridate::year(day))

# making into dataframe
FS_df <- phen_weather %>% 
  left_join(., peak_dates, by = "year") %>%
  mutate(stDayMonth  = format(as.Date(startDtEarly), format = "%m-%d"),
         endDayMonth = format(as.Date(endDtLate),    format = "%m-%d")) %>%
  group_by(year) %>%
  mutate(ord = row_number(startDtEarly)) %>%
  ungroup() %>%
  dplyr::select(year, cgPlaId, startDtEarly, endDtLate, yrPlanted, peak_date, stDayMonth:ord, tmax, prcp, burn)

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

p2005_temp <- 
  FS_df %>%
  filter(year == 2005) %>%
  ggplot()+
  geom_line(aes(x = startDtEarly, y = prcp), color = "red", size = 2)+
  geom_vline(xintercept = as.Date("2005-07-11"))+
  scale_x_date(limits = c(as.Date("2005-06-19"), as.Date("2005-09-01")))+
  labs(x = NULL, y = "Precip")+
  theme(axis.text.x = element_blank())

p2005_temp / p2005

p2008 <- FS_df %>%
  filter(year == 2008) %>%
  ggplot()+
  geom_segment(aes(x = startDtEarly, xend = endDtLate,
                   y = ord, yend = ord))+
  labs(x = NULL, y = NULL)+
  geom_vline(xintercept = as.Date("2008-07-26"))+
  geom_point(data = count_dates %>% filter(year == 2008),
             aes(day, count))+
  scale_x_date(limits = c(as.Date("2008-06-19"), as.Date("2008-09-01")))+
  scale_y_continuous(limits = c(0,300))

p2008_temp <- 
  FS_df %>%
  filter(year == 2008) %>%
  ggplot()+
  geom_line(aes(x = startDtEarly, y = prcp), color = "red", size = 2)+
  geom_vline(xintercept = as.Date("2008-07-26"))+
  scale_x_date(limits = c(as.Date("2008-06-19"), as.Date("2008-09-01")))+
  labs(x = NULL, y = "Maximum daily temperature")+
  theme(axis.text.x = element_blank())

p2008_temp / p2008



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
    geom_line(aes(x = startDtEarly, y = tmax/.35), color = "red", size = 2)+
    scale_x_date(limits = c(as.Date(paste(i,"-06-19", sep = "")), as.Date(paste(i,"-09-01", sep = ""))))+
    scale_y_continuous(limits = c(0,300),
                       sec.axis = sec_axis(~.*.35, name = "Maximum Daily Temp"))+
    facet_grid(. ~ mytitle)+
    theme(strip.background = element_rect(color = "black",
                                          fill = "white"),
          axis.text.x      = element_text(angle = 45, hjust = 1))
}

library(patchwork)
# figure out a more efficient way to do this....
all_FS <- wrap_plots(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]],
                     plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]], plot_list[[12]], plot_list[[13]])
all_FS



weather_avg %>%
  ggplot(aes(year, min_low))+
  geom_point(aes(color = burn))+
  facet_wrap(~month)+
  guides(color = FALSE)


g1 <- lm(startNum ~ year, data = phen_19967)
summary(g1)

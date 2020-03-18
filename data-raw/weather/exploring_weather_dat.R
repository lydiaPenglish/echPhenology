### Looking at relationship between maximum daily temperature and FFD 

library(tidyverse)
data("phen_dataset")
theme_set(theme_bw())
burn_years <- c(2006, 2008, 2011, 2013, 2015)

weather <- read_csv("data-raw/weather/weather_dat.csv") %>%
  dplyr::rename_all(tolower) %>%
  mutate(date = as.Date(date, "%m/%d/%Y")) %>%
  select(-c(station, name))

phen_weather <- 
  left_join(phen_19967, weather, by = c("startDtEarly" = "date")) %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "nb"))


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

# trying to look at spring temperatures

weather_avg <- 
  weather %>%
  mutate(month = lubridate::month(date, label = TRUE),
         year  = lubridate::year(date),
         tavg  = (tmax + tmin)/2) %>%
  group_by(month, year) %>%
  summarize(avg_avg  = mean(tavg, na.rm = TRUE),
            max_high = max(tmax, na.rm = TRUE),
            min_low  = min(tmin, na.rm = TRUE)) %>%
  filter(month %in% c("Mar", "Apr", "May", "Jun")) %>%
  mutate(burn = if_else(year %in% burn_years, "burned", "nb")) %>%
  # will add in later
  filter(year != 2019 & year != 2018)

weather_avg %>%
  ggplot(aes(year, min_low))+
  geom_point(aes(color = burn))+
  facet_wrap(~month)+
  guides(color = FALSE)


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
  select(date, GDD_50, GDD_40)

gdd_plant_dat <- gdds %>%
  mutate(month = lubridate::month(date, label = TRUE),
         year  = as.factor(lubridate::year(date))) %>%
  filter(month %in% c("May", "Jun", "Jul", "Aug")) %>%
  left_join(phen_19967, by = c("date" = "startDtEarly", "year"))

tt <-
  gdd_plant_dat %>%
  select(date, GDD_50, GDD_40, month, year, medianDate) %>%
  filter(date == medianDate) %>%
  distinct() %>%
  mutate(burn  = if_else(year %in% burn_years, "burn", "not_burned"))

t_mod <- lm(GDD_50 ~ burn, tt)  
anova(t_mod)
summary(t_mod)

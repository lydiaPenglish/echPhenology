# # Looking at the correlations in standardized flowering dates between different years # # 

# load packages and dataset
library(tidyverse)
library(psych)
data("phen_dataset")

# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# correlatios between (standardized) FFDs of different years

# long matrix of cgPlaId x Year and populated with DAM
cor_tab <- phen_19967 %>%
  select(cgPlaId, year, dam) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = dam) %>%
  column_to_rownames("cgPlaId")

cor_tab_ffd <- phen_19967 %>%
  select(cgPlaId, year, startNum) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = startNum) %>%
  column_to_rownames("cgPlaId")

# Going to use pearson correlation 
c1_holm  <- corr.test(cor_tab, use = "pairwise", method = "pearson")
ffd_holm  <- corr.test(cor_tab_ffd, use = "pairwise", method = "pearson")
r_vals <- round(c1_holm$r, 3)
r_vals_ffd <- round(ffd_holm$r, 3)
p_vals <- round(c1_holm$p, 4) # this doesn't look great.... 
p_vals_ffd <- round(ffd_holm$p, 4)
n_vals <- c1_holm$n
n_vals_ffd <- ffd_holm$n
 
# just want half of matrices
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

lower_r <- get_lower_tri(r_vals)
lower_n <- get_lower_tri(n_vals)
lower_p <- get_lower_tri(p_vals)
melt_r <- reshape2::melt(lower_r) %>%
  rename("r" = "value")
melt_n <- reshape2::melt(lower_n)%>%
  rename("n" = "value")
melt_p <- reshape2::melt(lower_p)%>%
  rename("p" = "value")

lower_r_ffd <- get_lower_tri(r_vals_ffd)
lower_n_ffd <- get_lower_tri(n_vals_ffd)
upper_n_ffd <- get_upper_tri(n_vals_ffd)
lower_p_ffd <- get_lower_tri(p_vals_ffd)
upper_p_ffd <- get_upper_tri(p_vals_ffd)

melt_r_ffd <- reshape2::melt(lower_r_ffd) %>%
  rename("r" = "value")
melt_n_ffd2 <- reshape2::melt(upper_n_ffd)%>%
  rename("n" = "value")
melt_n_ffd <- reshape2::melt(lower_n_ffd)%>%
  rename("n" = "value")
melt_p_ffd <- reshape2::melt(lower_p_ffd)%>%
  rename("p" = "value")
melt_p_ffd2 <- reshape2::melt(upper_p_ffd)%>%
  rename("p" = "value")

# table with everything...
cor_info_ffd <- left_join(melt_n_ffd, melt_r_ffd, by = c("Var1", "Var2")) %>%
  left_join(., melt_p_ffd, by = c("Var1", "Var2")) %>%
  left_join(., melt_n_ffd2, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)

# visualization version 1
c1 <- cor_info_ffd %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns")) %>%
  ggplot(aes(Var1, Var2))+
  geom_point(aes(size = r, alpha = samp_thres, color = p_value)) +
  geom_abline(slope = 1, intercept = nlevels(cor_info_ffd$Var1))+
  scale_y_discrete("", limits = levels(cor_info_ffd$Var1)) +
  scale_color_manual(values = c(my_cols[12], my_cols[10]))+
  scale_size_area(breaks = c(0, 0.1, 0.2, 0.4, 0.8))+
  labs(x = NULL, y = NULL, size = "Pearson R")+
  guides(alpha = FALSE, color = FALSE)+
  theme(legend.justification = c(-1,-2), legend.position=c(0,0),
        legend.background    = element_rect(color = "black"),
        legend.text          = element_text(size = rel(1.25)),
        legend.title         = element_text(size = rel(1.3)),
        axis.text            = element_text(size = rel(1.1)))
c1

# visualization version 2

cor_viz <- cor_info_ffd %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n.x < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns"),
         n.y        = replace(n.y, Var1 == Var2, NA))
cor_viz_sub <- cor_viz %>%
  filter(samp_thres == "yes")
  
c2 <- cor_viz %>%
  ggplot(aes(Var1, Var2))+
  geom_point(aes(size = r, color = p_value), pch = 21) +
  geom_point(data = cor_viz_sub, aes(size = r, fill = p_value), pch = 21)+
  geom_text(aes(label = n.y))+
  geom_abline(slope = 1, intercept = nlevels(cor_info_ffd$Var1))+
  scale_y_discrete("", limits = levels(cor_info_ffd$Var1)) +
  scale_color_manual(values = c(my_cols[12], my_cols[2]))+
  scale_fill_manual(values = c(my_cols[12], my_cols[2]))+
  scale_size(breaks = c(0, 0.1, 0.2, 0.4, 0.8), range = c(0,10))+
  labs(x = NULL, y = NULL, size = "Pearson R")+
  guides(color = FALSE, fill = FALSE)+
  theme(legend.background    = element_rect(color = "black"),
        legend.text          = element_text(size = rel(1.25)),
        legend.title         = element_text(size = rel(1.3)),
        axis.text            = element_text(size = rel(1.1)))
c2

# table with number times years were correlated with other years

year <- row.names(p_vals_ffd)

tab <- p_vals_ffd %>%
  as_tibble(.) %>%
  mutate_if(is.numeric, ~if_else(. <= 0.05, 1, 0)) %>%
  mutate(times_corr = rowSums(.))



tab <- cor_viz %>%
  filter(!(Var1 == Var2) & !(is.na(p_value))) %>%
  select(Var1, Var2, p_value) 
         
         
         & !(is.na(p_value))) %>%
  group_by(Var1, p_value) %>%
  summarize(n = n())

# correlation of duration 

cor_tab_dur <- phen_19967 %>%
  select(cgPlaId, year, dur) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = dur) %>%
  column_to_rownames("cgPlaId")

dur_holm  <- corr.test(cor_tab_dur, use = "pairwise", method = "pearson")
r_vals_dur <- round(dur_holm$r, 3)
p_vals_dur <- round(dur_holm$p, 3) 
n_vals_dur <- dur_holm$n

lower_r_dur <- get_lower_tri(r_vals_dur)
lower_n_dur <- get_lower_tri(n_vals_dur)
lower_p_dur <- get_lower_tri(p_vals_dur)
melt_r_dur <- reshape2::melt(lower_r_dur) %>%
  rename("r" = "value")
melt_n_dur <- reshape2::melt(lower_n_dur)%>%
  rename("n" = "value")
melt_p_dur <- reshape2::melt(lower_p_dur)%>%
  rename("p" = "value")

# table with everything...
cor_info_dur <- left_join(melt_n_dur, melt_r_dur, by = c("Var1", "Var2")) %>%
  left_join(., melt_p_dur, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)


d1 <- cor_info_dur %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns")) %>%
  ggplot(aes(Var1, Var2))+
  geom_point(aes(size = r, alpha = samp_thres, color = p_value)) +
  geom_abline(slope = 1, intercept = nlevels(cor_info_dur$Var1))+
  scale_y_discrete("", limits = levels(cor_info_dur$Var1)) +
  scale_color_manual(values = c(my_cols[12], my_cols[10]))+
  scale_size_area(breaks = c(0, 0.1, 0.2, 0.4, 0.8))+
  labs(x = NULL, y = NULL, size = "Pearson R")+
  guides(alpha = FALSE, color = FALSE)+
  theme(legend.justification = c(-1,-2), legend.position=c(0,0),
        legend.background    = element_rect(color = "black"),
        legend.text          = element_text(size = rel(1.25)),
        legend.title         = element_text(size = rel(1.3)),
        axis.text            = element_text(size = rel(1.1)))
d1

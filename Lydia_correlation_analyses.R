# # Looking at the correlations in standardized flowering dates between different years # # 

# load packages and dataset
library(tidyverse)
library(psych)
library(patchwork)
data("phen_dataset")

# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# correlatios between FFDs of different years - had previously been standardizing 
# FFDs beforehand but not am not

cor_tab_ffd <- phen_19967 %>%
  select(cgPlaId, year, startNum) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = startNum) %>%
  column_to_rownames("cgPlaId")

# Doing both pearson and spearman - going to use spearman bc then it ranks
ffd_spear  <- corr.test(cor_tab_ffd, use = "pairwise", method = "spearman")
r_spear_ffd <- round(ffd_spear$r, 3)
p_spear_ffd <- round(ffd_spear$p, 4)
n_spear_ffd <- ffd_spear$n

ffd_pears  <- corr.test(cor_tab_ffd, use = "pairwise", method = "pearson")
r_pears_ffd <- round(ffd_pears$r, 3)
p_pears_ffd <- round(ffd_pears$p, 4)
n_pears_ffd <- ffd_pears$n

# functions to get just half of matrices - note need UPPER half for p-values, that
# is adjusted for multiple testing correction
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

lower_r_spear <- get_lower_tri(r_spear_ffd)
lower_n_spear <- get_lower_tri(n_spear_ffd)
upper_p_spear <- get_upper_tri(p_spear_ffd)
upper_r_spear <- get_upper_tri(r_spear_ffd)
upper_n_spear <- get_upper_tri(n_spear_ffd)

lower_r_pears <- get_lower_tri(r_pears_ffd)
lower_n_pears <- get_lower_tri(n_pears_ffd)
upper_p_pears <- get_upper_tri(p_pears_ffd)
upper_r_pears <- get_upper_tri(r_pears_ffd)
upper_n_pears <- get_upper_tri(n_pears_ffd)

# make into a dataframe 
melt_r_spear <- reshape2::melt(upper_r_spear) %>%
  rename("r" = "value")
melt_n_spear <- reshape2::melt(upper_n_spear)%>%
  rename("n" = "value")
melt_p_spear <- reshape2::melt(upper_p_spear)%>%
  rename("p" = "value")
melt_n_spear_low <- reshape2::melt(lower_n_spear)%>%
  rename("n" = "value")

melt_r_pears <- reshape2::melt(upper_r_pears) %>%
  rename("r" = "value")
melt_n_pears <- reshape2::melt(upper_n_pears)%>%
  rename("n" = "value")
melt_p_pears <- reshape2::melt(upper_p_pears)%>%
  rename("p" = "value")

# table with everything...
cor_ffd_spear <- left_join(melt_n_spear, melt_r_spear, by = c("Var1", "Var2")) %>%
  left_join(., melt_p_spear, by = c("Var1", "Var2")) %>%
  left_join(., melt_n_spear_low, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)

cor_ffd_pears <- left_join(melt_n_pears, melt_r_pears, by = c("Var1", "Var2")) %>%
  left_join(., melt_p_pears, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)

# visualization spearman
cor_viz_spear <- cor_ffd_spear %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n.x < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns"),
         n.y        = replace(n.y, Var1 == Var2, NA))
  
cor_viz_sub <- cor_viz_spear %>%
  filter(samp_thres == "yes") %>%
  filter(!(is.na(r)))

numb_cor <- data.frame(
  Var1 = c("2005", "2006", "2007", "2008", "2009", "2010", 
           "2011", "2012", "2013", "2014", "2015", "2016", 
           "2017"),
  Var2 = c("2005", "2006", "2007", "2008", "2009", "2010", 
           "2011", "2012", "2013", "2014", "2015", "2016", 
           "2017"),
  Num_cor = c(0, 4, 4, 3, 3, 3, 6, 0, 1, 1, 1, 0, "NA")
)

ffd_plot <-
  cor_viz_spear %>%
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = r, color = p_value), pch = 21) +
  geom_point(data = cor_viz_sub, aes(size = r, fill = p_value), pch = 21) +
  geom_abline(slope = 1, intercept = nlevels(cor_ffd_spear$Var1), lty = 2) +
  geom_label(data = numb_cor, aes(Var1, Var2, label = Num_cor)) +
  geom_text(aes(label = n.y), fontface = "italic", size = 3) +
  scale_y_discrete("", limits = levels(cor_ffd_spear$Var1)) +
  scale_color_grey(start = 0.7, end = 0.3) +
  scale_fill_grey(start = 0.7, end = 0.3) +
  scale_size_area(breaks = c(0, 0.1, 0.2, 0.4, 0.8), max_size = 10) +
  labs(x = NULL, y = NULL, size = "Spearman R") +
  guides(fill = FALSE, color = FALSE) +
  coord_fixed() +
  ggtitle("FFD") +
  theme(
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = rel(1.25)),
    legend.title = element_text(size = rel(1.3)),
    axis.text = element_text(size = rel(1.1)),
    plot.title = element_text(size = rel(1.7), hjust = 0.5)
  )
ffd_plot

scale_size(breaks = c(0, 0.1, 0.2, 0.4, 0.8), range = c(0,10))+

# wanting to plot without p-values
cor_viz_spear %>%
  ggplot(aes(Var1, Var2))+
  geom_point(aes(color = r), size = 7, pch = 21) +
  geom_point(data = cor_viz_sub, aes(fill = r), size = 7, pch = 21)+
  geom_text(aes(label = n.y), fontface = "italic", size =3) +
  scale_color_gradient(na.value = NA, breaks = c(-0.1, 0.9)) +
  scale_fill_gradient(na.value = NA, breaks = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.9)) +
  guides(color = FALSE)+
  labs(x = NULL,
       y = NULL)


# visualizing pearson
cor_viz_pears <- cor_ffd_pears %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns"))
cor_viz_sub_p <- cor_viz_pears %>%
  filter(samp_thres == "yes")

p_pears <- cor_viz_pears %>%
  ggplot(aes(Var1, Var2)) +
  geom_point(aes(size = r, color = p_value), pch = 21) +
  geom_point(data = cor_viz_sub_p, aes(size = r, fill = p_value), pch = 21) +
  geom_abline(slope = 1, intercept = nlevels(cor_ffd_pears$Var1)) +
  scale_y_discrete("", limits = levels(cor_ffd_pears$Var1)) +
  scale_color_manual(values = c(my_cols[12], my_cols[10])) +
  scale_fill_manual(values = c(my_cols[12], my_cols[10])) +
  guides(fill = FALSE, color = FALSE) +
  scale_size(breaks = c(0, 0.1, 0.2, 0.4, 0.8), range = c(0, 10)) +
  labs(x = NULL, y = NULL, size = "Pearson R") +
  theme(
    legend.background = element_rect(color = "black"),
    legend.text = element_text(size = rel(1.25)),
    legend.title = element_text(size = rel(1.3)),
    axis.text = element_text(size = rel(1.1))
  )
p_pears

p_spear + p_pears

# correlation of duration - also will use spearman 

cor_tab_dur <- phen_19967 %>%
  select(cgPlaId, year, dur) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = dur) %>%
  column_to_rownames("cgPlaId")

dur_holm  <- corr.test(cor_tab_dur, use = "pairwise", method = "spearman")
r_vals_dur <- round(dur_holm$r, 3)
p_vals_dur <- round(dur_holm$p, 3) 
n_vals_dur <- dur_holm$n

upper_r_dur <- get_upper_tri(r_vals_dur)
upper_n_dur <- get_upper_tri(n_vals_dur)
upper_p_dur <- get_upper_tri(p_vals_dur)
lower_n_dur <- get_lower_tri(n_vals_dur)

# make into a dataframe 
melt_r_dur <- reshape2::melt(upper_r_dur) %>%
  rename("r" = "value")
melt_n_dur <- reshape2::melt(upper_n_dur)%>%
  rename("n" = "value")
melt_p_dur <- reshape2::melt(upper_p_dur)%>%
  rename("p" = "value")
melt_n_dur_low <- reshape2::melt(lower_n_dur)%>%
  rename("n" = "value")
# table with everything...
cor_info_dur <- left_join(melt_n_dur, melt_r_dur, by = c("Var1", "Var2")) %>%
  left_join(., melt_p_dur, by = c("Var1", "Var2")) %>%
  left_join(., melt_n_dur_low, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)

# visualizing
cor_viz_dur <- cor_info_dur %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n.x < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns"),
         n.y        = replace(n.y, Var1 == Var2, NA))

cor_sub_dur <- cor_viz_dur %>%
  filter(samp_thres == "yes")

numb_dur_cor <- data.frame(
  Var1 = c("2005", "2006", "2007", "2008", "2009", "2010", 
           "2011", "2012", "2013", "2014", "2015", "2016", 
           "2017"),
  Var2 = c("2005", "2006", "2007", "2008", "2009", "2010", 
           "2011", "2012", "2013", "2014", "2015", "2016", 
           "2017"),
  Num_cor = c(0, 0, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, "NA")
)


dur_plot <- 
  cor_viz_dur %>%
  ggplot(aes(Var1, Var2))+
  # geoms
  geom_point(aes(size = r, color = p_value), pch = 21) +
  geom_point(data = cor_sub_dur, aes(size = r, fill = p_value), pch = 21)+
  geom_abline(slope = 1, intercept = nlevels(cor_info_dur$Var1))+
  geom_text(aes(label = n.y), fontface = "italic", size =3)+
  geom_label(data = numb_dur_cor, aes(Var1, Var2, label = Num_cor))+
  # scales
  scale_y_discrete("", limits = levels(cor_info_dur$Var1)) +
  scale_color_grey(start = 0.7, end = 0.3)+
  scale_fill_grey(start = 0.7, end = 0.3)+
  scale_size_area(breaks = c(0, 0.1, 0.2, 0.4, 0.8), max_size = 10)+
  labs(x = NULL, 
       y = NULL, 
       size = "Spearman R")+
  guides(fill  = FALSE, 
         color = FALSE,
         size  = FALSE)+
  coord_fixed()+
  ggtitle("Duration")+
  theme(legend.background    = element_rect(color = "black"),
        legend.text          = element_text(size = rel(1.25)),
        legend.title         = element_text(size = rel(1.3)),
        axis.text.x          = element_text(size = rel(1.4)),
        axis.text.y          = element_blank(),
        plot.title           = element_text(size = rel(1.7), hjust = 0.5))
dur_plot
ffd_plot + dur_plot + plot_layout(guides = 'collect')


legend.background = element_rect(color = "black"),
legend.text = element_text(size = rel(1.25)),
legend.title = element_text(size = rel(1.3)),
axis.text = element_text(size = rel(1.1)),
plot.title = element_text(size = rel(1.7), hjust = 0.5)

# ---- older stuff ----

# long matrix of cgPlaId x Year and populated with DAM
cor_tab <- phen_19967 %>%
  select(cgPlaId, year, dam) %>%
  #filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = dam) %>%
  column_to_rownames("cgPlaId")

c1_holm  <- corr.test(cor_tab, use = "pairwise", method = "pearson")
r_vals <- round(c1_holm$r, 3)
p_vals <- round(c1_holm$p, 4) # this doesn't look great.... 
n_vals <- c1_holm$n
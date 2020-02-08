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

# Going to use pearson correlation 
c1_holm  <- corr.test(cor_tab, use = "pairwise", method = "pearson")
r_vals <- round(c1_holm$r, 3)
p_vals <- round(c1_holm$p, 3) # this doesn't look great.... 
n_vals <- c1_holm$n
 
# just want half of matrices
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
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

# table with everything...
cor_info <- left_join(melt_n, melt_r, by = c("Var1", "Var2")) %>%
  left_join(., melt_p, by = c("Var1", "Var2")) %>%
  mutate_at(vars(Var1, Var2), as.character)

# visualization 
cor_info %>%
  # creating extra columns for labeling variables
  mutate(r          = na_if(r, 1),
         samp_thres = if_else(n < 10 , "no", "yes"),
         p_value    = if_else(p <= 0.05, "sig", "ns")) %>%
  ggplot(aes(Var1, Var2))+
  geom_point(aes(size = r, alpha = samp_thres, color = p_value)) +
  geom_abline(slope = 1, intercept = nlevels(cor_info$Var1))+
  scale_y_discrete("", limits = levels(cor_info$Var1)) +
  scale_color_manual(values = c(my_cols[12], my_cols[10]))+
  scale_size_area(breaks = c(0, 0.1, 0.2, 0.4, 0.8))+
  labs(x = NULL, y = NULL, size = "Pearson R")+
  guides(alpha = FALSE, color = FALSE)+
  theme(legend.justification = c(-1,-2), legend.position=c(0,0),
        legend.background    = element_rect(color = "black"),
        legend.text          = element_text(size = rel(1.25)),
        legend.title         = element_text(size = rel(1.3)),
        axis.text            = element_text(size = rel(1.1)))







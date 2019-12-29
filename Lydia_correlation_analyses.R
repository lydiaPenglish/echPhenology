# # Looking at the correlations in standardized flowering dates between different years # # 

# load packages and dataset
library(tidyverse)
data("phen_dataset")

# correlatios between (standardized) FFDs of different years

# long table of cgPlaId x Year and populated with DAM
cor_tab <- phen_19967 %>%
  select(cgPlaId, year, dam) %>%
  filter(year != "2016" & year != "2017") %>%
  pivot_wider(names_from = year, values_from = dam) %>%
  column_to_rownames("cgPlaId")

# pearson vs spearman
cor_sp   <- round(cor(cor_tab, use = "pairwise.complete.obs", method = "spearman"), 2)
cor_pear <- round(cor(cor_tab, use = "pairwise.complete.obs", method = "pearson"), 2)

library(psych)
c1_holm  <- corr.test(cor_tab, use = "pairwise", method = "spearman")
c1_unadj <- corr.test(cor_tab, use = "pairwise", method = "spearman", adjust = "none")
print(c1_holm) # after adjusting, this doesn't look great. 


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
up_cormat <- get_upper_tri(cormat)
up_melt <- reshape2::melt(up_cormat, na.rm = TRUE)

corP <- round(cor.test(cor_tab, use = "pairwise.complete.obs", method = "spearman"), 2)







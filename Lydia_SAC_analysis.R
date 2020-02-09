# # Analysis to understand if there is any spatial autocorrelation in the data # #
# # In other words, is flowering time and duration associated with position in the CG? 
# #     Are earlier flowering plants located in a certain section while later flowering plants 
# #     in a different location?
# # I think the best way to do this is to run a separate model for each year and assess if any years have SAC...

library(tidyverse)
library(ape)
data("phen_dataset")
# General graphics parameters
theme_set(theme_bw())
my_cols <-c("#E6AB02", "#D95F02", "#74C476","#238B45", "#00441B", 
            "#6BAED6", "#08519C", "#D0D1E6", "#7570B3", "#F781BF",
            "#A6761D","#666666")
scales::show_col(my_cols) # To see colors 

# will use Moran's I from ape package - calculates using Euclidean distances

# example with one year (2014)
p14 <- phen_19967 %>%
  filter(year == 2014) %>%
  select(cgPlaId, row, pos) %>%
  column_to_rownames("cgPlaId")
phen14 <- phen_19967 %>%
  filter(year == 2014)
cgdists <- as.matrix(dist(as.matrix(p14))) # euclidean distance matrix 
cgdist_inv <- 1/cgdists # take the inverse
diag(cgdist_inv) <- 0
m2 <- ape::Moran.I(phen14$startNum, cgdist_inv, alternative = 'greater')
print(m2)

# Now going to run the same thing on all years using a loop
years <- phen_19967 %>% select(year) %>% distinct() %>% unlist() %>% unname() %>% as.character(.)
output <- list()
for(i in years){
  p <- phen_19967 %>%
    filter(year == i) %>%
    select(cgPlaId, row, pos) %>%
    column_to_rownames("cgPlaId")
  pp <- phen_19967 %>%
    filter(year == i)
  cgdists <- as.matrix(dist(as.matrix(p))) # euclidean distance matrix 
  cgdist_inv <- 1/cgdists # take the inverse
  diag(cgdist_inv) <- 0
  output[[i]] <-  ape::Moran.I(pp$startNum, cgdist_inv, alternative = 'greater')
}

# extracting output from the list and calculating each years z-score
moran_output <- output %>%
  transpose() %>%
  as_tibble() %>%
  unnest(cols = c(observed, expected, sd, p.value)) %>%
  add_column(years, .before = TRUE) %>%
  mutate(z_score = (observed-expected)/sd,
         sig     = if_else(p.value < 0.05, "yes", "no"))
# result there is Spatial autocorrelation in some, but not all, years. (roughly half)

# graph of z-scores
moran_z <- moran_output %>%
  ggplot(aes(years, z_score))+
  geom_point(aes(color = sig), size = 4)+
  geom_hline(yintercept = 0, lty = 2)+
  labs(y = "Z-score", x = NULL,
       color = "Significant at 0.05?")+
  theme_bw()+
  scale_colour_manual(values = c(my_cols[12], my_cols[10]))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color = "black"),
        axis.title.y = element_text(size = rel(1.3)),
        axis.text    = element_text(size = rel(1.25)),
        legend.text  = element_text(size = rel(1.25)),
        legend.title = element_text(size = rel(1.3)))
moran_z
#ggsave("moran_z-scores.png", plot = moran_z)

# ----- doing the same thing for flowering duration ----- # 

# example with one year (2014)
m3 <- ape::Moran.I(phen14$dur, cgdist_inv, alternative = 'greater')
print(m3)

# Now going to run the same thing on all years using a loop
output2 <- list()
for(i in years){
  p <- phen_19967 %>%
    filter(year == i) %>%
    select(cgPlaId, row, pos) %>%
    column_to_rownames("cgPlaId")
  pp <- phen_19967 %>%
    filter(year == i)
  cgdists <- as.matrix(dist(as.matrix(p))) # euclidean distance matrix 
  cgdist_inv <- 1/cgdists # take the inverse
  diag(cgdist_inv) <- 0
  output2[[i]] <-  ape::Moran.I(pp$dur, cgdist_inv, alternative = 'greater')
}

# extracting output from the list and calculating each years z-score
moran_output2 <- output2 %>%
  transpose() %>%
  as_tibble() %>%
  unnest(cols = c(observed, expected, sd, p.value)) %>%
  add_column(years, .before = TRUE) %>%
  mutate(z_score = (observed-expected)/sd,
         sig     = if_else(p.value < 0.05, "yes", "no"))
# result: 


# ----- old method of nesting dataframe ------
# isn't reliable!

# nested data.frame to calculate all the Moran's
nest_phen <- phen_19967 %>%
  remove_rownames()%>%
  dplyr::select(year, cgPlaId, row, pos, startNum) %>%
  nest(data = c(cgPlaId, startNum, row, pos)) %>%
  mutate(startNums = purrr::map(data, ~ select(., startNum) %>% 
                                  unlist(.) %>% unname(.))) %>%
  mutate(data = purrr::map(data, ~ select(., -startNum) %>% column_to_rownames(., var = "cgPlaId")
                           %>% as.matrix(.))) %>%
  mutate(dists = purrr::map(data, ~ as.matrix(1/dist(.)))) %>%
  mutate(morans = purrr::map2(startNums, dists,
                              ~ Moran.I(.x, .y)))

nest_phen$morans[[10]] 
print(m2) ### YAY they look the same!

# extracting all the Moran info 
nest_phen$morans

# some of these are significant and some of them are not....now what...?




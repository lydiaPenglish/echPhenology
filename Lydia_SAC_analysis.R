# # Analysis to understand if there is any spatial autocorrelation in the data # #
# # In other words, is flowering time associated with position in the CG? 
# #     Are earlier flowering plants located in a certain section while later flowering plants 
# #     in a different location?
# # I think the best way to do this is to run a separate model for each year and assess if any years have SAC

library(tidyverse)
library(ape)
library(vegan)
library(spdep)
data("phen_dataset")

# will use Moran's I from ape package - calculates using Euclidean distances

# example of how to do one....
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

years <- phen_19967 %>% select(year) %>% distinct() %>% unlist() %>% unname() %>% as.character(.)
output <- list()
# loop 
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

# result there is Spatial autocorrelation in some, but not all, years. (roughly half)

moran_output <- output %>%
  transpose() %>%
  as_tibble() %>%
  unnest(cols = c(observed, expected, sd, p.value)) %>%
  add_column(years, .before = TRUE) %>%
  mutate(z_score = (observed-expected)/sd,
         sig     = if_else(p.value < 0.05, "yes", "no"))

moran_output %>%
  ggplot(aes(years, z_score))+
  geom_point(aes(color = sig), size = 3)+
  geom_hline(yintercept = 0, lty = 2)+
  labs(y = "Z-score", x = NULL,
       color = "Significant at 0.05?")+
  theme_bw()+
  scale_colour_manual(values = c("black", "red"))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color = "black"))


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




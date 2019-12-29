# # Analysis to understand if there is any spatial autocorrelation in the data # #
# # In other words, is flowering time associated with position in the CG? 
# #     Are earlier flowering plants located in a certain section while later flowering plants 
# #     in a different location?
# # I think the best way to do this is to run a separate model for each year and assess if any years have SAC

library(tidyverse)
library(ape)
library(vegan)
data("phen_dataset")

# try with mantel for a given year
# testing with 2011 - NS. Could make a function to do this for every year...
p14 <- phen_19967 %>%
  filter(year == 2014) %>%
  select(cgPlaId, row, pos) %>%
  column_to_rownames("cgPlaId")
phen14 <- phen_19967 %>%
  filter(year == 2014)
mantel(dist(p14), dist(phen_19967 %>%
         filter(year == 2014) %>%
         select(startNum)))

# Moran's I - two different ways. They seem to get similar ish p-values. 
W <-tri2nb(as.matrix(p13)) #weights with Delauney tesselation
m1 <- moran.test(phen13$startNum, nb2listw(W))
print(m1)

# nested data.frame
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
                              ~ Moran.I(.x, .y, alternative = 'greater')))
nest_phen$startNums[[10]]

nest_moran$morans # none of the years have significant Moran p-values....
# why don't they match this...?

as.matrix(p14)
    
cgdists <- as.matrix(dist(as.matrix(p14)))
cgdist_inv <- 1/cgdists
cgdist_inv
diag(cgdist_inv) <- 0
ape::Moran.I(phen14$startNum, cgdist_inv, alternative = 'greater')
print(m2)
nest_phen$startNums

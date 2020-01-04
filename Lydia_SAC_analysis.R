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
# version A - confusing and I'm not sure I understand
library(spdep)
p14[sample(nrow(p14)),]
p14r <- as.matrix(p14[sample(nrow(p14)),])

W <-tri2nb(p14r) #weights with Delauney tesselation
m1 <- moran.test(phen14$startNum, nb2listw(W))
print(m1)

# version B - easier to understand, using ape package
cgdists <- as.matrix(dist(as.matrix(p14)))
cgdist_inv <- 1/cgdists
diag(cgdist_inv) <- 0
m2 <- ape::Moran.I(phen14$startNum, cgdist_inv, alternative = 'greater')
print(m2)


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


p11 <- phen_19967 %>%
  filter(year == 2011) %>%
  select(cgPlaId, row, pos) %>%
  column_to_rownames("cgPlaId")
phen11 <- phen_19967 %>%
  filter(year == 2011)

# version B - easier to understand, using ape package
cgdists11 <- as.matrix(dist(as.matrix(p11)))
cgdist_inv11 <- 1/cgdists11
diag(cgdist_inv11) <- 0
m11 <- ape::Moran.I(phen11$startNum, cgdist_inv11)
print(m11)

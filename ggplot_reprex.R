library(reprex)



reprex({
  library(tidyverse)
  
  # I've set up this dataframe
  ex <- data.frame(x1 = c(rep("2010", 3), rep("2011", 3), rep("2012", 3)),
                   x2 = rep(c("2010", "2011", "2012"), 3),
                   r = c(NA, 0.5, 0.2, NA, NA, 0.1, NA, NA, NA),
                   n = c(NA, NA, NA, 25, NA, NA, 70, 50, NA))
  
  ex
  
  # when I run this graph, with size as an aesthetic, NAs are automatically removed
  ex %>%
      ggplot(aes(x1, x2))+
      geom_point(aes(size = r))+
      geom_text(aes(label = n))
  
  # when I run this graph, with color as an aesthetic, NAs are included
  ex %>%
      ggplot(aes(x1, x2))+
      geom_point(aes(color = r), size = 6)+
      geom_text(aes(label = n))
  }, venue = "so")

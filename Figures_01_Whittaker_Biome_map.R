library(plotbiomes)
whittaker_base_plot()
library(ggplot2)

MAT_MAP <- read.csv(file = 'D:/Data/Global Thermoregulation/BESS/siteInfo_MAT_MAP_v2.csv')
mydata <- MAT_MAP[, c(10,11)]

whittaker_base_plot() +
  # add the temperature - precipitation data points
  geom_point(data = mydata, 
             aes(x = MAT, 
                 y = MAP), 
             size   = 2.5,
             shape  = 21,
             colour = "gray95",
             fill   = "brown1",
             stroke = 0,
             alpha  = 0.7) +
  theme_bw()


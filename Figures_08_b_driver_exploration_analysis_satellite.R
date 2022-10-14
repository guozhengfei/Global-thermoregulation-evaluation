# 
climate_LAI_PFT <- read.csv(file = 'Y:/project 2/RF analysis for global dT variation/dT_drivers_Satellite.csv')

mydata <- climate_LAI_PFT[, c(1:10)]
# mydata$VPD <- 0.611*exp(17.502*mydata$Tair/(mydata$Tair +240.97))*(1-mydata$RH/100);

mydata$igbp <- factor(mydata$igbp)
mydata$igbp <- as.numeric(mydata$igbp)

# collinear test
mydata <- climate_LAI_PFT[, c(1,3,4,6,2,9,10,5)]
names(mydata)[names(mydata) == "igbp"] <- "PFT"
names(mydata)[names(mydata) == "Tair"] <- "Ta"

head(mydata)
cormat <- round(cor(mydata),2)
head(cormat)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
 
upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)

# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

# distribution
library(Hmisc)
hist.data.frame(mydata)

library(lme4)
library(performance)
mydata<-mydata[-c(20430, 22571, 22785), ]
mydata_scaled <- scale(mydata)
mydata_scaled <- as.data.frame(mydata_scaled)
head(mydata_scaled)
dT<- climate_LAI_PFT[, c(11)]
dT<-dT[-c(20430, 22571, 22785)]


m1 = glm(dT~ Ta + PAR_org + RH + rainfall + windspeed + LAI + PFT + elevation + PAR_org * RH +rainfall*LAI + PAR_org * LAI + Ta*LAI +
            Ta*PFT + PAR_org*PFT + rainfall*PFT + rainfall*PFT, data = mydata_scaled)
model_performance(m1)
step(m1)
plot(m1)

m1_2 =  glm(formula = dT ~ Ta + PAR_org + RH + rainfall + windspeed + 
              LAI + PFT + elevation + PAR_org:RH + rainfall:LAI + PAR_org:LAI + 
              Ta:LAI + Ta:PFT + PAR_org:PFT, data = mydata_scaled)
model_performance(m1_2)
summary(m1_2)

plot(m1_2)
hist(residuals(m1_2), breaks = seq(-7, 7, 0.4))

d1<-summary(m1_2)
d2<-d1$coefficients

write.csv(d2,'Y:/project 2/RF analysis for global dT variation/dT_drivers_Satellite_importance.csv')

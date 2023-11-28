## Script to calculate stress from calcite twin density following Ryback et al., 2013 paleopiezometer
## Created by Madison Frank, updated 2023/11/28

library(openxlsx)
library(plotly)
library(ggplot2)

### Load raw data ### -------------------------------------------------------------------------

# Read in raw twin count data and calculate 
ca_twin_raw <- read.xlsx(filepath = "BLM_Ca_TwinDensity.xlsx", sheet = "Ca_Twins", startRow = 2, colnames = T)

### Calculate differential stress ###----------------------------------------------------------
 
# Create a data frame holding the number of twins measured over a given distance. Only keep rows with data (complete cases)
twin_data <- ca_twin_raw[,c('twin_count', 'distance_um')]
twin_data <- twin_data[complete.cases(twin_data),]
  
colnames(twin_data) <- c("twins", "distance_um")
str(twin_data)
  
# Convert um to mm
twin_data$distance_mm <- twin_data$distance_um/1000
  
# Calculate twins/mm
twin_data$density <- twin_data$twins/twin_data$distance_mm
  
# Calculate differential stress following Rybacki et al., 2013 equation 7
twin_data$diffStress <- 19.5 * sqrt(twin_data$density)
  
twin_data$error_high <- (19.5 + 9.8) * sqrt(twin_data$density)
twin_data$error_low <- (19.5 - 9.8) * sqrt(twin_data$density)
  
### Generate plot ### -------------------------------------------------------------------------
# plot with ggplot for exporting
ggRyb <- ggplot(data = twin_data)+
  geom_point(aes(x = density, y = diffStress), size = 2)+
  geom_errorbar(aes(x = density, ymax = error_high, ymin = error_low), color = "grey", width = 1.5)+
  theme_light()

ggRyb  

# export plot
# ggsave(ggRyb, 
#        filename = "Rybacki_TwinDensity_Piez.svg", 
#        device = "svg")


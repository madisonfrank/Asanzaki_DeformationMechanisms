## Script to calculate stress from calcite twin density following Ryback et al., 2013 paleopiezometer
## Created by Madison Frank, updated 2023/05/08

# load in functions
source("Ca_Twinning_Functions.R")

library(plotly)
library(ggplot2)

### Load raw data ### -------------------------------------------------------------------------

# Read in raw twin count data and calculate 
twin_data <- getRawTwins(filepath = "Ca_TwinDensity.xlsx", 
                         sheet = "Ca_TwinDensity", startRow = 2, 
                         colnames = c('Thin.Section', 'Shear.localisation', 'Distance.(um)', 'Twin.Count'))

# change sample name
twin_data$name <- gsub(pattern = "ESZL-1", replacement = "SZLE-1", x = twin_data$name)


### Boxplots to visualise outliers ###---------------------------------------------------------

fig <- plot_ly(data = twin_data, x = ~name, y = ~density, type = "box", color = ~localisation, hoverinfo = 'y') %>%
  layout(title = 'Twin densities', 
         yaxis = list(title = 'Twin density (# of twins/mm)'))
fig

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


## Script to calculate stress from calcite twin density following Ryback et al., 2013 paleopiezometer
## Created by Madison Frank, updated 2023/05/08

# load in functions
source("Ca_Twinning_Functions.R")

library(plotly)
library(ggplot2)

### Load raw data ### -------------------------------------------------------------------------

# Read in raw twin count data and calculate 
twin_data <- getRawTwins(filepath = "C:/Users/madho/Documents/University/PhD/R Scripts/Asanzaki/Ca_TwinDensity.xlsx", 
                         sheet = "Ca_TwinDensity", startRow = 2, 
                         colnames = c('Thin.Section', 'Shear.localisation', 'Distance.(um)', 'Twin.Count'))

# change sample name
twin_data$name <- gsub(pattern = "ESZL-1", replacement = "SZLE-1", x = twin_data$name)


### Boxplots to visualise outliers ###---------------------------------------------------------

fig <- plot_ly(data = twin_data, x = ~name, y = ~density, type = "box", color = ~localisation, hoverinfo = 'y') %>%
  layout(title = 'Twin densities', 
         yaxis = list(title = 'Twin density (# of twins/mm)'))
fig

# Remove outliers
no_outliers <- removeTwinOutliers(twin_data)

fig_no_outliers <- plot_ly(data = no_outliers, x = ~name, y = ~density, type = "box", color = ~localisation, hoverinfo = 'y') %>%
  layout(title = 'Twin densities', 
         yaxis = list(title = 'Twin density (# of twins/mm)'))

fig_no_outliers


# Split twin dataset by thin section 
ts_twins <- split(no_outliers, no_outliers$name)

# colour by proportion of basalt in sample
colorPal = data.frame(name = unique(no_outliers$name), 
                      colors = c("#ffa07a", "#89CFF0", "#0000ff", "#6495ed", 
                                 "#ffdab9", "#00A300", "#ff0000"))
no_outliers <- merge(no_outliers, colorPal)

# plot cleaned dataset as box plots
plot_ly(data = no_outliers, x = ~density, y = ~name, type = 'box', showlegend = FALSE) %>%
  layout(title = "Twin Density", 
         yaxis = list(title = "Sample", standoff = 320), 
         xaxis = list(title = "Twin Density (twin no/mm)"), 
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10))


# plot differential stress vs twin density (Rybacki et al,, 2013 piezometer)
ryb <- plot_ly(data = no_outliers, x = ~density, y = ~diffStress, mode = 'markers', type = 'scatter', 
               marker = list(color = ~colors),
        name = ~name, error_y = ~list(type = 'data', array = error, color = colors, opacity = 20)) %>%
  layout(title = "Calcite Twin Piezometer (Rybacki et al., 2013)", 
         yaxis = list(title = "Differential Stress (MPa)", standoff = 120), 
         xaxis = list(title = "Density (twin no/mm)"), 
         font = list(size = 13), 
         margin = list(t = 50, b = 70, pad = 10))

ryb

# plot with ggplot for exporting
ggRyb <- ggplot(data = no_outliers)+
  geom_point(aes(x = density, y = diffStress), size = 2)+
  geom_errorbar(aes(x = density, ymax = diffStress+error, ymin = diffStress-error), color = "grey", width = 1.5)+
  theme_light()

ggRyb  

# export plot
# ggsave(ggRyb, 
#        filename = "C:/Users/madho/Documents/University/PhD/Manuscipts/Asanzaki_Manuscript/Inkscape Files/images/RybackiOutput1.svg", 
#        device = "svg")


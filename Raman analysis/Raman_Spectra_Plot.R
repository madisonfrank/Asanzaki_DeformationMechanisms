## Script to plot carbonaceous material raman spectra from PeakFit solution 
## Created by Madison Frank, updated 2023/06/08

library(SciViews)
library(ggplot2)
library(gridExtra)

# load spectra file
spectra <- read.csv("MS-1_001_Solution.csv", skip = 1)

# create plot of raw spectra data
ggRaman_raw <- ggplot(data = spectra)+
  geom_line(aes(x = X.Observed, y = Y.Observed), size = 1)+
  geom_line(aes(x = X.Generated, y = Baseline..Linear.Bg), color = "grey", linetype = 2, size = 1) +
  theme_minimal()

# create plot decomped spectra
ggRaman_mod <- ggplot(data = spectra)+
  geom_line(aes(x = X.Generated, y = Y.Generated - Baseline..Linear.Bg), color = "black", size = 1)+#, linetype = 2)+
  geom_line(aes(x = X.Generated, y = Peak.1..Gauss.Lor.Amp), color = "violet", size = 1)+
  geom_line(aes(x = X.Generated, y = Peak.2..Gauss.Lor.Amp), color = "red", size = 1) +
  geom_line(aes(x = X.Generated, y = Peak.3..Gauss.Lor.Amp), color = "orange", size = 1)+
  geom_line(aes(x = X.Generated, y = Peak.4..Gauss.Lor.Amp), color = "darkgreen", size = 1) +
  geom_line(aes(x = X.Generated, y = Peak.5..Gauss.Lor.Amp), color = "blue", size = 1) +
  theme_minimal()

# plot raw and decomp solution as subplots
ggRaman <- grid.arrange(ggRaman_raw, ggRaman_mod, nrow = 2)
ggRaman

# save image as svg
# ggsave(ggRaman,
#        filename = "RamanSpectraOutput.svg",
#        device = "svg")


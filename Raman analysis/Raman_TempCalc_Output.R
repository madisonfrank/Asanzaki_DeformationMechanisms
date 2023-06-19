## Script to calculate temperature from carbonaceous materia raman spectra following methods of Kouketsu et al., 2014
## Rewritten from CMcalculator_Kouketsu2014 matlab script
## Created by Madison Frank, updated 2023/06/08

library(plotly)

#---------------------------------------------------------------------------------------------------------------------------
# FUNCTION to calculate temperatures

getRamanTemp <- function(filepath, targetBand){
  
  # Average peak location of each band and +/- range to search within 
  ramanBands <- data.frame(band = c("G", "D1", "D2", "D3", "D4"),
                           loc = c(1580, 1350, 1620, 1510, 1245))
  bandRange <- 20
  
  # Find the location of the requested peak   
  targetRange <- switch (targetBand,
                       "G" = c(ramanBands[1,2]-bandRange, ramanBands[1,2]+bandRange), 
                       "D1" = c(ramanBands[2,2]-bandRange, ramanBands[2,2]+bandRange), 
                       "D2" = c(ramanBands[3,2]-bandRange, ramanBands[3,2]+bandRange), 
                       "D3" = c(ramanBands[4,2]-bandRange, ramanBands[4,2]+bandRange), 
                       "D4" = c(ramanBands[5,2]-bandRange, ramanBands[5,2]+bandRange))

  # Get list of csv file names
  ramanResults <- list.files(path = filepath, pattern = ".txt")
  
  # Create an empty data frame to load each temp result into
  tempResultsDF <- data.frame(FileNo = integer(), 
                              Band = character(), 
                              Intensity = numeric(),
                              Center = numeric(), 
                              FWHM = numeric(), 
                              Temp_D1 = numeric(),
                              Temp_D2 = numeric(),
                              Filepath = character())
  
  # Loop through each file
  for (i in 1:length(ramanResults)) {
    
    # Create the filepath for the selected file
    tempfilepath <- paste(filepath, ramanResults[i], sep = "//")

    # Read in csv file, assumes no column names
    temp <- read.table(tempfilepath, header = TRUE, nrows = 5)
    
    # Find the row number which has a peak center within the target range
    targetRow <- which(temp$Center > targetRange[1] & temp$Center < targetRange[2])
    
    if (length(targetRow) > 0) {
      tempC_D1 <- -2.15 * temp[targetRow, 5] + 478; #Kouketsu et al., 2014_Eq1
      tempC_D2 <- -6.78 * temp[targetRow, 5] + 535; #Kouketsu et al., 2014_Eq2
      
      tempResultsDF <- rbind(tempResultsDF, 
                             data.frame(FileNo = i, 
                                        Band = targetBand, 
                                        Intensity = temp[targetRow, 3],
                                        Center = temp[targetRow, 4], 
                                        FWHM = temp[targetRow, 5], 
                                        Temp_D1 = tempC_D1,
                                        Temp_D2 = tempC_D2,
                                        Filepath = tempfilepath))
    }else{
      tempResultsDF <- rbind(tempResultsDF, 
                             data.frame(FileNo = i, 
                                        Band = targetBand, 
                                        Intensity = NA,
                                        Center = NA, 
                                        FWHM = NA, 
                                        Temp_D1 = NA,
                                        Temp_D2 = NA,
                                        Filepath = tempfilepath))
    }

    
    
    
  }
  
  # Remove temperature calculated from the wrong band
  if(targetBand == "D1"){
    tempResultsDF$Temp_D2 <- NA
  }else if(targetBand == "D2"){
    tempResultsDF$Temp_D1 <- NA
  }else if(targetBand == "G"){
    tempResultsDF$Temp_D1 <- NA
    tempResultsDF$Temp_D2 <- NA
  }
  
  return(tempResultsDF)

}

#---------------------------------------------------------------------------------------------------------------------------
# Temperature Calculation #

# Calculate raman temperature using D1 (eq. 1) and D2 (eq. 2) formulas
MS_1_DF_D1 <- getRamanTemp(filepath = "MS-1//R_CSVs", targetBand = "D1")
MS_1_DF_D2 <- getRamanTemp(filepath = "MS-1//R_CSVs", targetBand = "D2")
MS_1_DF_G <- getRamanTemp(filepath = "MS-1//R_CSVs", targetBand = "G")

MS_3_DF_D1 <- getRamanTemp(filepath = "MS-3//R_CSVs", targetBand = "D1")
MS_3_DF_D2 <- getRamanTemp(filepath = "MS-3//R_CSVs", targetBand = "D2")
MS_3_DF_G <- getRamanTemp(filepath = "MS-3//R_CSVs", targetBand = "G")

MS_6_DF_D1 <- getRamanTemp(filepath = "MS-6//R_CSVs", targetBand = "D1")
MS_6_DF_D2 <- getRamanTemp(filepath = "MS-6//R_CSVs", targetBand = "D2")
MS_6_DF_G <- getRamanTemp(filepath = "MS-6//R_CSVs", targetBand = "G")


#---------------------------------------------------------------------------------------------------------------------------
# Mean & standard deviations summary table #

# Create an empty summary table following Kouketsu table 2
ramanSummary <- data.frame(Sample = character(), 
                           D1_Center = numeric(),
                           D1_Center_sd = numeric(),
                           D1_FWHM = numeric(),
                           D1_FWHM_sd = numeric(),
                           D1_Temp = numeric(),
                           D1_Temp_sd = numeric(), 
                           G_Center = numeric(),
                           G_Center_sd = numeric(),
                           G_FWHM = numeric(),
                           G_FWHM_sd = numeric(),
                           G_Temp = numeric(),
                           G_Temp_sd = numeric(), 
                           D2_Center = numeric(),
                           D2_Center_sd = numeric(),
                           D2_FWHM = numeric(),
                           D2_FWHM_sd = numeric(),
                           D2_Temp = numeric(),
                           D2_Temp_sd = numeric(),
                           D1_Intensity = numeric(),
                           D1_Intensity_sd = numeric(),
                           G_Intensity = numeric(),
                           G_Intensity_sd = numeric(),
                           D2_Intensity = numeric(),
                           D2_Intensity_sd = numeric())

# sample numbers
samples <- c(1,3,6)

# bands of interest
bands <- c("D1", "G", "D2")

# to fill summary table, loop through each band
for (b in 1:3) {
  
  bandNo <- bands[b]
  
  # loop through each sample
  for (t in 1:3) {
    
    sampleNo <- samples[t]
    
    # retrieve the table for the current sample and band
    tempTable <- get(paste("MS", sampleNo, "DF", bandNo, sep = "_"))
  
    # enter current sample into master table
    ramanSummary[t, "Sample"] <- paste("MS", sampleNo, sep = "-")
    
    # calculate all mean and sd to fill master summary table for current band
    ramanSummary[t, paste(bandNo, "Center", sep = "_")] <- mean(tempTable$Center, na.rm = T)
    ramanSummary[t, paste(bandNo, "Center_sd", sep = "_")] <- sd(tempTable$Center, na.rm = T)
  
    ramanSummary[t, paste(bandNo, "FWHM", sep = "_")] <- mean(tempTable$FWHM, na.rm = T)
    ramanSummary[t, paste(bandNo, "FWHM_sd", sep = "_")] <- sd(tempTable$FWHM, na.rm = T)
    
    ramanSummary[t, paste(bandNo, "Temp", sep = "_")] <- mean(tempTable$Temp_D1, na.rm = T)
    ramanSummary[t, paste(bandNo, "Temp_sd", sep = "_")] <- sd(tempTable$Temp_D1, na.rm = T)
    
    ramanSummary[t, paste(bandNo, "Intensity", sep = "_")] <- mean(tempTable$Intensity, na.rm = T)
    ramanSummary[t, paste(bandNo, "Intensity_sd", sep = "_")] <- sd(tempTable$Intensity, na.rm = T)
  }
  
}

# Remove G band temperatures
ramanSummary$G_Temp <- NA 
ramanSummary$G_Temp_sd <- NA

# Calculate intensity ratios
ramanSummary$D2D1Intensity <- ramanSummary$D2_Intensity/ramanSummary$D1_Intensity
ramanSummary$D1GIntensity <- ramanSummary$D1_Intensity/ramanSummary$G_Intensity

# Write summary table to csv
# write.csv(ramanSummary, "RamanSummaryResults.csv", row.names = T)

#---------------------------------------------------------------------------------------------------------------------------
# Plotting results #

# Create full data frame with all D1 results to calculate average from
allTemps_D1 <- rbind(MS_1_DF_D1, MS_3_DF_D1, MS_6_DF_D1)


# Plot box plot of temps based on D1 formula
fig_D1 <- plot_ly(data = MS_1_DF_D1, type = "box", hoverinfo = 'y') %>%
  add_trace(y = ~Temp_D1, name = "MS-1") %>%
  add_boxplot(data = MS_3_DF_D1, y = ~Temp_D1, name = "MS-3") %>%
  add_boxplot(data = MS_6_DF_D1, y = ~Temp_D1, name = "MS-6") %>%
  add_boxplot(data = allTemps_D1, y = ~Temp_D1, name = "Average") %>%
  layout(title = 'Raman Temperatures (D1)', 
         yaxis = list(title = 'Temperature (C)'))

fig_D1

# Plot histogram of temps based on D1 formula with mean (dotted line) and median (dashed line) temperatures. Grey box is 2x sd.
vline <- function(x = 0, color = "black", dash = NA) {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(color = color, dash=dash)
  )
}


# Remove outliers greater than 2x sd from dataframes
MS_1_2SDOut <-  MS_1_DF_D1[MS_1_DF_D1$Temp_D1 < (mean(MS_1_DF_D1$Temp_D1) + 2*sd(MS_1_DF_D1$Temp_D1)) &
                            MS_1_DF_D1$Temp_D1 > (mean(MS_1_DF_D1$Temp_D1) - 2*sd(MS_1_DF_D1$Temp_D1)),]

MS_3_2SDOut <-  MS_3_DF_D1[MS_3_DF_D1$Temp_D1 < (mean(MS_3_DF_D1$Temp_D1) + 2*sd(MS_3_DF_D1$Temp_D1)) &
                            MS_3_DF_D1$Temp_D1 > (mean(MS_3_DF_D1$Temp_D1) - 2*sd(MS_3_DF_D1$Temp_D1)),]

MS_6_2SDOut <-  MS_6_DF_D1[MS_6_DF_D1$Temp_D1 < (mean(MS_6_DF_D1$Temp_D1) + 2*sd(MS_6_DF_D1$Temp_D1)) &
                            MS_6_DF_D1$Temp_D1 > (mean(MS_6_DF_D1$Temp_D1) - 2*sd(MS_6_DF_D1$Temp_D1)),]

MS_allT_2SDOut <-  rbind(MS_1_2SDOut, MS_3_2SDOut, MS_6_2SDOut)


# create histograms with mean temp for each sample and an overall average
hist_MS1_2SDOut <- plot_ly(data = MS_1_DF_D1, type = "histogram", hoverinfo = 'x', x = ~Temp_D1, name = "MS-1", nbinsx = 15) %>%
  layout(shapes = list(vline(mean(MS_1_2SDOut$Temp_D1), dash = "dot"), vline(median(MS_1_DF_D1$Temp_D1), dash = "dash"),
                       list(type = "rect", fillcolor = "grey", line = list(color = "grey"), opacity = 0.5, y0 = 0, y1 = 10,
                            x0 = mean(MS_1_2SDOut$Temp_D1) - 2*sd(MS_1_2SDOut$Temp_D1), 
                            x1 = mean(MS_1_2SDOut$Temp_D1) + 2*sd(MS_1_2SDOut$Temp_D1))))

hist_MS3_2SDOut <- plot_ly(data = MS_3_DF_D1, type = "histogram", hoverinfo = 'x', x = ~Temp_D1, name = "MS-3", nbinsx = 15) %>%
  layout(shapes = list(vline(mean(MS_3_2SDOut$Temp_D1), dash = "dot"), vline(median(MS_3_DF_D1$Temp_D1), dash = "dash"),
                       list(type = "rect", fillcolor = "grey", line = list(color = "grey"), opacity = 0.5, y0 = 0, y1 = 10,
                            x0 = mean(MS_3_2SDOut$Temp_D1) - 2*sd(MS_3_2SDOut$Temp_D1), 
                            x1 = mean(MS_3_2SDOut$Temp_D1) + 2*sd(MS_3_2SDOut$Temp_D1))))

hist_MS6_2SDOut <- plot_ly(data = MS_6_DF_D1, type = "histogram", hoverinfo = 'x', x = ~Temp_D1, name = "MS-6") %>%
  layout(shapes = list(vline(mean(MS_6_2SDOut$Temp_D1), dash = "dot"), vline(median(MS_6_DF_D1$Temp_D1), dash = "dash"),
                       list(type = "rect", fillcolor = "grey", line = list(color = "grey"), opacity = 0.5, y0 = 0, y1 = 12,
                            x0 = mean(MS_6_2SDOut$Temp_D1) - 2*sd(MS_6_2SDOut$Temp_D1), 
                            x1 = mean(MS_6_2SDOut$Temp_D1) + 2*sd(MS_6_2SDOut$Temp_D1))))

hist_Avg_2SDOut <- plot_ly(data = allTemps_D1, type = "histogram", hoverinfo = 'x', x = ~Temp_D1, name = "Avg") %>%
  layout(shapes = list(vline(mean(MS_allT_2SDOut$Temp_D1), dash = "dot"), vline(median(allTemps_D1$Temp_D1), dash = "dash"),
                       list(type = "rect", fillcolor = "grey", line = list(color = "grey"), opacity = 0.5, y0 = 0, y1 = 30,
                            x0 = mean(MS_allT_2SDOut$Temp_D1) - 2*sd(MS_allT_2SDOut$Temp_D1), 
                            x1 = mean(MS_allT_2SDOut$Temp_D1) + 2*sd(MS_allT_2SDOut$Temp_D1))))

# subplot all histograms
hist_D1_2SDOut <- subplot(hist_MS1_2SDOut, hist_MS3_2SDOut, hist_MS6_2SDOut, hist_Avg_2SDOut, nrows = 2, titleX = T, shareX = T) %>%
  layout(title = list(text = 'Raman Temperatures (D1) - 2SD outliers'), xaxis = list(title = 'Temperature (C)'),
         annotations = list(x = mean(MS_allT_2SDOut$Temp_D1), y = 30, xref = "x2", yref = "y4",
                            bgcolor = "white", text = paste("mean =", round(mean(MS_allT_2SDOut$Temp_D1),1))))

hist_D1_2SDOut


## Script containing functions used to calculate stress from calcite twin density following Ryback et al., 2013 paleopiezometer
## Created by Madison Frank, updated 2023/05/08

# calculate differential stress from calcite twin density using Rybacki et al., 2013 equation
getRawTwins <- function(filepath, sheet, startRow, colnames){
  library(openxlsx)
  
  # Read in excel file
  ca_twin_raw <- read.xlsx(filepath, sheet = sheet, startRow = startRow, colNames = T)

  
  #---------------------------------------------------------------------------------------------
  ### Thin section data and conversion to differential stress
  
  # Create a data frame holding the number of twins measured over a given distance. Only keep rows with data (complete cases)
  twin_data <- ca_twin_raw[,colnames]
  twin_data <- twin_data[complete.cases(twin_data),]
  
  colnames(twin_data) <- c("name", "localisation", "distance_um", "twins")
  str(twin_data)
  
  # Convert um to mm
  twin_data$distance_mm <- twin_data$distance_um/1000
  
  # Calculate twins/mm
  twin_data$density <- twin_data$twins/twin_data$distance_mm
  
  # Calculate differential stress following Rybacki et al., 2013
  twin_data$diffStress <- 19.5 * sqrt(twin_data$density)
  
  twin_data$error_high <- (19.5 + 9.8) * sqrt(twin_data$density)
  twin_data$error_low <- (19.5 - 9.8) * sqrt(twin_data$density)
  
  return(twin_data)
}

# Determine outliers using IQR method
removeTwinOutliers <- function(dataset){
  
  no_outliers <- data.frame()
  
  for (i in 1:length(unique(dataset$name))) {
    
    temp_data <- dataset[dataset$name == unique(dataset$name)[i],]
    
    # define quantiles defining outliers
    Q1 <- quantile(temp_data$density, .25)
    Q3 <- quantile(temp_data$density, .75)
    IQR <- IQR(temp_data$density)
    
    # only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
    no_outliers <- rbind(no_outliers,
                         subset(temp_data, temp_data$density> (Q1 - 1.5*IQR) & temp_data$density< (Q3 + 1.5*IQR)))
    
  }
  
  return(no_outliers)
}


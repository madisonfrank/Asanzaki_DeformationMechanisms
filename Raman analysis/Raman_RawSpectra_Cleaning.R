### SCRIPT TO REMOVE NAN LINES FROM RAMAN ANALYSIS CSV FILES ###
## Created by Madison Frank, updated 2023/06/08


# Select folder with RAMAN data
filepath <- "Raman//MS-1"

# Get list of csv file names
ramanFiles <- list.files(path = filepath, pattern = ".CSV")

# Loop through each file
for (i in 1:length(ramanFiles)) {
  
  # Create the filepath for the selected file
  tempfilepath <- paste(filepath, ramanFiles[i], sep = "//")
  
  # Read in csv file, assumes no column names
  temp <- read.csv(tempfilepath, header = FALSE)
  
  # Find any NaN lines
  missingData <- grepl(pattern = "#NaN", x = temp[,2])
  
  # Remove NaN lines
  temp <- temp[!missingData,]

  # Create the export filepath for a new file (R_old file name) - ASSUMES A FOLDER CALLED R_CSVs IN THE SELECTED FOLDER!!
  r_filepath <- paste(filepath, "R_CSVs", paste("R", ramanFiles[i], sep = "_"), sep = "//")
  
  # Write new csv file - use comma separator, no row or column names, no commas around values
  write.table(x = temp, file = r_filepath, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}




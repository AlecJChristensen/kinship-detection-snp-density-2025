#Colony.format.file 600K QC 15.20
#Wild Cervid Project
#Obj.One 
#Last Updated 7.26.24

# Load the data.table package
library(data.table)
library(tibble)
library(dplyr)
# Specify the file path
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/data_output/600.15.20_numeric_7.11.24"

# Read the file into a data.table
genotype_data <- fread(file_path)

###################################################################################################################
#colony
#DATA WRANGLIN!

#Bring in metadata
#to merge together with VCF
# in order to pull different sets of individuals 

#rename to make collective
gdf<-genotype_data

# Convert "ind_id" to row names
gdf <- gdf %>% column_to_rownames(var = "ind_id")


#taking the individuals ID's and putting in new data frame
row_names <- rownames(gdf)
# Create a new data frame with the row names
row_names_df <- data.frame(row_names)
#rename the column in the new data frame
colnames(row_names_df) <- "WHP_sample_id"

#loading in Emily Spreadsheet of individuals
metadata <- read.csv("C:/Users/alec1/Desktop/SNP_data/MNDNR_SamplesForSNPsProject_10182023_230pm (8).csv", na.strings=c("","NA"))

#merging the data frame followed by adding in duplicated individuals
new.merged <- merge(row_names_df, metadata, by = "WHP_sample_id", sort = FALSE)


#merging the data frame followed by adding in duplicated individuals
new.merged <- merge(row_names_df, metadata, by = "WHP_sample_id", sort = FALSE)

new_row <- data.frame(WHP_sample_id = "DUP_MN184503",
                      year = "2023",
                      sample_category = "Targeted",
                      sample_acquisition = "Agency Culled",
                      date_harvested = "1/23/2023",
                      date_collected = "1/24/2023",
                      permit_area = "648",
                      town = "103",
                      range = "10",
                      section = "33",
                      sex = "Male",
                      age = "Yearling",
                      Genetic.Tissue.Available = "LN homogenate",
                      Sample.location..bag. = "CWD+ (1 bag)",
                      Known_Age = "1.5",
                      CWD_status = "Positive",
                      lab_notes = "OD value 0.64",
                      original_data_resolution = "UTM",
                      UTME = "577220",
                      UTMN = "483660",
                      microsat.work.done. = NA,
                      Obj1Grouping = "GroupA",
                      Obj1 = "X",
                      Obj2 = NA,
                      Obj3 = "X",
                      Obj4 = "X",
                      Number.of.Objectives = "3")
new.merged <- merge(row_names_df, metadata, by = "WHP_sample_id", sort = FALSE)

# Append the new row to merged_df
merged.df <- rbind(new.merged, new_row)

###############################
####################
##########
######

#merge genetic data
gdf <- rownames_to_column(gdf, var = "WHP_sample_id")



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# offspring

#taking offspring
offspring <- subset(merged.df, age %in% c("Fawn", "Fetus", "Yearling"))
offspringID <- subset(offspring, select = WHP_sample_id)

# Merge the data frames based on the common column "WHP_sample_id"
OffspringC <- merge(offspringID, gdf, by = "WHP_sample_id", all.x = TRUE)
# 'all.x = TRUE' keeps all rows from the left data frame (offspringID)
# If you want to keep all rows from the right data frame (gdf), use 'all.y = TRUE' instead

##########################
##########################
###################################################################################################################

# male adults

#number of adult male 
adults <- subset(merged.df, age == "Adult")
maleadults <- subset(adults, sex == "Male")
maleadultsID <- subset(maleadults, select = WHP_sample_id)
maleadults <- merge(maleadultsID, gdf, by = "WHP_sample_id", all.x = TRUE)


###################################################################################################################
##########################
##########################

#female adults

adults <- subset(merged.df, age == "Adult")
femaleadults <- subset(adults, sex == "Female")
femaleadultsID <- subset(femaleadults, select = WHP_sample_id)
femaleadults <- merge(femaleadultsID, gdf, by = "WHP_sample_id", all.x = TRUE)


###################################################################################################################
##########################
##########################
###################################################################################################################

#little switcharoo
#MN164234 is a fawn with offspring. With this, needs to be switched groups
# Find the row index of the individual MN164234 in offspringc
index <- which(OffspringC$WHP_sample_id == "MN164234")

# If the individual exists in offspringc
if (length(index) > 0) {
  # Extract the row corresponding to MN164234
  individual_info <- OffspringC[index, ]
  
  # Add the individual to femaleadults
  femaleadults <- rbind(femaleadults, individual_info)
  
  # Remove the individual from offspringc
  OffspringC <- OffspringC[-index, ]
}

num_rows <- nrow(femaleadults)
femaleadults <- femaleadults[1:(num_rows-1), ]


save.image(file = "my_R_colony600K.1520.environment.RData")

###################################################################################################################
###################################################################################################################
# Define the directory path
directory_path <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\COLONY\\Input_600KQC15.20"

# Define the full file path with the desired file name
file_path <- file.path(directory_path, "femaleadult600k15.20.txt")

# Export the data frame to a text file
write.table(femaleadults, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#
# Define the directory path
directory_path <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\COLONY\\Input_600KQC15.20"

# Define the full file path with the desired file name
file_path <- file.path(directory_path, "maleadult600k15.20.txt")
# Export the data frame to a text file
write.table(maleadults, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



# Define the directory path
directory_path <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\COLONY\\Input_600KQC15.20"

# Define the full file path with the desired file name
file_path <- file.path(directory_path, "offspring600k15.20.txt")

# Export the data frame to a text file
write.table(OffspringC, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




# Drop the "WHP_sample_id" column and create a new data frame "error.rate"
error.rate <- maleadults[, !names(maleadults) %in% "WHP_sample_id"]

# Drop all rows from "error.rate"
error.rate <- error.rate[0, ]

# Drop columns whose names end with "_1"
error.rate <- error.rate[, !grepl("_1$", names(error.rate))]

# Add three rows of "0" for each column
error.rate <- rbind(error.rate, matrix(0, nrow = 3, ncol = ncol(error.rate)))

# Replace the last row with "0.175" for all columns
error.rate[nrow(error.rate), ] <- 0.175

# Convert matrix to data frame and set column names back
error.rate <- as.data.frame(error.rate)
colnames(error.rate) <- names(maleadults[, !names(maleadults) %in% c("WHP_sample_id", grep("_1$", names(maleadults), value = TRUE))])



# Define the directory path
directory_path <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\COLONY\\Input_600KQC15.20"

# Define the full file path with the desired file name
file_path <- file.path(directory_path, "error.rate0175.600k15.20.txt")

# Export the data frame to a text file
write.table(error.rate, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
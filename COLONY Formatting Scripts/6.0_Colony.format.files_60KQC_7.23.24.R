#Colony format for 60KQC
#Wild Cervid Project
#Obj.One 
#Last Updated 7.23.24

# Making colony files
#4.23.24
###################################################################################################################
##########################
##########################

#loading in packages
library(tidyverse)
library(vcfR)
library(adegenet)

vcf_path <- "C:\\Users\\alec1\\Desktop\\SNP_data\\data\\LowThresholds_SNPQC\\167001_LowThresholds_LowSNPQCs\\1670010002_LowThreshold_LowSNPQCVCF.vcf"
# Read the VCF file
vcf60 <- read.vcfR(vcf_path)
#make into genind
OvSNP60_genind <- vcfR2genind(vcf60)
#turn into Data frame
df60<-genind2df(OvSNP60_genind)

##########################
##########################
###################################################################################################################
#DATA WRANGLIN!

#Bring in metadata
#to merge together with VCF
# in order to pull different sets of individuals 

#loading in Emily Spreadsheet of individuals
metadata <- read.csv("C:/Users/alec1/Desktop/SNP_data/MNDNR_SamplesForSNPsProject_10182023_230pm (8).csv", na.strings=c("","NA"))
#renaming the data frame that was switched from genind
gdf<-df60
#taking the individuals ID's and putting in new data frame
row_names <- rownames(gdf)
# Create a new data frame with the row names
row_names_df <- data.frame(row_names)
#rename the column in the new data frame
colnames(row_names_df) <- "WHP_sample_id"


# Extract and edit the names in the "IID" column
row_names_df$WHP_sample_id <- str_replace(row_names_df$WHP_sample_id, "^[^_]*_(.*?)\\.CEL$", "\\1")
# Remove everything before the second underscore
row_names_df$WHP_sample_id <- str_extract(row_names_df$WHP_sample_id, "(?<=_).*")

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

###

#merge genetic data
gdf <- rownames_to_column(gdf, var = "WHP_sample_id")

# Extract and edit the names in the "IID" column
gdf$WHP_sample_id <- str_replace(gdf$WHP_sample_id, "^[^_]*_(.*?)\\.CEL$", "\\1")
# Remove everything before the second underscore
gdf$WHP_sample_id <- str_extract(gdf$WHP_sample_id, "(?<=_).*")


###################################################################################################################
##########################
##########################

#loci error rate
# Create a new data frame with column names of df60
new_df <- data.frame(matrix(0, ncol = ncol(df60), nrow = 3))
colnames(new_df) <- colnames(df60)

# Set the last row to 0.175
new_df[nrow(new_df), ] <- 0.175

file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.error.SNP.input.175"
## Export the data frame to a text file
#save as tab delim
write.table(new_df, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
##write as table
#new_dfsubset <- new_df[, 1:5]
#file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.error.SNP"
#write.table(new_dfsubset, file = file_path, sep = "\t", row.names = FALSE, col.names = TRUE,quote =FALSE)
##46861

##########################
##########################
###################################################################################################################
# offspring

#taking offspring
offspring <- subset(merged.df, age %in% c("Fawn", "Fetus", "Yearling"))
offspringID <- subset(offspring, select = WHP_sample_id)

# Merge the data frames based on the common column "WHP_sample_id"
OffspringC <- merge(offspringID, gdf, by = "WHP_sample_id", all.x = TRUE)
# 'all.x = TRUE' keeps all rows from the left data frame (offspringID)
# If you want to keep all rows from the right data frame (gdf), use 'all.y = TRUE' instead

####

#making into COLONY Friendly
#replace 1 with 2...
OffspringC[,c(2:46862)] <- data.frame(lapply(OffspringC[,c(2:46862)], function(x) {
  gsub("1", "2", x)}))
# ...and 0 with 1
OffspringC[,c(2:46862)] <- data.frame(lapply(OffspringC[,c(2:46862)], function(x) {
  gsub("0", "1", x)}))  

#replace NA with 00 (still need to split it)
OffspringC[is.na(OffspringC)] <- "00"

#separate each column
OffspringC <- OffspringC %>% gather(snp, value, -WHP_sample_id) %>% separate(col = value, into = c("SNP1.1","SNP1.2"), sep=1)

#get back into wide
OffspringC <- OffspringC |> pivot_wider(names_from = snp, values_from=c("SNP1.1", "SNP1.2"))

#get the neighbors back together
#get just the non"SNP1.1" and "SNP1.2" parts of the name
col_no_snp <- names(OffspringC)
ending_chars <- substr(col_no_snp, 8, nchar(col_no_snp))

#order the column names based on the characters after removing the first 7 characters
sorted_col_names <- col_no_snp[order(ending_chars)]

#reorder the data frame based on the sorted column names
OffspringC <- OffspringC[, sorted_col_names]

#move sample ID to first column
#get the last column
last_col <- OffspringC[, ncol(OffspringC)]

#remove the last column from the dataframe
OffspringC <- OffspringC[, -ncol(OffspringC)]

#insert the last column as the first column
OffspringC <- cbind(last_col, OffspringC)


file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.offspring.input"
## Export the data frame to a text file
#save as tab delim
write.table(OffspringC, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

#subset
#subsetoffspring <- OffspringC[, 1:9]
#file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.offspring.subset"
#write.table(subsetoffspring , file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
##########################
##########################
###################################################################################################################

# male adults

#number of adult male 
adults <- subset(merged.df, age == "Adult")
maleadults <- subset(adults, sex == "Male")
maleadultsID <- subset(maleadults, select = WHP_sample_id)
maleadults <- merge(maleadultsID, gdf, by = "WHP_sample_id", all.x = TRUE)

maleadults[,c(2:46862)] <- data.frame(lapply(maleadults[,c(2:46862)], function(x) {
  gsub("1", "2", x)}))
# ...and 0 with 1
maleadults[,c(2:46862)] <- data.frame(lapply(maleadults[,c(2:46862)], function(x) {
  gsub("0", "1", x)}))  

#replace NA with 00 (still need to split it)
maleadults[is.na(maleadults)] <- "00"

#separate each column
maleadults <- maleadults %>% gather(snp, value, -WHP_sample_id) %>% separate(col = value, into = c("SNP1.1","SNP1.2"), sep=1)

#get back into wide
maleadults <- maleadults|> pivot_wider(names_from = snp, values_from=c("SNP1.1", "SNP1.2"))

#get the neighbors back together
#get just the non"SNP1.1" and "SNP1.2" parts of the name
col_no_snp <- names(maleadults)
ending_chars <- substr(col_no_snp, 8, nchar(col_no_snp))

#order the column names based on the characters after removing the first 7 characters
sorted_col_names <- col_no_snp[order(ending_chars)]

#reorder the data frame based on the sorted column names
maleadults <- maleadults[, sorted_col_names]

#move sample ID to first column
#get the last column
last_col <- maleadults[, ncol(maleadults)]

#remove the last column from the dataframe
maleadults <- maleadults[, -ncol(maleadults)]

#insert the last column as the first column
maleadults <- cbind(last_col, maleadults)

# Export the data frame to a text file
file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/male.adults.input"
write.table(maleadults, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

#subsetmaleadults <- maleadults[, 1:9]
#file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.offspring.subsetmaleadults"
#write.table(subsetmaleadults , file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
###################################################################################################################
##########################
##########################

#female adults

adults <- subset(merged.df, age == "Adult")
femaleadults <- subset(adults, sex == "Female")
femaleadultsID <- subset(femaleadults, select = WHP_sample_id)
femaleadults <- merge(femaleadultsID, gdf, by = "WHP_sample_id", all.x = TRUE)


femaleadults[,c(2:46862)] <- data.frame(lapply(femaleadults[,c(2:46862)], function(x) {
  gsub("1", "2", x)}))
# ...and 0 with 1
femaleadults[,c(2:46862)] <- data.frame(lapply(femaleadults[,c(2:46862)], function(x) {
  gsub("0", "1", x)}))  

#replace NA with 00 (still need to split it)
femaleadults[is.na(femaleadults)] <- "00"

#separate each column
femaleadults <- femaleadults %>% gather(snp, value, -WHP_sample_id) %>% separate(col = value, into = c("SNP1.1","SNP1.2"), sep=1)

#get back into wide
femaleadults <- femaleadults|> pivot_wider(names_from = snp, values_from=c("SNP1.1", "SNP1.2"))

#get the neighbors back together
#get just the non"SNP1.1" and "SNP1.2" parts of the name
col_no_snp <- names(femaleadults)
ending_chars <- substr(col_no_snp, 8, nchar(col_no_snp))

#order the column names based on the characters after removing the first 7 characters
sorted_col_names <- col_no_snp[order(ending_chars)]

#reorder the data frame based on the sorted column names
femaleadults <- femaleadults[, sorted_col_names]

#move sample ID to first column
#get the last column
last_col <- femaleadults[, ncol(femaleadults)]

#remove the last column from the dataframe
femaleadults <- femaleadults[, -ncol(femaleadults)]

#insert the last column as the first column
femaleadults <- cbind(last_col, femaleadults)

file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/female.adults.input"
# Export the data frame to a text file
write.table(femaleadults, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

#subsetfemaleadults <- femaleadults[, 1:9]
#file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.offspring.subsetfemaleadults"
#write.table(subsetfemaleadults , file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###################################################################################################################
##########################
##########################
###################################################################################################################

#known Paternal sibs = NA/unknown. 


###################################################################################################################
##########################
##########################
###################################################################################################################

#known Maternal sibs (done off of R)

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

file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/female.adults.input.plusMN164234"
# Export the data frame to a text file
write.table(femaleadults, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

file_path <- "C:/Users/alec1/Desktop/SNP_data/output/COLONY_SNP_60K/colony.offspring.input.minusMN164234"
## Export the data frame to a text file
#save as tab delim
write.table(OffspringC, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
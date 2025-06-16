#Coancestry input file
#Wild Cervid Project
#Obj.One 
#Last Updated 7.23.24

library(vcfR)
library(adegenet)
library(stringr)
library(tibble)
library(tidyr)
library(data.table)
library(tibble)
library(dplyr)
vcf_path60 <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/Random600.pilot.vcf"
vcf60 <- read.vcfR(vcf_path60)

OvSNP60_genind <- vcfR2genind(vcf60)
genotype60K <-genind2df(OvSNP60_genind)

# Specify the path where you want to save the CSV file
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/gdf.csv"
# Export the DataFrame as a CSV file

###################################################################################################################


###################################################################################################################
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
gdf<-genotype60K
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

################################################################################################################
gdf[,c(2:601)] <- data.frame(lapply(gdf[,c(2:601)], function(x) {
  gsub("1", "2", x)}))
# ...and 0 with 1
gdf[,c(2:601)] <- data.frame(lapply(gdf[,c(2:601)], function(x) {
  gsub("0", "1", x)}))  

#replace NA with 00 (still need to split it)
gdf[is.na(gdf)] <- "00"

#separate each column
gdf <- gdf %>% gather(snp, value, -WHP_sample_id) %>% separate(col = value, into = c("SNP1.1","SNP1.2"), sep=1)

#get back into wide
gdf <- gdf |> pivot_wider(names_from = snp, values_from=c("SNP1.1", "SNP1.2"))

#get the neighbors back together
#get just the non"SNP1.1" and "SNP1.2" parts of the name
col_no_snp <- names(gdf)
ending_chars <- substr(col_no_snp, 8, nchar(col_no_snp))

#order the column names based on the characters after removing the first 7 characters
sorted_col_names <- col_no_snp[order(ending_chars)]

#reorder the data frame based on the sorted column names
gdf <- gdf[, sorted_col_names]

#move sample ID to first column
#get the last column
last_col <- gdf[, ncol(gdf)]

#remove the last column from the dataframe
gdf <- gdf[, -ncol(gdf)]

#insert the last column as the first column
gdf <- cbind(last_col, gdf)


# Define the file path where you want to save the R object
save_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY/600locicoancestry.RData"

# Save the dataframe as an R object
save(gdf, file = save_path)
###################################################################################################################
###################################################################################################################

write.table(gdf, file = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY/600loci.coancestryinput.file", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save the entire environment
#save.image(file = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY/600K.15.20_numeric.dataframe_7.2.24_environment.RData")

# Save only the combined_numeric_df dataframe
#save(combined_numeric_df, file = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY600K.15.20_numeric.dataframe_7.2.24_dataframe.RData")

# Save combined_numeric_df as an RData file
#save(combined_numeric_df, file = "combined_numeric_df.RData")

#save.image(file = "my_R_environment.RData")

# Alternatively, save as an .rds file
#saveRDS(combined_numeric_df, file = "combined_numeric_df.rds")
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
vcf_path60 <- "Objective.One_60Kvs.600K_Relatedness\\Data_raw\\LOWThresholdsSNPQCs\\1670010002_LowThreshold_LowSNPQCVCF.vcf"
vcf60 <- read.vcfR(vcf_path60)

OvSNP60_genind <- vcfR2genind(vcf60)
genotype60K <-genind2df(OvSNP60_genind)

# Specify the path where you want to save the CSV file
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Coancestry/Data/gdf.csv"
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
gdf[,c(2:46862)] <- data.frame(lapply(gdf[,c(2:46862)], function(x) {
  gsub("1", "2", x)}))
# ...and 0 with 1
gdf[,c(2:46862)] <- data.frame(lapply(gdf[,c(2:46862)], function(x) {
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
save_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Saved_R_Objects/Coancestry/60Kcoancestry.RData"

# Save the dataframe as an R object
save(gdf, file = save_path)
###################################################################################################################
###################################################################################################################

#600KQC
library(data.table)

# Specify the paths to your PLINK files
bim_file <- "600QC_filtered15.20.bim"
fam_file <- "600QC_filtered15.20.fam"
ped_file <- "600QC_filtered15.20.ped"

# Read the .bim file
bim_df <- fread(bim_file, header = FALSE)
setnames(bim_df, c("chr", "snp_id", "gen_dist", "bp_pos", "allele1", "allele2"))

# Read the .fam file
fam_df <- fread(fam_file, header = FALSE)
setnames(fam_df, c("fam_id", "ind_id", "pat_id", "mat_id", "sex", "phenotype"))

# Read the .ped file
ped_df <- fread(ped_file, header = FALSE)

# Extract the genotype data from the .ped file (starting from the 7th column)
genotype_df <- ped_df[, 7:ncol(ped_df), with = FALSE]

# Combine the .fam data with the genotype data
combined_df <- cbind(fam_df, genotype_df)

# Generate column names for the genotype data
snp_ids <- bim_df$snp_id
genotype_colnames <- c(rbind(snp_ids, paste0(snp_ids, "_1")))

# Assign column names to the combined data frame
setnames(combined_df, old = names(combined_df)[7:ncol(combined_df)], new = genotype_colnames)

# Create a function to map alleles to numeric values based on the reference alleles
map_genotypes_to_numeric <- function(genotype_data, bim_data) {
  numeric_genotypes <- as.data.table(matrix(nrow = nrow(genotype_data), ncol = ncol(genotype_data)))
  
  for (i in seq(1, nrow(bim_data))) {
    ref_allele1 <- bim_data$allele1[i]
    ref_allele2 <- bim_data$allele2[i]
    
    allele1_col <- genotype_data[[2 * (i - 1) + 1]]
    allele2_col <- genotype_data[[2 * (i - 1) + 2]]
    
    numeric_genotypes[[2 * (i - 1) + 1]] <- sapply(allele1_col, function(a) {
      if (is.na(a)) {
        return(NA)
      } else if (a == ref_allele1) {
        return(0)
      } else if (a == ref_allele2) {
        return(1)
      } else {
        return(NA)
      }
    })
    
    numeric_genotypes[[2 * (i - 1) + 2]] <- sapply(allele2_col, function(a) {
      if (is.na(a)) {
        return(NA)
      } else if (a == ref_allele1) {
        return(0)
      } else if (a == ref_allele2) {
        return(1)
      } else {
        return(NA)
      }
    })
  }
  
  setnames(numeric_genotypes, old = names(numeric_genotypes), new = genotype_colnames)
  return(numeric_genotypes)
}

# Apply the mapping function to the genotype data
genotype_numeric <- map_genotypes_to_numeric(combined_df[, 7:ncol(combined_df), with = FALSE], bim_df)

# Combine the fam_df with the numeric genotype data
combined_numeric_df <- cbind(fam_df, genotype_numeric)

# Check the dimensions to ensure no columns are dropped
print(dim(combined_df))
print(dim(combined_numeric_df))

# Save the entire environment
save.image(file = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Saved_R_Objects/600K.15.20_numeric.dataframe_7.2.24_environment.RData")

# Save only the combined_numeric_df dataframe
save(combined_numeric_df, file = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Saved_R_Objects/600K.15.20_numeric.dataframe_7.2.24_dataframe.RData")

# Save combined_numeric_df as an RData file
save(combined_numeric_df, file = "combined_numeric_df.RData")

save.image(file = "my_R_environment.RData")

# Alternatively, save as an .rds file
saveRDS(combined_numeric_df, file = "combined_numeric_df.rds")
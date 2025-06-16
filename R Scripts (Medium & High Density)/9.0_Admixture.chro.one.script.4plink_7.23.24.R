#Admixture: text file script to turn all chromosomes to one
#Wild Cervid Project
#Obj.One 
#Last Updated 7.23.24

library(vcfR)
library(adegenet)
#trying to alter 60K chromosomes than turn back into a vcf
#loading in data 
# Specify the path to your VCF file
vcf_path <- "Objective.One_60Kvs.600K_Relatedness\\Data_raw\\LOWThresholdsSNPQCs\\1670010002_LowThreshold_LowSNPQCVCF.vcf"
vcf_data <- read.vcfR(vcf_path)# Read the VCF file

deergenind <- vcfR2genind(vcf_data)
deerdf <- genind2df(deergenind)

######################################
######################################
######################################

# Step 1: Extract the column names
column_names <- colnames(deerdf)

# Step 2: Create a new data frame with the column names in the first column
# and "1" in the second column
chromosome.is.one <- data.frame(
  Column_Name = column_names,
  Value = 1
)

# Define the file path
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/PLINK/Files/Chromosome.is.one.txt"

# Export the data frame to a text file
write.table(chromosome.is.one, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE)



###################################################################################################################
###################################################################################################################
###################################################################################################################

#600K
load("C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Saved_R_Objects/com.df.RData")

# Step 1: Extract the column names
column_names600 <- colnames(com.df)

# Step 2: Create a new data frame with the column names in the first column
# and "1" in the second column
chromosome.is.one600 <- data.frame(
  Column_Name600 = column_names600,
  Value = 1
)
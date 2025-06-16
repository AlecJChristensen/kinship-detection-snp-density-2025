library(readxl)
library(dplyr)

# Define the file path
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/MissQC.xlsx"

# Read the Excel file into a dataframe
missing <- read_excel(file_path, sheet = 2)

# Check the first few rows
head(missing)


# Create a new dataframe with rows where F_MISS == 0
missing.of.0 <- missing %>% filter(F_MISS == 0)

Random600.pilot <- missing.of.0 %>% sample_n(600)
head(Random600.pilot)


# Create a new dataframe with the "ID" column and add a "chromosome" column
Random600.pilot.chromosome <- Random600.pilot %>%
  select(ID) %>%
  mutate(chromosome = 1)

# View the new dataframe
head(Random600.pilot.chromosome)

# Reorder columns so that "ID" comes first and "chromosome" comes second
Random600.pilot.chromosome <- Random600.pilot.chromosome %>%
  select(ID, chromosome)

Random600.pilot.chromosome1 <- data.frame(
  chromosome = Random600.pilot.chromosome$chromosome,
  ID = Random600.pilot.chromosome$ID
)

Random600.pilot.chromosome0 <- Random600.pilot.chromosome1 %>%
  mutate(chromosome = 0)

Random600.pilot.chromosome00 <- Random600.pilot.chromosome %>%
  mutate(chromosome = 0)

Random600.pilotloci <- Random600.pilot.chromosome00[, "ID", drop = FALSE]

# Define the file path
file_path <- "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/Random600.pilotloci.txt"

# Write the dataframe to a text file without headers
write.table(Random600.pilotloci, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

###################################################################################################################
###################################################################################################################
###################################################################################################################

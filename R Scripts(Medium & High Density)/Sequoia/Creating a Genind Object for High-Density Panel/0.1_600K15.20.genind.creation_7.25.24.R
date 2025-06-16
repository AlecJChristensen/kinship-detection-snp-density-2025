#creating 600K 15.20 genind
#Wild Cervid Project
#Obj.One 
#Last Updated 7.24.24

#packages
library(adegenet)
library(data.table)
library(dplyr)

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

# Drop the specified columns
combined_numeric_df <- combined_numeric_df %>%
  select(-fam_id, -pat_id, -mat_id, -sex, -phenotype)

################################################################################
# Function to combine allele columns
combine_alleles <- function(df) {
  # Extract individual IDs
  individual_ID <- df[[1]]
  
  # Get column names
  column_names <- colnames(df)
  
  # Initialize an empty list to store new columns
  new_columns <- list()
  
  # Loop through the columns in pairs
  for (i in seq(2, ncol(df), by = 2)) {
    # Combine alleles
    combined_allele <- paste0(df[[i]], df[[i + 1]])
    
    # Get the base column name (remove "_1" if it exists)
    base_colname <- sub("_1$", "", column_names[i])
    
    # Add combined allele column to the list
    new_columns[[base_colname]] <- combined_allele
  }
  
  # Create a new data frame with the combined columns
  combined_df <- cbind(individual_ID, as.data.frame(new_columns))
  
  # Set the column names correctly
  colnames(combined_df)[1] <- "individual_ID"
  
  return(combined_df)
}

# Apply the function data frame
com.df <- combine_alleles(combined_numeric_df)


################################################################################

fix_colnames <- function(names) {
  sapply(names, function(name) {
    parts <- unlist(strsplit(name, "\\.")) # Split by period
    if (length(parts) > 2) {
      # If more than two parts, concatenate first two with a single period
      paste(parts[1], parts[2], sep = ".")
    } else {
      name # If exactly two parts, keep as is
    }
  })
}

# Apply the function to rename columns
new_colnames <- fix_colnames(colnames(com.df))
colnames(com.df) <- new_colnames

# Clean up and ensure no extra spaces or characters
clean_colnames <- function(names) {
  names <- gsub("\\s+", "", names) # Remove any whitespace
  names <- gsub("_", ".", names) # Replace underscores with periods
  names <- gsub("\\.+", ".", names) # Replace multiple periods with a single one
  sapply(names, function(name) {
    parts <- unlist(strsplit(name, "\\.")) # Split by period
    if (length(parts) > 2) {
      paste(parts[1], parts[2], sep = ".") # Concatenate first two parts
    } else if (length(parts) == 2) {
      paste(parts[1], parts[2], sep = ".") # Ensure the correct format
    } else {
      paste(name, "allele", sep = ".") # Append 'allele' if only one part
    }
  })
}

# Apply the function to rename columns
colnames(com.df) <- clean_colnames(colnames(com.df))

# Replace periods with dashes in column names
colnames(com.df) <- gsub("\\.", "-", colnames(com.df))

#create genind
genind.600K <- df2genind(com.df, sep = "\t")



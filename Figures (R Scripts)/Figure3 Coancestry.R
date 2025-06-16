# Load necessary libraries
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
# Read in the data
relatedness_estimates_600loci <- read_excel("C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY/600loci_relatedness.estimates.xlsx")
known_parentage_siblings_600loci <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")

# Process 60K dataset
known_parentage_siblings_600loci_flipped <- known_parentage_siblings_600loci
known_parentage_siblings_600loci_flipped <- known_parentage_siblings_600loci[, c("Individual2", "Individual1", setdiff(names(known_parentage_siblings_600loci), c("Individual1", "Individual2")))]
names(known_parentage_siblings_600loci_flipped)[1:2] <- c("Individual1", "Individual2")

estimates600loci_original <- merge(known_parentage_siblings_600loci, relatedness_estimates_600loci, by = c("Individual1", "Individual2"), all.x = TRUE)
estimates600loci_flipped <- merge(known_parentage_siblings_600loci_flipped, relatedness_estimates_600loci, by = c("Individual1", "Individual2"), all.x = TRUE)

estimates_600.loci <- rbind(estimates600loci_original, estimates600loci_flipped)
estimates_600.loci <- unique(estimates_600.loci)
estimates_600.loci <- na.omit(estimates_600.loci)

###################################################################################################################
###################################################################################################################
###################################################################################################################


# Read in the data
relatedness_estimates_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_relatedness.estimates.csv")
known_parentage_siblings_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")
relatedness_estimates_600K <- read_csv("C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Coancestry\\Output\\600QC15.20_relatedness.estimate.results_7.15.24.csv")
known_parentage_siblings_600K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")


# Process 60K dataset
known_parentage_siblings_60K_flipped <- known_parentage_siblings_60K
known_parentage_siblings_60K_flipped <- known_parentage_siblings_60K_flipped[, c("Individual2", "Individual1", setdiff(names(known_parentage_siblings_60K), c("Individual1", "Individual2")))]
names(known_parentage_siblings_60K_flipped)[1:2] <- c("Individual1", "Individual2")

estimates60_original <- merge(known_parentage_siblings_60K, relatedness_estimates_60K, by = c("Individual1", "Individual2"), all.x = TRUE)
estimates60_flipped <- merge(known_parentage_siblings_60K_flipped, relatedness_estimates_60K, by = c("Individual1", "Individual2"), all.x = TRUE)

estimates_60K <- rbind(estimates60_original, estimates60_flipped)
estimates_60K <- unique(estimates_60K)
estimates_60K <- na.omit(estimates_60K)

# Process 600K dataset
known_parentage_siblings_600K_flipped <- known_parentage_siblings_600K
known_parentage_siblings_600K_flipped <- known_parentage_siblings_600K_flipped[, c("Individual2", "Individual1", setdiff(names(known_parentage_siblings_600K), c("Individual1", "Individual2")))]
names(known_parentage_siblings_600K_flipped)[1:2] <- c("Individual1", "Individual2")

estimates600K_original <- merge(known_parentage_siblings_600K, relatedness_estimates_600K, by = c("Individual1", "Individual2"), all.x = TRUE)
estimates600K_flipped <- merge(known_parentage_siblings_600K_flipped, relatedness_estimates_600K, by = c("Individual1", "Individual2"), all.x = TRUE)

estimates_600K <- rbind(estimates600K_original, estimates600K_flipped)
estimates_600K <- unique(estimates_600K)
estimates_600K <- na.omit(estimates_600K)

# Combine 60K, GTseq, and 600K datasets
combine_data <- function(data, dataset_name) {
  data %>%
    pivot_longer(cols = c("TrioML", "Wang", "LynchLi", "LynchRd", "Ritland", "QuellerGt", "DyadML"), 
                 names_to = "Estimator_Type", 
                 values_to = "Genetic_Relatedness_Estimate") %>%
    mutate(Dataset = dataset_name)
}

snp60_long <- combine_data(estimates_60K, "60K")
GTseq_long <- combine_data(estimates_600.loci, "GTseq")
snp600_long <- combine_data(estimates_600K, "600K")

# Combine all datasets into one
all_data_long <- bind_rows(snp60_long, GTseq_long, snp600_long)

# Function to identify outliers
is_outlier <- function(x) {
  if (length(x) > 0) {
    return(x < quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) | 
             x > quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE))
  } else {
    return(rep(FALSE, length(x)))
  }
}

# Add a column to identify outliers
all_data_long <- all_data_long %>%
  group_by(Estimator_Type, Known.relationship, Dataset) %>%
  mutate(is_outlier = is_outlier(Genetic_Relatedness_Estimate))

# Create plots for each Estimator Type with lines for each dataset
ggplot(all_data_long, aes(x = Known.relationship, y = Genetic_Relatedness_Estimate, fill = Known.relationship)) +
  geom_boxplot() +
  geom_text(data = filter(all_data_long, is_outlier),
            aes(label = Number), 
            position = position_jitter(width = 0.25, height = 0), 
            color = "black", 
            size = 3) +
  facet_grid(Estimator_Type ~ Dataset, scales = "free_y") +
  labs(title = "Genetic Relatedness Estimates by Estimator Type", 
       x = "Known Relationship", 
       y = "Genetic Relatedness Estimate") +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  theme_bw()



###################################################################################################################
###################################################################################################################
###################################################################################################################

# Reorder the 'Dataset' column
#all_data_long$Dataset <- factor(all_data_long$Dataset, 
#   levels = c("600 Loci", "Medium-Density", "High-Density"))

# Plot the data
ggplot(all_data_long, aes(x = Known.relationship, y = Genetic_Relatedness_Estimate, fill = Known.relationship)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +  # Keeps outlier dots but removes labels
  facet_grid(Estimator_Type ~ Dataset, scales = "fixed") +  # Ensures consistent Y-axis scaling
  scale_fill_manual(values = c("PO" = "darkgreen", "Siblings" = "lightgreen"),
                    labels = c("PO" = "Doe-Fetus", "Siblings" = "Fetus-Fetus")) +  
  scale_x_discrete(labels = c("PO" = "Doe-Fetus", "Siblings" = "Fetus-Fetus")) +  # Update x-axis labels
  labs(x = "Known Relationship", y = "Genetic Relatedness Estimate") +  # No title
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    strip.text.x = element_text(size = 24),
    strip.text.y = element_text(size = 24)
  ) +
  theme_bw()

####
library(dplyr)

# Relevel Dataset in desired order
all_data_long <- all_data_long %>%
  mutate(Dataset = factor(Dataset, levels = c("GTseq", "60K", "600K")))

# Plot
ggplot(all_data_long, aes(x = Known.relationship, y = Genetic_Relatedness_Estimate, fill = Known.relationship)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  facet_grid(
    Estimator_Type ~ Dataset, 
    scales = "fixed",
    labeller = labeller(
      Dataset = c("GTseq" = "600-Loci", "60K" = "Medium-Density", "600K" = "High-Density")
    )
  ) +
  scale_fill_manual(
    values = c("PO" = "darkgreen", "Siblings" = "lightgreen"),
    labels = c("PO" = "Female-Fetus", "Siblings" = "Fetus-Fetus")
  ) +
  scale_x_discrete(
    labels = c("PO" = "Female-Fetus", "Siblings" = "Fetus-Fetus")
  ) +
  guides(fill = guide_legend(title = "Known Relationship")) +
  labs(x = "Known Relationship", y = "Genetic Relatedness Estimate") +
  theme(
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    strip.text.x = element_text(size = 24),
    strip.text.y = element_text(size = 24),
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  theme_bw()


#######################
library(dplyr)

# Relevel Dataset in desired order
all_data_long <- all_data_long %>%
  mutate(Dataset = factor(Dataset, levels = c("GTseq", "60K", "600K")))

# Plot
Coancestry <- ggplot(all_data_long, aes(x = Known.relationship, y = Genetic_Relatedness_Estimate, fill = Known.relationship)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  facet_grid(
    Estimator_Type ~ Dataset, 
    scales = "fixed",
    labeller = labeller(
      Dataset = c("GTseq" = "600-Loci", "60K" = "Medium-Density", "600K" = "High-Density")
    )
  ) +
  scale_fill_manual(
    values = c("PO" = "darkgreen", "Siblings" = "lightgreen"),
    labels = c("PO" = "Female-Fetus", "Siblings" = "Fetus-Fetus")
  ) +
  scale_x_discrete(
    labels = c("PO" = "Female-Fetus", "Siblings" = "Fetus-Fetus")
  ) +
  guides(fill = guide_legend(title = "Known Relationship")) +
  labs(x = "Known Relationship", y = "Genetic Relatedness Estimate") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),   # Keep X facet labels large
    strip.text.y = element_text(size = 14),   # Reduce size of Y facet labels
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.box = "horizontal"
  )

#ggsave(
# filename = "Figure_Name.tif",
#plot = Coancestry,
#device = "tiff",
#path = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/Publication.Figures/Scripts/G3 Code & Figures",
#dpi = 300,       # High resolution for publication
#width = 12,       # Adjust width as needed
#height = 9,      # Adjust height as needed
# units = "in"
#)

###################################################################################################################

# Summarize relatedness ranges
relatedness_summary <- all_data_long %>%
  group_by(Dataset, Estimator_Type, Known.relationship) %>%
  summarize(
    Min = min(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Max = max(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Mean = mean(Genetic_Relatedness_Estimate, na.rm = TRUE),
    SD = sd(Genetic_Relatedness_Estimate, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

# Preview summary
print(relatedness_summary)





#################################################################################################################3
# Define the problematic pairs
problem_pairs <- list(
  c("MN164214_F1", "MN164214"),
  c("MN164225_F1", "MN164225")
)

# Create a function to check if a row contains a problematic pair
is_problem_pair <- function(ind1, ind2, pairs) {
  any(sapply(pairs, function(p) {
    (ind1 == p[1] & ind2 == p[2]) | (ind1 == p[2] & ind2 == p[1])
  }))
}

# Apply the function to identify problematic rows
all_data_long$problem_pair <- mapply(is_problem_pair, 
                                     all_data_long$Individual1, 
                                     all_data_long$Individual2, 
                                     MoreArgs = list(pairs = problem_pairs))

# Subset into two dataframes
Problematicpairs <- all_data_long[all_data_long$problem_pair, ]
noprobpair_all_data_long <- all_data_long[!all_data_long$problem_pair, ]

# Drop the helper column
Problematicpairs$problem_pair <- NULL
noprobpair_all_data_long$problem_pair <- NULL

###################################################################################################################

# Summarize relatedness ranges
relatedness_summaryF_Noprob <- noprobpair_all_data_long %>%
  group_by(Dataset, Estimator_Type, Known.relationship) %>%
  summarize(
    Min = min(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Max = max(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Mean = mean(Genetic_Relatedness_Estimate, na.rm = TRUE),
    SD = sd(Genetic_Relatedness_Estimate, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )


# Summarize relatedness ranges
relatedness_summaryF_prob <- Problematicpairs %>%
  group_by(Dataset, Estimator_Type, Known.relationship) %>%
  summarize(
    Min = min(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Max = max(Genetic_Relatedness_Estimate, na.rm = TRUE),
    Mean = mean(Genetic_Relatedness_Estimate, na.rm = TRUE),
    SD = sd(Genetic_Relatedness_Estimate, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )





###################################################################################################################
####################################################################################################################

#unknown relatedness estimates

###################################################################################################################
####################################################################################################################
#600-Loci

# Read data
relatedness_estimates_600loci <- read_excel("C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/600.Individual.loci/COANCESTRY/600loci_relatedness.estimates.xlsx")
known_parentage_siblings_600loci <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")

# Create both pair orders (A-B and B-A) as character keys for matching
known_pairs <- known_parentage_siblings_600loci %>%
  mutate(pair_key_1 = paste(Individual1, Individual2, sep = "_"),
         pair_key_2 = paste(Individual2, Individual1, sep = "_")) %>%
  select(pair_key_1, pair_key_2) %>%
  pivot_longer(cols = everything(), values_to = "pair_key") %>%
  distinct(pair_key)

# Create pair keys in estimates dataframe
relatedness_estimates_600loci <- relatedness_estimates_600loci %>%
  mutate(pair_key = paste(Individual1, Individual2, sep = "_"))

# Filter to keep only unknown pairs
relatedness_estimates_unknown_600loci <- relatedness_estimates_600loci %>%
  filter(!(pair_key %in% known_pairs$pair_key))
###################################################################################################################
#############################################################################################################3######
# Read in the data
relatedness_estimates_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_relatedness.estimates.csv")
known_parentage_siblings_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")
relatedness_estimates_600K <- read_csv("C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Coancestry\\Output\\600QC15.20_relatedness.estimate.results_7.15.24.csv")
known_parentage_siblings_600K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")

### 60K ###
# Read data
relatedness_estimates_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_relatedness.estimates.csv")
known_parentage_siblings_60K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv")

# Create both pair orders as keys
known_pairs_60K <- known_parentage_siblings_60K %>%
  mutate(pair_key_1 = paste(Individual1, Individual2, sep = "_"),
         pair_key_2 = paste(Individual2, Individual1, sep = "_")) %>%
  pivot_longer(cols = starts_with("pair_key"), values_to = "pair_key") %>%
  distinct(pair_key)

# Add pair keys and filter
relatedness_estimates_unknown_60K <- relatedness_estimates_60K %>%
  mutate(pair_key = paste(Individual1, Individual2, sep = "_")) %>%
  filter(!(pair_key %in% known_pairs_60K$pair_key))


### 600K ###
# Read data
relatedness_estimates_600K <- read_csv("C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Coancestry/Output/600QC15.20_relatedness.estimate.results_7.15.24.csv")
known_parentage_siblings_600K <- read_csv("Objective.One_60Kvs.600K_Relatedness/Coancestry/60K_known.parentage.siblings.csv") # still using 60K-known pairs here?

# Create both pair orders as keys
known_pairs_600K <- known_parentage_siblings_600K %>%
  mutate(pair_key_1 = paste(Individual1, Individual2, sep = "_"),
         pair_key_2 = paste(Individual2, Individual1, sep = "_")) %>%
  pivot_longer(cols = starts_with("pair_key"), values_to = "pair_key") %>%
  distinct(pair_key)

# Add pair keys and filter
relatedness_estimates_unknown_600K <- relatedness_estimates_600K %>%
  mutate(pair_key = paste(Individual1, Individual2, sep = "_")) %>%
  filter(!(pair_key %in% known_pairs_600K$pair_key))

###############################################################################################################3###
###################################################################################################################
# Columns to evaluate
relatedness_methods <- c("TrioML", "Wang", "LynchLi", "LynchRd", "Ritland", "QuellerGt", "DyadML")

# Function to filter out rows with any value > 0.90
filter_below_threshold <- function(df, columns, threshold = 0.90) {
  df %>% filter(if_all(all_of(columns), ~ .x <= threshold | is.na(.x)))
}

# Function to get ranges
get_ranges <- function(df, columns) {
  df %>%
    select(any_of(columns)) %>%
    summarise(across(everything(), ~ paste0(round(min(.x, na.rm = TRUE), 4), " to ", round(max(.x, na.rm = TRUE), 4))))
}

# Apply to each dataset
filtered_600loci <- filter_below_threshold(relatedness_estimates_unknown_600loci, relatedness_methods)
filtered_60K     <- filter_below_threshold(relatedness_estimates_unknown_60K, relatedness_methods)
filtered_600K    <- filter_below_threshold(relatedness_estimates_unknown_600K, relatedness_methods)

range_600loci <- get_ranges(filtered_600loci, relatedness_methods)
range_60K     <- get_ranges(filtered_60K, relatedness_methods)
range_600K    <- get_ranges(filtered_600K, relatedness_methods)

# Print results
cat("600 Loci Ranges (after removing values > 0.90):\n")
print(range_600loci)

cat("\n60K Ranges (after removing values > 0.90):\n")
print(range_60K)

cat("\n600K Ranges (after removing values > 0.90):\n")
print(range_600K)




###################################################################################################################
#getting average
get_summary_stats <- function(df, columns) {
  df %>%
    select(any_of(columns)) %>%
    summarise(across(everything(), list(
      Min = ~ round(min(.x, na.rm = TRUE), 4),
      Max = ~ round(max(.x, na.rm = TRUE), 4),
      Mean = ~ round(mean(.x, na.rm = TRUE), 4)
    ), .names = "{.col}_{.fn}"))
}

summary_600loci <- get_summary_stats(filtered_600loci, relatedness_methods)
summary_60K     <- get_summary_stats(filtered_60K, relatedness_methods)
summary_600K    <- get_summary_stats(filtered_600K, relatedness_methods)

cat("600 Loci Summary Stats (after removing values > 0.90):\n")
print(summary_600loci)

cat("\n60K Summary Stats (after removing values > 0.90):\n")
print(summary_60K)

cat("\n600K Summary Stats (after removing values > 0.90):\n")
print(summary_600K)


###################################################################################################################

get_summary_stats <- function(df, columns) {
  df %>%
    select(any_of(columns)) %>%
    summarise(across(everything(), list(
      Min = ~ round(min(.x, na.rm = TRUE), 4),
      Max = ~ round(max(.x, na.rm = TRUE), 4),
      Mean = ~ round(mean(.x, na.rm = TRUE), 4)
    ), .names = "{.col}_{.fn}"))
}

# Only retain unrelated pairs based on low relatedness estimates
unrelated_600loci <- filter_below_threshold(relatedness_estimates_unknown_600loci, relatedness_methods)
unrelated_60K     <- filter_below_threshold(relatedness_estimates_unknown_60K, relatedness_methods)
unrelated_600K    <- filter_below_threshold(relatedness_estimates_unknown_600K, relatedness_methods)

unrelated_summary_600loci <- get_summary_stats(unrelated_600loci, relatedness_methods)
unrelated_summary_60K     <- get_summary_stats(unrelated_60K, relatedness_methods)
unrelated_summary_600K    <- get_summary_stats(unrelated_600K, relatedness_methods)



cat("600 Loci Summary Stats for Unrelated Pairs (<= 0.90):\n")
print(unrelated_summary_600loci)

cat("\n60K Summary Stats for Unrelated Pairs (<= 0.90):\n")
print(unrelated_summary_60K)

cat("\n600K Summary Stats for Unrelated Pairs (<= 0.90):\n")
print(unrelated_summary_600K)



###################################################################################################################

# All pairs
relatedness_means <- relatedness_summary %>%
  select(Dataset, Estimator_Type, Known.relationship, Mean)

# Non-problematic pairs
relatedness_means_noprob <- relatedness_summaryF_Noprob %>%
  select(Dataset, Estimator_Type, Known.relationship, Mean)

# Problematic pairs
relatedness_means_prob <- relatedness_summaryF_prob %>%
  select(Dataset, Estimator_Type, Known.relationship, Mean)

# Preview
print(relatedness_means)
print(relatedness_means_noprob)
print(relatedness_means_prob)
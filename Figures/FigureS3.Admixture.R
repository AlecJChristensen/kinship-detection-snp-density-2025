# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the data
Loci600 <- read.table("C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\600.Individual.loci\\Admixture\\Random600.pilotA.2.Q")
MediumDensity <- read.table("C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\data_output\\Admixture\\60K\\AMix60QC6.11.24.2.Q")
HighDensity <- read.table("C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/data_output/Admixture/600K/AMix600QC6.20.24.2.Q")

# Function to rename columns dynamically
rename_admixture_data <- function(df, dataset_name) {
  colnames(df) <- paste0("Cluster", seq_len(ncol(df)))
  df$Dataset <- dataset_name
  df$Individual <- seq_len(nrow(df))
  return(df)
}

# Rename and prepare datasets
Loci600 <- rename_admixture_data(Loci600, "600 Loci")
MediumDensity <- rename_admixture_data(MediumDensity, "Medium-Density")
HighDensity <- rename_admixture_data(HighDensity, "High-Density")

# Combine datasets and sort individuals by Cluster2 proportion within each dataset
combined_df <- bind_rows(Loci600, MediumDensity, HighDensity) %>%
  mutate(Dataset = factor(Dataset, levels = c("600 Loci", "Medium-Density", "High-Density"))) %>%
  group_by(Dataset) %>%
  arrange(desc(Cluster2), .by_group = TRUE) %>%
  mutate(Sorted_Individual = row_number()) %>%
  ungroup()

# Convert to long format for plotting
combined_df_long <- combined_df %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion")

# Order clusters by mean proportion
cluster_order <- combined_df_long %>%
  group_by(Cluster) %>%
  summarise(mean_prop = mean(Proportion)) %>%
  arrange(mean_prop) %>%
  pull(Cluster)

combined_df_long$Cluster <- factor(combined_df_long$Cluster, levels = cluster_order)

# Define color scheme for clusters
cluster_colors <- c("Cluster1" = "#DB6912",
                    "Cluster2" = "#1B9E77")

# Plot
admixture <- ggplot(combined_df_long, aes(x = Sorted_Individual, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Dataset, scales = "free_x",
             labeller = labeller(Dataset = c(
               "600 Loci" = "(a) 600-Loci",
               "Medium-Density" = "(b) Medium-Density",
               "High-Density" = "(c) High-Density"
             ))) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 20) +
  labs(x = "Individual #", y = "Ancestry", fill = "Cluster") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 20),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

ggsave(
  filename = "Figure_Name.tif",
  plot = admixture,
  device = "tiff",
  path = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/Publication.Figures/Scripts/G3 Code & Figures",
  dpi = 300,       # High resolution for publication
  width = 11,       # Adjust width as needed
  height = 6,      # Adjust height as needed
  units = "in"
)

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create the dataset
data <- data.frame(
  Panel = rep(c("600 Loci", "Medium-Density", "High-Density"), each = 2),
  Relationship = rep(c("Doe-Fetus", "Fetus-Fetus"), times = 3),
  PLINK = c(0.929, 1, 0.5, 1, 0.536, 1),
  Sequoia = c(0.929, 1, 0.929, 1, 0.929, 1),
  COLONY = c(0.929, 1, 1, 1, 0.929, 1)
)

# Reshape data for ggplot
long_data <- pivot_longer(data, cols = PLINK:COLONY, names_to = "Method", values_to = "Proportion")

# Rename PLINK method
long_data$Method <- factor(long_data$Method, levels = c("PLINK", "Sequoia", "COLONY"),
                           labels = c("PLINK (KING-robust)", "Sequoia", "COLONY"))

# Order the Panel categories correctly
long_data$Panel <- factor(long_data$Panel, levels = c("600 Loci", "Medium-Density", "High-Density"))

# Define colors for Relationship
relationship_colors <- c("Doe-Fetus" = "darkgreen", "Fetus-Fetus" = "lightgreen")

# Create the plot
ggplot(long_data, aes(x = Panel, y = Proportion, fill = Relationship)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Method, nrow = 1) + 
  scale_fill_manual(values = relationship_colors) +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "black", linewidth = 1) +
  labs(x = "Density of Panel", y = "Proportion of Correct Known Pairs") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),  # Vertical labels
    legend.position = "bottom"  # Move legend underneath
  )



###################################################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create the dataset with proportion and matches
data <- data.frame(
  Panel = rep(c("600 Loci", "Medium-Density", "High-Density"), each = 2),
  Relationship = rep(c("Doe-Fetus", "Fetus-Fetus"), times = 3),
  PLINK = c(0.929, 1, 0.5, 1, 0.536, 1),
  Sequoia = c(0.929, 1, 0.929, 0, 0.929, 0),
  COLONY = c(0.929, 1, 1, 1, 0.929, 1),
  Matches_PLINK = c(26, 15, 14, 15, 15, 15),
  Matches_Sequoia = c(26, 15, 26, 15, 26, 15),
  Matches_COLONY = c(26, 15, 28, 15, 26, 15)
)

# Reshape data for ggplot
long_data <- pivot_longer(data, cols = c(PLINK, Sequoia, COLONY), 
                          names_to = "Method", values_to = "Proportion")

long_matches <- pivot_longer(data, cols = c(Matches_PLINK, Matches_Sequoia, Matches_COLONY), 
                             names_to = "Matches_Method", values_to = "Matches")

# Extract method names to match proportions
long_matches$Matches_Method <- gsub("Matches_", "", long_matches$Matches_Method)

# Merge proportion and matches
long_data <- left_join(long_data, long_matches, by = c("Panel", "Relationship", "Method" = "Matches_Method"))

# Rename PLINK method
long_data$Method <- factor(long_data$Method, levels = c("PLINK", "Sequoia", "COLONY"),
                           labels = c("PLINK (KING-robust)", "Sequoia", "COLONY"))

# Order the Panel categories correctly
long_data$Panel <- factor(long_data$Panel, levels = c("600 Loci", "Medium-Density", "High-Density"))

# Define colors for Relationship
relationship_colors <- c("Doe-Fetus" = "darkgreen", "Fetus-Fetus" = "lightgreen")

# Create the plot
ggplot(long_data, aes(x = Panel, y = Proportion, fill = Relationship)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Method, nrow = 1) + 
  scale_fill_manual(values = relationship_colors) +
  geom_text(aes(label = Matches), position = position_dodge(width = 0.9), 
            vjust = 1.5, color = "white", size = 5) +  # White text for matches
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "black", linewidth = 1) +
  labs(x = "Density of Panel", y = "Proportion of Correct Known Pairs") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),  # Vertical labels
    legend.position = "bottom"  # Move legend underneath
  )

###################################################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create the dataset with proportion and matches
data <- data.frame(
  Panel = rep(c("600-Loci", "Medium-Density", "High-Density"), each = 2),
  Relationship = rep(c("Female-Fetus", "Fetus-Fetus"), times = 3),
  PLINK = c(0.929, 1, 0.5, 1, 0.536, 1),
  Sequoia = c(0.929, 1, 0.929, 1, 0.929, 1),
  COLONY = c(0.929, 1, 1, 1, 0.929, 1),
  Matches_PLINK = c(26, 15, 14, 15, 15, 15),
  Matches_Sequoia = c(26, 15, 26, 15, 26, 15),
  Matches_COLONY = c(26, 15, 28, 15, 26, 15)
)

# Reshape data for ggplot
long_data <- pivot_longer(data, cols = c(PLINK, Sequoia, COLONY), 
                          names_to = "Method", values_to = "Proportion")

long_matches <- pivot_longer(data, cols = c(Matches_PLINK, Matches_Sequoia, Matches_COLONY), 
                             names_to = "Matches_Method", values_to = "Matches")

# Extract method names to match proportions
long_matches$Matches_Method <- gsub("Matches_", "", long_matches$Matches_Method)

# Merge proportion and matches
long_data <- left_join(long_data, long_matches, by = c("Panel", "Relationship", "Method" = "Matches_Method"))

# Rename methods with prefixes
long_data$Method <- factor(long_data$Method, levels = c("PLINK", "Sequoia", "COLONY"),
                           labels = c("(a) PLINK (KING-robust)", "(b) Sequoia", "(c) COLONY"))

# Order the Panel categories correctly
long_data$Panel <- factor(long_data$Panel, levels = c("600-Loci", "Medium-Density", "High-Density"))

# Define colors for Relationship
relationship_colors <- c("Female-Fetus" = "darkgreen", "Fetus-Fetus" = "lightgreen")

# Create the plot
Sequoia <- ggplot(long_data, aes(x = Panel, y = Proportion, fill = Relationship)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Method, nrow = 1) + 
  scale_fill_manual(values = relationship_colors) +
  geom_text(aes(label = Matches), position = position_dodge(width = 0.9), 
            vjust = 1.5, color = "white", size = 5) +  # White text for matches
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "black", linewidth = 1) +
  labs(x = "Density of Panel", y = "Proportion of Correct Known Pairs") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),  # Vertical labels
    legend.position = "bottom"  # Move legend underneath
  )


ggsave(
  filename = "Figure_Name.tif",
  plot = Sequoia,
  device = "tiff",
  path = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/Publication.Figures/Scripts/G3 Code & Figures",
  dpi = 300,       # High resolution for publication
  width = 12,       # Adjust width as needed
  height = 9,      # Adjust height as needed
  units = "in"
)

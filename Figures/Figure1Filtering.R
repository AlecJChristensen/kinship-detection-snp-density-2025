
###################################################################################################################
# Load libraries
library(ggplot2)
library(dplyr)

# Create dataset
data <- data.frame(
  Panel = rep(c("Medium-Density", "High-Density"), each = 9),
  SNP_Missingness = rep(c("15%", "15%", "15%", "20%", "20%", "20%", "25%", "25%", "25%"), 2),
  Individual_Missingness = rep(c("10%", "15%", "20%"), 6),
  SNPs_Retained = c(65594, 65594, 65594, 67457, 67457, 67457, 67099, 67099, 67099,
                    507864, 507864, 507864, 600866, 600866, 600866, 574217, 574217, 574217),
  Retained_Individuals = c(94, 95, 96, 94, 95, 96, 94, 95, 96,
                           91, 95, 96, 86, 95, 96, 79, 95, 96),
  Known_Pairs = c(26, 26, 26, 26, 26, 26, 26, 26, 26,
                  22, 25, 26, 2, 24, 24, 6, 24, 24),
  Proportion = c(0.929, 0.929, 0.929, 0.929, 0.929, 0.929, 0.929, 0.929, 0.929,
                 0.786, 0.893, 0.929, 0.071, 0.857, 0.857, 0.214, 0.857, 0.857)
)

# Modify Filtering_Threshold to remove prefixes "Medium-Density" and "High-Density"
data$Filtering_Threshold <- paste0(gsub("Medium-Density", "", data$Panel), " ",
                                   gsub("%", "", data$SNP_Missingness), ".", 
                                   gsub("%", "", data$Individual_Missingness), " (", 
                                   format(data$SNPs_Retained, big.mark = ","), ")")

# Reorder the Filtering_Threshold factor levels to ensure correct x-axis order
data$Filtering_Threshold <- factor(data$Filtering_Threshold, 
                                   levels = unique(data$Filtering_Threshold))

# Calculate Retained Percentage
data$Retained_Percentage <- paste0(round((data$Retained_Individuals / 96) * 100, 0), "%")

# Define colors for panels
panel_colors <- c("Medium-Density" = "#66C2A5", "High-Density" = "#1B4D3E")

# Ensure the "Medium-Density" panel appears first in the facet order
data$Panel <- factor(data$Panel, levels = c("Medium-Density", "High-Density"))

# Create plot
ggplot(data, aes(x = Filtering_Threshold, y = Proportion, fill = Panel)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = panel_colors) +
  # Known Pairs label inside the column
  geom_text(aes(label = Known_Pairs), vjust = 1.2, hjust = 0.5, color = "white", size = 3, fontface = "bold") +  
  # Retained Percentage label inside the column
  geom_text(aes(label = Retained_Percentage), vjust = 2.5, hjust = 0.5, color = "white", size = 3) +  
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "black", linewidth = 1) +
  facet_wrap(~ Panel, scales = "free_x", strip.position = "top") +  # Flip facet order
  labs(x = "Filtering Threshold (SNPs Retained)", 
       y = "Proportion of Correct Doe-Fetus Pairs") +
  theme_bw() +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),  # Adjust X-axis labels size
    strip.text.x = element_text(size = 16),  # Adjust facet title size
    axis.title.x = element_text(size = 18),  # Adjust X-axis title size
    axis.title.y = element_text(size = 18),  # Adjust Y-axis title size
    legend.position = "none",  # Remove legend since we are faceting
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines for clarity
    panel.grid.minor.x = element_blank()  # Remove vertical grid lines for clarity
  )


# Modify Filtering_Threshold to remove prefixes "Medium-Density" and "High-Density" where applicable
data$Filtering_Threshold <- ifelse(data$Panel == "Medium-Density", 
                                   paste0(gsub("%", "", data$SNP_Missingness), ".", 
                                          gsub("%", "", data$Individual_Missingness), 
                                          " (", format(data$SNPs_Retained, big.mark = ","), ")"),
                                   paste0(gsub("High-Density", "", data$Panel), " ", 
                                          gsub("%", "", data$SNP_Missingness), ".", 
                                          gsub("%", "", data$Individual_Missingness), 
                                          " (", format(data$SNPs_Retained, big.mark = ","), ")"))

# Reorder the Filtering_Threshold factor levels to ensure correct x-axis order
data$Filtering_Threshold <- factor(data$Filtering_Threshold, 
                                   levels = unique(data$Filtering_Threshold))
# Update panel labels
data$Panel <- recode(data$Panel, 
                     "Medium-Density" = "(a) Medium-Density", 
                     "High-Density" = "(b) High-Density")

# Update colors for new facet labels
panel_colors <- c("(a) Medium-Density" = "#66C2A5", "(b) High-Density" = "#1B4D3E")

# Create plot
filter <- ggplot(data, aes(x = Filtering_Threshold, y = Proportion, fill = Panel)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = panel_colors) +
  geom_text(aes(label = Known_Pairs), vjust = 1.2, hjust = 0.5, color = "white", size = 6, fontface = "bold") +  
  geom_text(aes(label = Retained_Percentage), vjust = 2.5, hjust = 0.5, color = "white", size = 6) +  
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "black", linewidth = 1) +
  facet_wrap(~ Panel, scales = "free_x", strip.position = "top") +  
  labs(x = "Filtering Threshold (SNPs Retained)", 
       y = "Proportion of Correct Female-Fetus Pairs") +
  theme_bw() +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),  
    strip.text.x = element_text(size = 16),  
    axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 18),  
    legend.position = "none",  
    panel.grid.major.x = element_blank(),  
    panel.grid.minor.x = element_blank()  
  )


###################################################################################################################
###################################################################################################################

###################################################################################################################


ggsave(
  filename = "Figure_Name.tif",
  plot = filter,
  device = "tiff",
  path = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/Publication.Figures/Scripts/G3 Code & Figures",
  dpi = 300,       # High resolution for publication
  width = 13,       # Adjust width as needed
  height = 10,      # Adjust height as needed
  units = "in"
)

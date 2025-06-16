# Load necessary library
library(ggplot2)

# Your data
Medium_Density <- c(0.5724, 0.57858, 0.58983, 0.61688, 0.63047, 0.65229, 0.66364, 0.69877, 0.73314, 0.74255)
High_Density <- c(0.58137, 0.57062, 0.59609, 0.61691, 0.63635, 0.65091, 0.67325, 0.6927, 0.73231, 0.76194)

# 600 Loci data
K_values <- 1:10
Loci_600 <- c(0.5375, 0.54484, 0.55298, 0.56361, 0.58077, 0.57343, 0.59863, 0.61723, 0.62621, 0.62548)

# Create a data frame for Medium-Density, High-Density, and 600 Loci
df <- data.frame(
  K = rep(K_values, 3),
  Cross_Validation_Error = c(Medium_Density, High_Density, Loci_600),
  Density_Type = rep(c("Medium-Density", "High-Density", "600-Loci"), each = length(K_values))
)

# Set the order of levels in the Density_Type factor
df$Density_Type <- factor(df$Density_Type, levels = c("600-Loci", "Medium-Density", "High-Density"))

# Plot using ggplot with bw theme
CV<- ggplot(df, aes(x = K, y = Cross_Validation_Error, color = Density_Type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = c("600-Loci" = "#B5D6C3", 
                                "Medium-Density" = "#66C2A5", 
                                "High-Density" = "#1B4D3E")) +
  labs(
    y = "Cross-Validation Error",
    x = "K value",
    color = "Density Type"
  ) +
  theme_bw() + # black and white theme
  theme(
    legend.position = "bottom",
    text = element_text(size = 20)  # Set minimum font size to 20
  )

ggsave(
  filename = "Figure_Name.tif",
  plot = CV,
  device = "tiff",
  path = "C:/Users/alec1/Desktop/CWDProject/Objective.One_60Kvs.600K_Relatedness/Scripts/Final_Scripts/List.for.Google.Doc/Publication.Figures/Scripts/G3 Code & Figures",
  dpi = 300,       # High resolution for publication
  width = 8,       # Adjust width as needed
  height = 6,      # Adjust height as needed
  units = "in"
)
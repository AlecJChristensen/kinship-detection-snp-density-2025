#Last Updated 7.23.24
library(vcfR)
library(adegenet)
library(stringr)
library(ggplot2)
library(dplyr)


###################################################################################################################
###################################################################################################################
###################################################################################################################
#60K
df <- read.csv("MNDNR_SamplesForSNPsProject_10182023_230pm.csv", na.strings=c("","NA"))

vcf_path60 <- "1670010002_LowThreshold_LowSNPQCVCF.vcf"
vcf60 <- read.vcfR(vcf_path60)
#loading in Emily Spreadsheet of individuals
Metadata<- read.csv("MNDNR_SamplesForSNPsProject_10182023_230pm.csv")
df$CWD_status <- toupper(df$CWD_status)  # Convert all values to uppercase

#turning Vcf to genind
OvSNP60_genind <- vcfR2genind(vcf60)
#Looking at Genind
OvSNP60_genind
gdf60K<-genind2df(OvSNP60_genind)

###################################

row_names60K <- rownames(gdf60K)

# Create a new data frame with the row names
row_names_df60K <- data.frame(row_names60K)

# Optionally, you can rename the column in the new data frame
colnames(row_names_df60K) <- "WHP_sample_id"

# Print or view the new data frame
print(row_names_df60K)

# Extract and edit the names in the "IID" column
row_names_df60K$WHP_sample_id <- str_replace(row_names_df60K$WHP_sample_id, "^[^_]*_(.*?)\\.CEL$", "\\1")
# Remove everything before the second underscore
row_names_df60K$WHP_sample_id <- str_extract(row_names_df60K$WHP_sample_id, "(?<=_).*")

new.merged60K <- merge(row_names_df60K, df, by = "WHP_sample_id", sort = FALSE)

new_row60K <- data.frame(WHP_sample_id = "Dup_MN84503",
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

# Append the new row to merged_df
merged.df60K <- rbind(new.merged60K, new_row60K)

OvSNP60_genind@pop <- as.factor(merged.df60K$permit_area)
#########################################################

#PCA
deer60K <- tab(OvSNP60_genind, freq=TRUE, NA.method="mean")
pca.deer60K <- dudi.pca(df = deer60K, center = TRUE, scale = FALSE, nf = 200, scannf = FALSE)
pca.deer60K
s.label(pca.deer60K$li)
s.class(pca.deer60K$li,fac=pop(OvSNP60_genind), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)

s.class(pca.deer60K$li,fac=pop(OvSNP60_genind), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)
add.scatter.eig(pca.deer60K$eig[1:50],3,1,2, ratio=.3)

#without labels
deer60K <- tab(OvSNP60_genind, freq=TRUE, NA.method="mean")
pca.deer60K <- dudi.pca(df = deer60K, center = TRUE, scale = FALSE, nf = 200, scannf = FALSE)
pca.deer60K
s.label(pca.deer60K$li, clabel=0)
s.class(pca.deer60K$li, fac=pop(OvSNP60_genind), col=transp(funky(15), .6), axesel=FALSE, cstar=0, clabel=0, cpoint=3)
add.scatter.eig(pca.deer60K$eig[1:50], 3, 1, 2, ratio=.3)
#######
# Generate the allele frequencies table
deer60K <- tab(OvSNP60_genind, freq = TRUE, NA.method = "mean")

# Perform PCA on the allele frequencies
pca.deer60K <- dudi.pca(df = deer60K, center = TRUE, scale = FALSE, nf = 200, scannf = FALSE)
pca.deer60K

# Specify custom colors for the groups
group_colors <- c("peachpuff2", "slateblue4", "seagreen", "darkgoldenrod3", "salmon4")  # Adjust these colors as needed

# PCA plot with labels
s.label(pca.deer60K$li)
s.class(pca.deer60K$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, cpoint = 3)
add.scatter.eig(pca.deer60K$eig[1:50], 3, 1, 2, ratio = 0.3)

# PCA plot without labels
s.label(pca.deer60K$li, clabel = 0)
s.class(pca.deer60K$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer60K$eig[1:50], 3, 1, 2, ratio = 0.3)


###################################################################################################################
###################################################################################################################
###################################################################################################################
#600K

#600KQC15.20
load("~/R.Work/com.df.RData")

colnames(com.df) <- gsub("\\.", "-", colnames(com.df))
genind.600K <- df2genind(com.df, sep = "\t")



row_names600 <- data.frame(individual_ID = com.df$individual_ID)

# Create a new data frame with the row names
row_names_df.600 <- data.frame(row_names600)

# Optionally, you can rename the column in the new data frame
colnames(row_names_df.600) <- "WHP_sample_id"

# Print or view the new data frame
print(row_names_df.600)

# Extract and edit the names in the "IID" column
row_names_df.600$WHP_sample_id <- str_replace(row_names_df.600$WHP_sample_id, "^[^_]*_(.*?)\\.CEL$", "\\1")
# Remove everything before the second underscore
row_names_df.600$WHP_sample_id <- str_extract(row_names_df.600$WHP_sample_id, "(?<=_).*")

new.merged.600 <- merge(row_names_df.600, df, by = "WHP_sample_id", sort = FALSE)

new_row <- data.frame(WHP_sample_id = "Dup_MN84503",
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

# Append the new row to merged_df
merged.df.600 <- rbind(new.merged.600, new_row)

genind.600K@pop <- as.factor(merged.df.600$permit_area)


# Specify colors for each group
group_colors <- c("peachpuff2", "slateblue4", "seagreen", "darkgoldenrod3", "salmon4")


#PCA
deer600 <- tab(genind.600K, freq=TRUE, NA.method="mean")
pca.deer.600 <- dudi.pca(df = deer600, center = TRUE, scale = FALSE, nf = 12000, scannf = FALSE)
pca.deer.600
s.label(pca.deer.600$li)
s.class(pca.deer.600$li,fac=pop(genind.600K), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)

s.class(pca.deer.600$li,fac=pop(genind.600K), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)
add.scatter.eig(pca.deer.600$eig[1:50],3,1,2, ratio=.3)

#without labels
deer600 <- tab(genind.600K, freq=TRUE, NA.method="mean")
pca.deer.600 <- dudi.pca(df = deer600, center = TRUE, scale = FALSE, nf = 1300, scannf = FALSE)
pca.deer.600
s.label(pca.deer.600$li, clabel=0)
s.class(pca.deer.600$li, fac=pop(genind.600K), col=transp(funky(15), .6), axesel=FALSE, cstar=0, clabel=0, cpoint=3)
add.scatter.eig(pca.deer.600$eig[1:50], 3, 1, 2, ratio=.3)





# Specify custom colors for the groups
group_colors <- c("peachpuff2", "slateblue4", "seagreen", "darkgoldenrod3", "salmon4")  # Adjust these colors as needed

# PCA plot with labels
s.label(pca.deer.600$li)
s.class(pca.deer.600$li, fac = pop(genind.600K), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, cpoint = 3)
add.scatter.eig(pca.deer.600$eig[1:50], 3, 1, 2, ratio = 0.3)

# PCA plot without labels
s.label(pca.deer.600$li, clabel = 0)
s.class(pca.deer.600$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer.600$eig[1:50], 3, 1, 2, ratio = 0.3)

###################################################################################################################
###################################################################################################################
###################################################################################################################

###################################################################################################################
###################################################################################################################
###################################################################################################################
#600 Loci
df <- read.csv("MNDNR_SamplesForSNPsProject_10182023_230pm.csv", na.strings=c("","NA"))
vcf_path600loci <-"Random600.pilot.vcf"
vcf600loci <- read.vcfR(vcf_path600loci)

#turning Vcf to genind
loci600_genind <- vcfR2genind(vcf600loci)
#Looking at Genind
loci600_genind
gdf600L<-genind2df(loci600_genind)


###################################

row_names.600L <- rownames(gdf600L)

# Create a new data frame with the row names
row_names_df.600L <- data.frame(row_names.600L)

# Optionally, you can rename the column in the new data frame
colnames(row_names_df.600L) <- "WHP_sample_id"

# Print or view the new data frame
print(row_names_df.600L)

# Extract and edit the names in the "IID" column
row_names_df.600L$WHP_sample_id <- str_replace(row_names_df.600L$WHP_sample_id, "^[^_]*_(.*?)\\.CEL$", "\\1")
# Remove everything before the second underscore
row_names_df.600L$WHP_sample_id <- str_extract(row_names_df.600L$WHP_sample_id, "(?<=_).*")

new.merged.600L<- merge(row_names_df.600L, df, by = "WHP_sample_id", sort = FALSE)

new_row.600L <- data.frame(WHP_sample_id = "Dup_MN84503",
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

# Append the new row to merged_df
merged.df.600L <- rbind(new.merged.600L, new_row)

loci600_genind@pop <- as.factor(merged.df.600L$permit_area)
#########################################################

#PCA
deer600loci <- tab(loci600_genind, freq=TRUE, NA.method="mean")
pca.deer600loci <- dudi.pca(df = deer600loci, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE)
pca.deer600loci
s.label(pca.deer600loci$li)
s.class(pca.deer600loci$li,fac=pop(loci600_genind), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)

s.class(pca.deer600loci$li,fac=pop(loci600_genind), col=transp(funky(15),.6), axesel=FALSE,cstar=0, cpoint=3)
add.scatter.eig(loci600_genind$eig[1:50],3,1,2, ratio=.3)

#without labels
deer600loci <- tab(loci600_genind, freq=TRUE, NA.method="mean")
pca.deer600loci <- dudi.pca(df = deer600loci, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE)
pca.deer600loci
s.label(pca.deer600loci$li, clabel=0)
s.class(pca.deer600loci$li, fac=pop(OvSNP60_genind), col=transp(funky(15), .6), axesel=FALSE, cstar=0, clabel=0, cpoint=3)
add.scatter.eig(loci600_genind$eig[1:50], 3, 1, 2, ratio=.3)
#######
# Generate the allele frequencies table
deer600loci <- tab(loci600_genind, freq = TRUE, NA.method = "mean")

# Perform PCA on the allele frequencies
pca.deer600loci <- dudi.pca(df = deer600loci, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE)
pca.deer600loci

# Specify custom colors for the groups
group_colors <- c("peachpuff2", "slateblue4", "seagreen", "darkgoldenrod3", "salmon4")  # Adjust these colors as needed

# PCA plot with labels
s.label(pca.deer600loci$li)
s.class(pca.deer600loci$li, fac = pop(pca.deer600loci), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, cpoint = 3)
add.scatter.eig(pca.deer600loci$eig[1:50], 3, 1, 2, ratio = 0.3)

# PCA plot without labels
s.label(pca.deer600loci$li, clabel = 0)
s.class(pca.deer600loci$li, fac = pop(loci600_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(loci600_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer600loci$eig[1:50], 3, 1, 2, ratio = 0.3)



################################################################################
save.image(file = "PCA_600loci.60K.600K.RData")
################################################################################

# PCA plot without labels
s.label(pca.deer60K$li, clabel = 0)
s.class(pca.deer60K$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer60K$eig[1:50], 3, 1, 2, ratio = 0.3)


################################################################################

#600K
# PCA plot without labels
s.label(pca.deer.600$li, clabel = 0)
s.class(pca.deer.600$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer.600$eig[1:50], 3, 1, 2, ratio = 0.3)
################################################################################
#600Loci
# PCA plot without labels
s.label(pca.deer600loci$li, clabel = 0)
s.class(pca.deer600loci$li, fac = pop(loci600_genind), col = transp(group_colors, 0.6), axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(loci600_genind)), fill = group_colors, title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer600loci$eig[1:50], 3, 1, 2, ratio = 0.3)
################################################################################




################################################################################
# 60K Panel
pca_60K_plot <- recordPlot()
s.label(pca.deer60K$li, clabel = 0)
s.class(pca.deer60K$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), 
        axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, 
       title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer60K$eig[1:50], 3, 1, 2, ratio = 0.3)
pca_60K_plot <- recordPlot()
################################################################################
# 600K Panel
pca_600K_plot <- recordPlot()
s.label(pca.deer.600$li, clabel = 0)
s.class(pca.deer.600$li, fac = pop(OvSNP60_genind), col = transp(group_colors, 0.6), 
        axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(OvSNP60_genind)), fill = group_colors, 
       title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer.600$eig[1:50], 3, 1, 2, ratio = 0.3)
pca_600K_plot <- recordPlot()
################################################################################
# 600 Loci Panel
pca_600Loci_plot <- recordPlot()
s.label(pca.deer600loci$li, clabel = 0)
s.class(pca.deer600loci$li, fac = pop(loci600_genind), col = transp(group_colors, 0.6), 
        axesel = FALSE, cstar = 0, clabel = 0, cpoint = 3)
legend("topright", legend = levels(pop(loci600_genind)), fill = group_colors, 
       title = "Groups", cex = 0.8, inset = c(-0.3, 0), xpd = TRUE)
add.scatter.eig(pca.deer600loci$eig[1:50], 3, 1, 2, ratio = 0.3)
pca_600Loci_plot <- recordPlot()


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
library(ggplot2)
library(dplyr)

# Convert PCA results into a single dataframe
pca_60K_df <- as.data.frame(pca.deer60K$li) %>%
  mutate(Group = pop(OvSNP60_genind), Dataset = "60K")

pca_600K_df <- as.data.frame(pca.deer.600$li) %>%
  mutate(Group = pop(OvSNP60_genind), Dataset = "600K")

pca_600Loci_df <- as.data.frame(pca.deer600loci$li) %>%
  mutate(Group = pop(loci600_genind), Dataset = "600 Loci")

# Combine into one dataframe
pca_combined_df <- bind_rows(pca_60K_df, pca_600K_df, pca_600Loci_df)
################################################################################
################################################################################
################################################################################
group_colors <- c("#984EA3",  # Red  #984EA3
                  "#377EB8",  # Blue  
                  "#4DAF4A",  # Green  
                  "#E41A1C",  # Purple  #E41A1C
                  "#A65628")  # Brown
################################################################################
# Rename and reorder Dataset values
pca_combined_df <- pca_combined_df %>%
  mutate(Dataset = recode(Dataset, "60K" = "Medium-Density", "600K" = "High-Density")) %>%
  mutate(Dataset = factor(Dataset, levels = c("600 Loci", "Medium-Density", "High-Density")))  # Set the order

# Plot with updated order
ggplot(pca_combined_df, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~Dataset) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "PC1", y = "PC2", color = "Groups")





library(ggplot2)
library(dplyr)

# Convert PCA results into a single dataframe
pca_60K_df <- as.data.frame(pca.deer60K$li) %>%
  mutate(Group = pop(OvSNP60_genind), Dataset = "60K")

pca_600K_df <- as.data.frame(pca.deer.600$li) %>%
  mutate(Group = pop(OvSNP60_genind), Dataset = "600K")

pca_600Loci_df <- as.data.frame(pca.deer600loci$li) %>%
  mutate(Group = pop(loci600_genind), Dataset = "600 Loci")

# Combine into one dataframe
pca_combined_df <- bind_rows(pca_60K_df, pca_600K_df, pca_600Loci_df)

################################################################################
################################################################################
################################################################################
group_colors <- c("#984EA3",  # Red  #984EA3
                  "#377EB8",  # Blue  
                  "#4DAF4A",  # Green  
                  "#E41A1C",  # Purple  #E41A1C
                  "#A65628")  # Brown
################################################################################
# Rename and reorder Dataset values
pca_combined_df <- pca_combined_df %>%
  mutate(Dataset = recode(Dataset, 
                          "60K" = "(b) Medium-Density", 
                          "600K" = "(c) High-Density", 
                          "600 Loci" = "(a) 600-Loci")) %>%
  mutate(Dataset = factor(Dataset, levels = c("(a) 600-Loci", "(b) Medium-Density", "(c) High-Density")))  # Set the order

# Plot with updated order and increased font sizes
ggplot(pca_combined_df, aes(x = Axis1, y = Axis2, color = Group)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~Dataset) +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 16),      # Increases overall text size
        axis.title = element_text(size = 18), # Increases axis titles
        axis.text = element_text(size = 14),  # Increases axis labels
        strip.text = element_text(size = 18), # Increases facet labels
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  labs(x = "PC1", y = "PC2", color = "Groups")




pca_herd <- factoextra::fviz_pca_ind(pca.deer60K,
                                     label = "none", # color by groups
                                     addEllipses = TRUE, # Concentration ellipses
                                     ellipse.type = "confidence",
                                     legend.title = "Herd",
                                     repel = TRUE)

pca_herd600K <- factoextra::fviz_pca_ind(pca.deer.600,
                                         label = "none", # color by groups
                                         addEllipses = TRUE, # Concentration ellipses
                                         ellipse.type = "confidence",
                                         legend.title = "Herd",
                                         repel = TRUE)

pca_herd600loci <- factoextra::fviz_pca_ind(pca.deer600loci,
                                            label = "none", # color by groups
                                            addEllipses = TRUE, # Concentration ellipses
                                            ellipse.type = "confidence",
                                            legend.title = "Herd",
                                            repel = TRUE)
###################################################
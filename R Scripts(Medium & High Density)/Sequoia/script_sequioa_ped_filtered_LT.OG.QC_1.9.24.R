#3)	Sequioa of the different filtered data sets
#Wild Cervid Project
#Obj.One 
#Last Updated 7.23.24

#libraries
library(sequoia)

###################################################################################################################
#og data file
#loading in data file

#60K
infile60 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\PLINKSNP\\OvSNP60ped.ped", InFormat="ped")
LH60K <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\LH60_RealMaster.csv")

#600K
LH600K <- read.csv("C:/Users/alec1/Desktop/SNP_data/scripts/kn_LH600.csv")
infile600 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\PLINKSNP\\OvSNP600ped.ped", InFormat="ped")

AgePriors60 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#600K
AgePriors600 <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#Regular run off first originial vcfs
Ped60.175KD <- GetMaybeRel(GenoM = infile60,
                           LifeHistData = LH60K,
                           AgePrior = AgePriors60,
                           Module = 'ped',
                           Err = 0.175)# 34Po

Ped600.175KD <- GetMaybeRel(GenoM = infile600,
                            LifeHistData = LH600K,
                            AgePrior = AgePriors600,
                            Module = 'ped',
                            Err = 0.175)# 24PO

###################################################################################################################
#LOW THRESHOLD #2
#loading in data
#60K
LH60K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\Low.Threshold60.csv")
infile60LT <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\OvSNP60LTped.ped", InFormat="ped")

#600K
LH600K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\Low.Threshold600.csv")
infile600LT <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\OvSNP600LTped.ped", InFormat="ped")
#60K
AgePriors60LT <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K96,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#600K
AgePriors600LT <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K96,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#60K
Ped60.175LT <- GetMaybeRel(GenoM = infile60LT,
                           LifeHistData = LH60K96,
                           AgePrior = AgePriors60LT,
                           Module = 'ped',
                           Err = 0.175)#34 parent-offspring

Ped600.175LT <- GetMaybeRel(GenoM = infile600LT,
                            LifeHistData = LH600K96,
                            AgePrior = AgePriors600LT,
                            Module = 'ped',
                            Err = 0.175)#33 parent-offspring
###################################################################################################################
#QUALITY CONTROL #3

#data
#60K
LH60K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\Low.Threshold60.csv")
infile60QC <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\3.0_seq.three.files.OG.LT.QC_7.22.24\\Input.Files\\OvSNP60QCped.ped", InFormat="ped")

#600K
LH600K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\Low.Threshold600.csv")
infile600QC <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\3.0_seq.three.files.OG.LT.QC_7.22.24\\Input.Files\\OvSNP600QCped.ped", InFormat="ped")

#age priors
#60K
AgePriors60QC <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K96,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#600K
AgePriors600QC <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K96,
  MinAgeParent = 2,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

#60K parental analysis
Ped60.175QC <- GetMaybeRel(GenoM = infile60QC,
                           LifeHistData = LH60K96,
                           AgePrior = AgePriors60QC,
                           Module = 'ped',
                           Err = 0.175)#31 PO
#600K parental analysis
Ped600.175QC <- GetMaybeRel(GenoM = infile600QC,
                            LifeHistData = LH600K96,
                            AgePrior = AgePriors600QC,
                            Module = 'ped',
                            Err = 0.175)#20 PO
###################################################################################################################
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.175KD$MaybeRel, file = paste0(output_dir, "Ped60.175KD.csv"), row.names = FALSE)
write.csv(Ped600.175KD$MaybeRel, file = paste0(output_dir, "Ped600.175KD.csv"), row.names = FALSE)
write.csv(Ped60.175LT$MaybeRel, file = paste0(output_dir, "Ped60.175LT.csv"), row.names = FALSE)
write.csv(Ped600.175LT$MaybeRel, file = paste0(output_dir, "Ped600.175LT.csv"), row.names = FALSE)
write.csv(Ped60.175QC$MaybeRel, file = paste0(output_dir, "Ped60.175QC.csv"), row.names = FALSE)
write.csv(Ped600.175QC$MaybeRel, file = paste0(output_dir, "Ped600.175QC.csv"), row.names = FALSE)




###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
#sequioa with different filters QC
#Wild Cervid Project
#Obj.One 
#Last Updated 7.23.24

library(sequoia)

###################################################################################################################
####################################
##################

#60K

############################

#Main split of the three vcf's
#lifehistory data
LH60K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\LH60_RealMaster.csv")
#ped file for low threshold
infile60LT <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\OvSNP60LTped.ped", InFormat="ped")
#60K
AgePriors60LT <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K96,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

Ped60.175LT <- GetMaybeRel(GenoM = infile60LT,
                           LifeHistData = LH60K96,
                           AgePrior = AgePriors60LT,
                           Module = 'ped',
                           Err = 0.175)
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.175LT$MaybeRel, file = paste0(output_dir, "Ped60.175LT.csv"), row.names = FALSE)


#########################
#Low SNP Quality Control
#60K
LH60K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\LH60_RealMaster.csv")
infile60QC <- GenoConvert(InFile = "C:/Users/alec1/Desktop/SNP_data/data/LowThresholds_SNPQC/OvSNP60QCped.ped", InFormat="ped")

#600K
LH600K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\LH600_RealMaster.csv")
infile600QC <- GenoConvert(InFile = "C:/Users/alec1/Desktop/SNP_data/data/LowThresholds_SNPQC/OvSNP600QCped.ped", InFormat="ped")

AgePriors60QC <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K96,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
Ped60.175QC <- GetMaybeRel(GenoM = infile60QC,
                           LifeHistData = LH60K96,
                           AgePrior = AgePriors60QC,
                           Module = 'ped',
                           Err = 0.175)

output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.175QC$MaybeRel, file = paste0(output_dir, "Ped60.175QC.csv"), row.names = FALSE)

#save 

#############################
#15%(10%)

#LH Input
LH60K15.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.10individuals\\LH60_15.10_4.16.24.csv")

#Ped Input file
infile60K15.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.10individuals\\60QC_filtered15.10.ped", InFormat="ped")

#age prior
AgePriors60.15.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K15.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.15.10 <- GetMaybeRel(GenoM = infile60K15.10,
                           LifeHistData = LH60K15.10,
                           AgePrior = AgePriors60.15.10,
                           Module = 'ped',
                           Err = 0.175)
#write csv
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.15.10$MaybeRel, file = paste0(output_dir, "Ped60.15.10.csv"), row.names = FALSE)


##############################

#############################
#15%(15%)

LH60K15.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.15individuals\\LH60_15.15_4.16.24.csv")

#Ped Input file
infile60K15.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.15individuals\\60QC_filtered15.15.ped", InFormat="ped")

#age prior
AgePriors60.15.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K15.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.15.15 <- GetMaybeRel(GenoM = infile60K15.15,
                           LifeHistData = LH60K15.15,
                           AgePrior = AgePriors60.15.15,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.15.15$MaybeRel, file = paste0(output_dir, "Ped60.15.15.csv"), row.names = FALSE)


##############################

#############################
#15%(20%)

LH60K15.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.20individuals\\LH60_15.20_4.16.24.csv")

#Ped Input file
infile60K15.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\15percentloci\\15.20individuals\\60QC_filtered15.20.ped", InFormat="ped")

#age prior
AgePriors60.15.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K15.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.15.20 <- GetMaybeRel(GenoM = infile60K15.20,
                           LifeHistData = LH60K15.20,
                           AgePrior = AgePriors60.15.20,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.15.20$MaybeRel, file = paste0(output_dir, "Ped60.15.20.csv"), row.names = FALSE)


##############################

#############################
#20%(10%)

LH60K20.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.10individuals\\LH60_20.10_4.16.24.csv")

#Ped Input file
infile60K20.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.10individuals\\60QC_filtered20.10.ped", InFormat="ped")

#age prior
AgePriors60.20.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K20.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.20.10 <- GetMaybeRel(GenoM = infile60K20.10,
                           LifeHistData = LH60K20.10,
                           AgePrior = AgePriors60.20.10,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.20.10$MaybeRel, file = paste0(output_dir, "Ped60.20.10.csv"), row.names = FALSE)


##############################
#############################
#20%(15%)

LH60K20.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.15individuals\\LH60_20.15_4.16.24.csv")

#Ped Input file
infile60K20.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.15individuals\\60QC_filtered20.15.ped", InFormat="ped")

#age prior
AgePriors60.20.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K20.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.20.15 <- GetMaybeRel(GenoM = infile60K20.15,
                           LifeHistData = LH60K20.15,
                           AgePrior = AgePriors60.20.15,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.20.15$MaybeRel, file = paste0(output_dir, "Ped60.20.15.csv"), row.names = FALSE)


##############################

#############################
#20%(20%)

LH60K20.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.20individuals\\LH60_20.20_4.16.24.csv")

#Ped Input file
infile60K20.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\20percentloci\\20.20individuals\\60QC_filtered20.20.ped", InFormat="ped")

#age prior
AgePriors60.20.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K20.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.20.20 <- GetMaybeRel(GenoM = infile60K20.20,
                           LifeHistData = LH60K20.20,
                           AgePrior = AgePriors60.20.20,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.20.20$MaybeRel, file = paste0(output_dir, "Ped60.20.20.csv"), row.names = FALSE)


##############################

#############################
#25%(10%)

LH60K25.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.10individuals\\LH60_25.10_4.16.24.csv")

#Ped Input file
infile60K25.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.10individuals\\60QC_filtered25.10.ped", InFormat="ped")

#age prior
AgePriors60.25.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K25.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.25.10 <- GetMaybeRel(GenoM = infile60K25.10,
                           LifeHistData = LH60K25.10,
                           AgePrior = AgePriors60.25.10,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.25.10$MaybeRel, file = paste0(output_dir, "Ped60.25.10.csv"), row.names = FALSE)


##############################

#############################
#25%(15%)

LH60K25.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.15individuals\\LH60_25.15_4.16.24.csv")

#Ped Input file
infile60K25.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.15individuals\\60QC_filtered25.15.ped", InFormat="ped")

#age prior
AgePriors60.25.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K25.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.25.15 <- GetMaybeRel(GenoM = infile60K25.15,
                           LifeHistData = LH60K25.15,
                           AgePrior = AgePriors60.25.15,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.25.15$MaybeRel, file = paste0(output_dir, "Ped60.25.15.csv"), row.names = FALSE)


##############################

#############################
#25%(20%)

LH60K25.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.20individuals\\LH60_25.20_4.16.24.csv")

#Ped Input file
infile60K25.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\60K\\25percentloci\\25.20individuals\\60QC_filtered25.20.ped", InFormat="ped")

#age prior
AgePriors60.25.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH60K25.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped60.25.20 <- GetMaybeRel(GenoM = infile60K25.20,
                           LifeHistData = LH60K25.20,
                           AgePrior = AgePriors60.25.20,
                           Module = 'ped',
                           Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.25.20$MaybeRel, file = paste0(output_dir, "Ped60.25.20.csv"), row.names = FALSE)


##############################



##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################


#600K



#############################
LH600K15.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.10individuals\\LH600_15.10_4.16.24.csv")

#Ped Input file
infile600K15.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.10individuals\\600QC_filtered15.10.ped", InFormat="ped")

#age prior
AgePriors600.15.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K15.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.15.10 <- GetMaybeRel(GenoM = infile600K15.10,
                            LifeHistData = LH600K15.10,
                            AgePrior = AgePriors600.15.10,
                            Module = 'ped',
                            Err = 0.175)
#write csv
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.15.10$MaybeRel, file = paste0(output_dir, "Ped600.15.10.csv"), row.names = FALSE)


##############################

#############################
#15%(15%)

LH600K15.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.15individuals\\LH600_15.15_4.16.24.csv")

#Ped Input file
infile600K15.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.15individuals\\600QC_filtered15.15.ped", InFormat="ped")

#age prior
AgePriors600.15.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K15.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.15.15 <- GetMaybeRel(GenoM = infile600K15.15,
                            LifeHistData = LH600K15.15,
                            AgePrior = AgePriors600.15.15,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.15.15$MaybeRel, file = paste0(output_dir, "Ped600.15.15.csv"), row.names = FALSE)


##############################

#############################
#15%(20%)

LH600K15.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.20individuals\\LH600_15.20_4.16.24.csv")

#Ped Input file
infile600K15.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\15percentloci\\15.20individuals\\600QC_filtered15.20.ped", InFormat="ped")

#age prior
AgePriors600.15.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K15.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.15.20 <- GetMaybeRel(GenoM = infile600K15.20,
                            LifeHistData = LH600K15.20,
                            AgePrior = AgePriors600.15.20,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.15.20$MaybeRel, file = paste0(output_dir, "Ped600.15.20.csv"), row.names = FALSE)


##############################

#############################
#20%(10%)

LH600K20.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.10individuals\\LH600_20.10_4.16.24.csv")

#Ped Input file
infile600K20.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.10individuals\\600QC_filtered20.10.ped", InFormat="ped")

#age prior
AgePriors600.20.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K20.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.20.10 <- GetMaybeRel(GenoM = infile600K20.10,
                            LifeHistData = LH600K20.10,
                            AgePrior = AgePriors600.20.10,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.20.10$MaybeRel, file = paste0(output_dir, "Ped600.20.10.csv"), row.names = FALSE)

##############################
#############################
#20%(15%)

LH600K20.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.15individuals\\LH600_20.15_4.16.24.csv")

#Ped Input file
infile600K20.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.15individuals\\600QC_filtered20.15.ped", InFormat="ped")

#age prior
AgePriors600.20.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K20.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.20.15 <- GetMaybeRel(GenoM = infile600K20.15,
                            LifeHistData = LH600K20.15,
                            AgePrior = AgePriors600.20.15,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.20.15$MaybeRel, file = paste0(output_dir, "Ped600.20.15.csv"), row.names = FALSE)

##############################

#############################
#20%(20%)

LH600K20.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.20individuals\\LH600_20.20_4.16.24.csv")

#Ped Input file
infile600K20.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\20percentloci\\20.20individuals\\600QC_filtered20.20.ped", InFormat="ped")

#age prior
AgePriors600.20.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K20.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.20.20 <- GetMaybeRel(GenoM = infile600K20.20,
                            LifeHistData = LH600K20.20,
                            AgePrior = AgePriors600.20.20,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.20.20$MaybeRel, file = paste0(output_dir, "Ped600.20.20.csv"), row.names = FALSE)

##############################

#############################
#25%(10%)

LH600K25.10 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.10individuals\\LH600_25.10_4.16.24.csv")

#Ped Input file
infile600K25.10 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.10individuals\\600QC_filtered25.10.ped", InFormat="ped")

#age prior
AgePriors600.25.10 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K25.10,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.25.10 <- GetMaybeRel(GenoM = infile600K25.10,
                            LifeHistData = LH600K25.10,
                            AgePrior = AgePriors600.25.10,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.25.10$MaybeRel, file = paste0(output_dir, "Ped600.25.10.csv"), row.names = FALSE)


##############################

#############################
#25%(15%)

LH600K25.15 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.15individuals\\LH600_25.15_4.16.24.csv")

#Ped Input file
infile600K25.15 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.15individuals\\600QC_filtered25.15.ped", InFormat="ped")

#age prior
AgePriors600.25.15 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K25.15,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.25.15 <- GetMaybeRel(GenoM = infile600K25.15,
                            LifeHistData = LH600K25.15,
                            AgePrior = AgePriors600.25.15,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.25.15$MaybeRel, file = paste0(output_dir, "Ped600.25.15.csv"), row.names = FALSE)


##############################

#############################
#25%(20%)

LH600K25.20 <- read.csv("C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.20individuals\\LH600_25.20_4.16.24.csv")

#Ped Input file
infile600K25.20 <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\SNP_data\\output\\QCFiltering\\PLINK_Sequioa_filtering.process\\600K\\25percentloci\\25.20individuals\\600QC_filtered25.20.ped", InFormat="ped")

#age prior
AgePriors600.25.20 <-MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K25.20,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
#getmayberel PO pairs
Ped600.25.20 <- GetMaybeRel(GenoM = infile600K25.20,
                            LifeHistData = LH600K25.20,
                            AgePrior = AgePriors600.25.20,
                            Module = 'ped',
                            Err = 0.175)
#write csv
# Define the file path
output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.25.20$MaybeRel, file = paste0(output_dir, "Ped600.25.20.csv"), row.names = FALSE)


##############################

#Main split of the three vcf's
##############################

#lowthresholds
LH600K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\Low.Threshold600.csv")
infile600LT <- GenoConvert(InFile = "C:/Users/alec1/Desktop/SNP_data/data/LowThresholds_SNPQC/OvSNP600QCped.ped", InFormat="ped")

AgePriors600LT <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K96,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

Ped600.175LT <- GetMaybeRel(GenoM = infile600LT,
                            LifeHistData = LH600K96,
                            AgePrior = AgePriors600LT,
                            Module = 'ped',
                            Err = 0.175)

output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.175LT$MaybeRel, file = paste0(output_dir, "Ped600.175LT.csv"), row.names = FALSE)


##########################

#Low SNP Quality Control

#600K
LH600K96 <- read.csv("C:\\Users\\alec1\\Desktop\\Updated.Sequioa\\LowThreshold\\LH600_RealMaster.csv")
infile600QC <- GenoConvert(InFile = "C:/Users/alec1/Desktop/SNP_data/data/LowThresholds_SNPQC/OvSNP600QCped.ped", InFormat="ped")

AgePriors600QC <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = LH600K96,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)
Ped600.175QC <- GetMaybeRel(GenoM = infile600QC,
                            LifeHistData = LH600K96,
                            AgePrior = AgePriors600QC,
                            Module = 'ped',
                            Err = 0.175)

output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.175QC$MaybeRel, file = paste0(output_dir, "Ped600.175QC.csv"), row.names = FALSE)

###################################################################



#original
#60k
kn_LH60k <- read.csv("C:/Users/alec1/Desktop/SNP_data/scripts/LH60_OG.csv")
infile60k <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\PLINKSNP\\OvSNP60ped.ped", InFormat="ped")

kn_LH600k <- read.csv("C:/Users/alec1/Desktop/SNP_data/scripts/LH600_OG.csv")
infile600k <- GenoConvert(InFile = "C:\\Users\\alec1\\Desktop\\PLINKSNP\\OvSNP600ped.ped", InFormat="ped")

AgePriors60K <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = kn_LH60k,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

AgePriors600K <- MakeAgePrior(
  Pedigree = NULL,
  LifeHistData = kn_LH600k,
  MinAgeParent = 1,
  MaxAgeParent = NULL,
  Discrete = NULL,
  Flatten = NULL,
  lambdaNW = -log(0.5)/100,
  Smooth = TRUE,
  Plot = TRUE,
  Return = "LR",
  quiet = FALSE
)

Ped60.175K <- GetMaybeRel(GenoM = kn_LH60k,
                          LifeHistData = LH60K96,
                          AgePrior = AgePriors60K,
                          Module = 'ped',
                          Err = 0.175)

Ped600.175K <- GetMaybeRel(GenoM = kn_LH600k,
                           LifeHistData = LH600K96,
                           AgePrior = AgePriors600K,
                           Module = 'ped',
                           Err = 0.175)

output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped60.175K$MaybeRel, file = paste0(output_dir, "Ped60.175K.csv"), row.names = FALSE)



output_dir <- "C:\\Users\\alec1\\Desktop\\CWDProject\\Objective.One_60Kvs.600K_Relatedness\\Scripts\\Final_Scripts\\List.for.Google.Doc\\Ped_LT.OG.QC_Filtered_Both.panels"

write.csv(Ped600.175K$MaybeRel, file = paste0(output_dir, "Ped600.175K.csv"), row.names = FALSE)





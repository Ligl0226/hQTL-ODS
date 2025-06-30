##==================================================================##
##--Description: QTL mapping using powerful hQTL-ODS method
##--Highlight: example for MPH effects, BPH effects, and phenotype itself 
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2024.10.14
##--Date: 2025.06.30 update:
##        (1) data integrity check to ensure that all parental lines referenced 
##            in hybrids are present in both phenotype and genotype files.
##        (2) scan the individual heterotic effect component of the detected hQTL.
##==================================================================##

source("./R/CalKinMat.R")
source("./R/MPH_hQTL_ODS.R")
source("./R/hQTL_component_effects.R")
source("./R/BPH_hQTL_ODS.R")
source("./R/QTL_ODS.R")

WorkDir <- "./"
path_data <- "./data/HyW_PHdata_GNdata_Demo.Rdata"
Memo <- "HyW_Demo"


##==================================================================##
##---Initialization settings
##==================================================================##
cat("\n##=======Initialization settings=======##\n")
setwd(WorkDir);getwd();## dir()
library(gaston)
library(stringr)
library(openxlsx)
library(data.table)


Sys.time()
##==================================================================##
##---loading the phenotype and genotype data
##==================================================================##
load(file = path_data,verbose = T)

## Genotype data
dim(GenoData);GenoData[1:5,1:5]

## Phenotype data
colnames(PhenoData)[1] <- "Genotype"
rownames(PhenoData) <- PhenoData$Genotype
head(PhenoData);tail(PhenoData);dim(PhenoData)


##==================================================================##
## data integrity check to ensure that all parental lines referenced in hybrids 
## are present in both phenotype and genotype files
##==================================================================##
if(all(rownames(GenoData)==PhenoData$Genotype)){
  print("The genotype order in GenoData and PhenoData is perfectly matched.")
}else{
  stop("Please double-check that the genotype order in GenoData and PhenoData is consistent.")
}

Linesindex <- grep("_x_",PhenoData$Genotype,value = FALSE,invert = TRUE)
Hybrdindex <- grep("_x_",PhenoData$Genotype,value = FALSE,invert = FALSE)

Lineslist <- PhenoData$Genotype[Linesindex]
Hybrdlist <- PhenoData$Genotype[Hybrdindex]

Femalelist <- str_split_fixed(Hybrdlist,"_x_",2)[,1]
Malelist <- str_split_fixed(Hybrdlist,"_x_",2)[,2]

if(all(Femalelist==PhenoData$Female[Hybrdindex]) &
   all(Malelist==PhenoData$Male[Hybrdindex])){
  print("All hybrid names meet the required format.")
}else{
  stop("Please double-check that hybrid names are formatted as female_x_male.")
}

FemaleMalelist <- unique(c(Femalelist,Malelist))
tempaa <- setdiff(FemaleMalelist,Lineslist)
if(length(tempaa)>0){
  stop(paste("The ",paste0(tempaa, collapse = ",")," are not listed in data !!!", sep = ""))
}

tempbb <- setdiff(Lineslist,FemaleMalelist)
if(length(tempbb)>0){
  stop(paste("The ",paste0(tempbb, collapse = ",")," were not involved in any hybrid combinations !!! !!!", sep = ""))
}


Sys.time()
#==================================================================##
##---calculating the marker-derived kinship matrices
##==================================================================##
cat("\n##===Calculat the marker-derived kinship matrices!===##\n")
KinMatlist <- CalKinMat(GenoData)
cat("KinMatlist$K_a[1:5,1:5]:\n");KinMatlist$K_a[1:5,1:5]
cat("KinMatlist$K_d[1:5,1:5]:\n");KinMatlist$K_d[1:5,1:5]
cat("KinMatlist$K_aa[1:5,1:5]:\n");KinMatlist$K_aa[1:5,1:5]
cat("KinMatlist$K_ad[1:5,1:5]:\n");KinMatlist$K_ad[1:5,1:5]
cat("KinMatlist$K_dd[1:5,1:5]:\n");KinMatlist$K_dd[1:5,1:5]


#==================================================================##
##---QTL mapping for MPH effect
##   WARNING: The genome-wide scan is time-consuming despite some acceleration algorithms having been applied
##   better distributing the task into parallel sessions on a server
##=================================================================##
MPH_hQTLres <- MPH_hQTL_ODS(PHdata = PhenoData,
                            GNdata = GenoData, 
                            KinMatlist = KinMatlist,
                            P3D = TRUE)
head(MPH_hQTLres);tail(MPH_hQTLres);dim(MPH_hQTLres)


#==================================================================##
##---test individual heterotic effect component of detected hQTL 
##=================================================================##
hQTL_Comp_d <- hQTL_CompEff(PHdata = PhenoData,
                            GNdata = GenoData,
                            KinMatlist = KinMatlist,
                            CompEff = "d")
head(hQTL_Comp_d);tail(hQTL_Comp_d);dim(hQTL_Comp_d)

hSNPslist <- sample(colnames(GenoData),5,replace = FALSE) ## just as a example
hQTL_Comp_aa <- hQTL_CompEff(PHdata = PhenoData,
                             GNdata = GenoData,
                             KinMatlist = KinMatlist,
                             hSNPslist = hSNPslist,
                             CompEff = "aa")
head(hQTL_Comp_aa);tail(hQTL_Comp_aa);dim(hQTL_Comp_aa)

hQTL_Comp_ad <- hQTL_CompEff(PHdata = PhenoData,
                             GNdata = GenoData,
                             KinMatlist = KinMatlist,
                             hSNPslist = hSNPslist,
                             CompEff = "ad")
head(hQTL_Comp_ad);tail(hQTL_Comp_ad);dim(hQTL_Comp_ad)

hQTL_Comp_da <- hQTL_CompEff(PHdata = PhenoData,
                             GNdata = GenoData,
                             KinMatlist = KinMatlist,
                             hSNPslist = hSNPslist,
                             CompEff = "da")
head(hQTL_Comp_da);tail(hQTL_Comp_da);dim(hQTL_Comp_da)

hQTL_Comp_dd <- hQTL_CompEff(PHdata = PhenoData,
                             GNdata = GenoData,
                             KinMatlist = KinMatlist,
                             hSNPslist = hSNPslist,
                             CompEff = "dd")
head(hQTL_Comp_dd);tail(hQTL_Comp_dd);dim(hQTL_Comp_dd)


#==================================================================##
##---QTL mapping for BPH effect
##=================================================================##
BPH_hQTLres <- BPH_hQTL_ODS(PHdata = PhenoData,
                            GNdata = GenoData, 
                            KinMatlist = KinMatlist,
                            P3D = TRUE)
head(BPH_hQTLres);tail(BPH_hQTLres);dim(BPH_hQTLres)


#==================================================================##
##---QTL mapping for trait itself: a+aa effect
##=================================================================##
QTLres <- QTL_ODS(PHdata = PhenoData,
                  GNdata = GenoData, 
                  KinMatlist = KinMatlist,
                  P3D = TRUE,
                  EffectsBeingTested = "a+aa")
head(QTLres);tail(QTLres);dim(QTLres)


#==================================================================##
##---QTL mapping for trait itself: aa+ad+da+dd effect
##=================================================================##
QTLres <- QTL_ODS(PHdata = PhenoData,
                  GNdata = GenoData, 
                  KinMatlist = KinMatlist,
                  P3D = TRUE,
                  EffectsBeingTested = "aa+ad+da+dd")
head(QTLres);tail(QTLres);dim(QTLres)





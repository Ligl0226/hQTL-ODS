##==================================================================##
##--Description: QTL mapping using powerful ODS method
##--Highlight: example for MPH effects, BPH effects, and phenotype itself 
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2024.10.14
##==================================================================##

source("./R/CalKinMat.R")
source("./R/MPH_hQTL_ODS.R")
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

colnames(PhenoData)[1] <- "Genotype"
rownames(PhenoData) <- PhenoData$Genotype

head(PhenoData);tail(PhenoData);dim(PhenoData)
dim(GenoData);GenoData[1:5,1:5]


Sys.time()
#==================================================================##
##---calculating the marker-derived kinship matrices
##==================================================================##
cat("\n##===Calculat the marker-derived kinship matrices!===##\n")
KinMatlist <- CalKinMat(GenoData)
cat("KinMatlist$K_a[1:5,1:5]:\n");KinMatlist$K_a[1:5,1:5]
cat("KinMatlist$K_d[1:5,1:5]:\n");KinMatlist$K_d[1:5,1:5]
cat("KinMatlist$K_ad[1:5,1:5]:\n");KinMatlist$K_ad[1:5,1:5]
cat("KinMatlist$K_da[1:5,1:5]:\n");KinMatlist$K_da[1:5,1:5]
cat("KinMatlist$K_dd[1:5,1:5]:\n");KinMatlist$K_dd[1:5,1:5]


#==================================================================##
##---QTL mapping for MPH effect
##   WARNING: The genome-wide scan is time-consuming despite some acceleration algorithms has been applied
##   better distributing the task into parellel sessions on a server
##=================================================================##
MPH_hQTLres <- MPH_hQTL_ODS(PHdata = PhenoData,
                            GNdata = GenoData, 
                            KinMatlist = KinMatlist,
                            P3D = TRUE)
head(MPH_hQTLres);tail(MPH_hQTLres);dim(MPH_hQTLres)


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





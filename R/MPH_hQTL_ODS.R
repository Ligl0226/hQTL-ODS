##==================================================================##
##--R function: Powerful one-dimensional scan to detect heterotic QTL
##--Highlight: MPH-hQTL-ODS
##--Highlight: update the PVE corrected by VC*(mean(diag(K))-mean(K))
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2023.11.10
##==================================================================##

## source("CalKinMat.R")
library(gaston)

##==================================================================##
##---transform the LRT (LR) to Pvalue based on the theory:
##---If the null hypothesis is true, LR will follow a mixture of 
##   Chi-square 0 and Chi-square 1 distributions with an equal weight.
##==================================================================##
LRtoPvalue <- function(lrt){
  # if (lrt<=1e-03){
  #   lrt.p<-0.5+pchisq(lrt,df=1,lower.tail=FALSE)*0.5
  # }else{
  #   lrt.p<-1-(0.5+pchisq(lrt,df=1)*0.5)
  # }
  lrt.p <- pchisq(lrt,df=1,lower.tail=FALSE)/2
  return(lrt.p)
}


MPH_hQTL_ODS <- function(PHdata, GNdata, KinMatlist=NULL, P3D=TRUE){
  ## PHdata: should be formulated as a data frame with 4 columns
  ##         the first column "Genotype", where the name of hybrids have to be composed of its two corresponding parents name linked by "_x_"
  ##         the second column "Female"
  ##         the third column "Male"
  ##         the fourth column "GrainYield" (target phenotype)
  
  ## GNdata should be formulated as a matrix with 
  ##        rows: genotypes including both parents and hybirds
  ##        colums: markers
  ##        the markers should be coded as 0,1,2 depending the number of minor alleles
  
  TraitName <- colnames(PHdata)[4]
  print(paste("##===The target trait is: ",TraitName,"===##",sep=""))
  
  cat("\n##===Check the individuals in phenotype and genotype data===##\n")
  if (all.equal(rownames(PHdata),rownames(GNdata))) {
    print("The individuals are perfectly matched in phenotype and genotype data.")
  }else{
    stop("Please check the individual order in phenotype and genotype data !!!")
  }
  
  if(is.null(KinMatlist)){
    print("\n##===Calculate the kinship matrix.")
    KinMatlist <- CalKinMat(GNdata)
  }
  
  Sys.time()
  ##==================================================================##
  ##--Generate Adata (-1,0,1) and Ddata (0,1)
  ##==================================================================##
  cat("\n##===Adata (-1,0,1) and Ddata (0,1)!===##\n")
  GenoData <- as.matrix(GNdata)
  Adata <- GenoData-1;
  # cat("dim(Adata):\n");dim(Adata)
  # cat("Adata[1:5,1:5]:\n");Adata[1:5,1:5]
  Ddata <- matrix(1,nrow(Adata),ncol(Adata))-abs(Adata)
  # cat("dim(Ddata):\n");dim(Ddata)
  # cat("Ddata[1:5,1:5]:\n");Ddata[1:5,1:5]
  
  
  Sys.time()
  ##==================================================================##
  ##--calculate AA' and DD' in advance
  ##==================================================================##
  cat("\n##===Calculate AA' and DD' in advance!===##\n")
  AtA <- Adata %*% t(Adata)
  # cat("dim(AtA):\n");dim(AtA)
  # cat("AtA[1:5,1:5]:\n");AtA[1:5,1:5]
  DtD <- Ddata %*% t(Ddata)
  # cat("dim(DtD):\n");dim(DtD)
  # cat("DtD[1:5,1:5]:\n");DtD[1:5,1:5]
  
  
  Sys.time()
  ##==================================================================##
  ##--Generate the transition matrix (Tmat) from trait to MPH
  ##==================================================================##
  cat("\n##===Generate the transition matrix from trait to MPH===##\n")
  ngen <- nrow(PHdata)
  hybrids.index <- grep("_x_",PHdata$Genotype)
  cat("length(hybrids.index):");print(length(hybrids.index))
  names_hyb <- PHdata$Genotype[hybrids.index]
  nhybrid <- length(names_hyb)
  
  Tmat <- matrix(0,nhybrid,ngen)
  rownames(Tmat) <- names_hyb
  colnames(Tmat) <- PHdata$Genotype
  # cat("dim(Tmat):\n");dim(Tmat)
  # cat("Tmat[1:5,1:5]:\n");Tmat[1:5,1:5];
  
  for(i in 1:nhybrid){
    ##i=10
    hybridi <- names_hyb[i]
    femalei <- str_split_fixed(names_hyb[i],"_x_",2)[,1]
    malei <- str_split_fixed(names_hyb[i],"_x_",2)[,2]
    
    Tmat[i,femalei] <- -0.5
    Tmat[i,malei] <- -0.5
    Tmat[i,hybridi] <- 1
  }
  ## cat("dim(Tmat):\n");dim(Tmat)
  ## cat("Tmat[1:5,1:5]:\n");Tmat[1:5,1:5]
  
  cat("Show Tmat randomly:\n")
  for(ra in 1:10){
    ra1=sample(names_hyb,1)
    ra1f=str_split_fixed(ra1,"_x_",2)[,1]
    ra1m=str_split_fixed(ra1,"_x_",2)[,2]
    print(Tmat[ra1,c(ra1f,ra1m,ra1)])
  }
  all(Tmat %*% matrix(rep(1, ngen), ngen, 1)==0)  ##TI==0    TRUE
  all(Tmat %*% as.matrix(Adata[,1])==0)           ##Tmiai==0 TRUE
  all(Tmat %*% as.matrix(Adata[,c(1:5)])==0)      ##TMa==0   TRUE
  
  
  Sys.time()
  ##==================================================================##
  ##--Taking eigen-decomposition of the matrix TT'
  ##==================================================================##
  cat("\n##===Taking eigen-decomposition of the matrix TT'===##\n")
  TTp <- tcrossprod(Tmat)
  eigTTp <- eigen(TTp)
  Dval <- eigTTp$values
  U <- eigTTp$vectors
  V <- U%*%diag(sqrt(Dval))
  Vinv <- solve(V)
  
  VinvTmat <- Vinv %*% Tmat
  tTmattVinv <- t(Tmat) %*% t(Vinv)
  
  # cat("dim(VinvTmat):\n");dim(VinvTmat)
  # cat("VinvTmat[1:10,1:5]:\n");VinvTmat[1:10,1:5]
  # 
  # cat("dim(tTmattVinv):\n");dim(tTmattVinv)
  # cat("tTmattVinv[1:10,1:5]:\n");tTmattVinv[1:10,1:5]
  
  
  Sys.time()
  ##==================================================================##
  ##---common input variables for null model and full model 
  ##==================================================================##
  cat("\n##===Common input variables for null model and full model!===##\n")
  ## before equivalent transformation
  y=PHdata[,TraitName]
  ## cat("summary of",TraitName,"with length =",length(y),":\n");summary(y)
  
  ## after equivalent transformation (variables with wave)
  y=matrix(y,length(y),1)
  MPH <- Tmat %*% y
  ## cat("Summary MPH of",TraitName,"with length =",length(MPH),":\n");summary(MPH[,1])
  y_wave <- VinvTmat %*% y;
  
  Kd_wave <- Vinv %*% Tmat %*% KinMatlist$K_d %*% t(Tmat) %*% t(Vinv)
  Kaa_wave <- Vinv %*% Tmat %*% KinMatlist$K_aa %*% t(Tmat) %*% t(Vinv)
  Kad_wave <- Vinv %*% Tmat %*% KinMatlist$K_ad %*% t(Tmat) %*% t(Vinv)
  Kdd_wave <- Vinv %*% Tmat %*% KinMatlist$K_dd %*% t(Tmat) %*% t(Vinv)
  
  #--Is the following operation necessary?? Yes!!
  Kd_wave2 <- Kd_wave/mean(diag(Kd_wave))
  Kaa_wave2 <- Kaa_wave/mean(diag(Kaa_wave))
  Kad_wave2 <- Kad_wave/mean(diag(Kad_wave))
  Kdd_wave2 <- Kdd_wave/mean(diag(Kdd_wave))
  
  
  Sys.time()
  ##==================================================================##
  ##---solve the null model by gaston
  ##==================================================================##
  cat("\n##===Fitting the null model using lmm.aireml() in gaston package!===##\n")
  gaston_soln0 <- lmm.aireml(Y = y_wave,X = NULL,
                             K = list(Kd_wave2,
                                      Kaa_wave2,
                                      Kad_wave2,
                                      Kdd_wave2),
                             verbose = F,
                             max_iter = 500)
  # str(gaston_soln0)
  # round(c(gaston_soln0$tau,gaston_soln0$sigma2)/sum(gaston_soln0$tau,gaston_soln0$sigma2),4)
  # cat("gaston_soln0$logL:   ");gaston_soln0$logL
  
  
  Sys.time()
  ##==================================================================##
  ##---P3D (population parameters previously determined) methods!
  ##---Estimate the variance components only once in this "null model", 
  ##---then fixed then throughout the testing procedure. 
  ##==================================================================##
  if(P3D==TRUE){
    cat("\n##===P3D (population parameters previously) methods!===##\n")
    lamdapar <- gaston_soln0$tau/gaston_soln0$sigma2
    # cat("The lamda estimated in the null model:\n");print(round(lamdapar,10))
    
    Kc_wave <- lamdapar[1]*Kd_wave2 + lamdapar[2]*Kaa_wave2 + lamdapar[3]*Kad_wave2 + lamdapar[4]*Kdd_wave2
    Kc_wave2 <- Kc_wave/mean(diag(Kc_wave))
    # cat("dim(Kc_wave2):\n");print(dim(Kc_wave2))
    # cat("Kc_wave2[1:5,1:5]:\n");print(Kc_wave2[1:5,1:5])
  }
  
  
  Sys.time()
  ##==================================================================##
  ##--log-likelihood value evaluated from the full model 
  ##--(null model + Heterotic Effect)=full model
  ##==================================================================##
  cat("\n##===Fitting the full model and get the log-likelihood value!===##\n")
  LRres <- as.data.frame(matrix(data = NA,nrow = ncol(Adata),ncol = 7))
  colnames(LRres) <- c("No","Marker","logLik_null","logLik_full","LR","pval","PVE")
  LRres$No <- c(1:ncol(Adata))
  LRres$Marker <- colnames(Adata)
  LRres$logLik_null <- gaston_soln0$logL
  
  Sys.time()
  ## loop for each SNPs
  for(i in 1:ncol(Adata)){
    ##i=1
    
    ## The Adata and Ddata for marker i
    Ai <- as.matrix(Adata[,i])
    Di <- matrix(1-abs(Ai))
    
    ##--Generate and calculate the MM' of the heterotic effect of marker i
    DitDi <- tcrossprod(Di)
    AiA <- tcrossprod(Ai)*AtA-tcrossprod(Ai*Ai)
    AiD <- tcrossprod(Ai)*DtD-tcrossprod(Ai*Di)
    DiA <- tcrossprod(Di)*AtA-tcrossprod(Di*Ai)
    DiD <- tcrossprod(Di)*DtD-tcrossprod(Di*Di)
    
    MtM <- DitDi+AiA+AiD+DiA+DiD
    Khi_wave <- VinvTmat %*% MtM %*% tTmattVinv  ### Khi_wave <- Vinv %*% Tmat %*% tcrossprod(Mmati) %*% t(Tmat) %*% t(Vinv)
    Khi_wave2 <- Khi_wave/mean(diag(Khi_wave))
    
    ##--estimate the log-likelihood value of full model by gaston
    # if(P3D==FALSE){
    #   gaston_solni <- lmm.aireml(Y = y_wave,X = NULL,
    #                              K = list(Khi_wave2,
    #                                       Kd_wave2,
    #                                       Kaa_wave2,
    #                                       Kad_wave2,
    #                                       Kdd_wave2),
    #                              verbose = F)
    # }else{
    #   gaston_solni <- lmm.aireml(Y = y_wave,X = NULL,
    #                              K = list(Khi_wave2,
    #                                       Kc_wave2),
    #                              verbose = F)
    # }
    # VarCompi <- c(gaston_solni$tau,gaston_solni$sigma2)
    
    ## update here for PVE only
    gaston_solni <- NULL
    if(P3D==FALSE){
      gaston_solni <- lmm.aireml(Y = y_wave,X = NULL,
                                 K = list(Khi_wave2,
                                          Kd_wave2,
                                          Kaa_wave2,
                                          Kad_wave2,
                                          Kdd_wave2),
                                 verbose = F)
      VarCompi <- c(gaston_solni$tau,gaston_solni$sigma2)
      VarCompi[1] <- VarCompi[1] * (mean(diag(Khi_wave2))-mean(Khi_wave2))
      VarCompi[2] <- VarCompi[2] * (mean(diag(Kd_wave2))-mean(Kd_wave2))
      VarCompi[3] <- VarCompi[3] * (mean(diag(Kaa_wave2))-mean(Kaa_wave2))
      VarCompi[4] <- VarCompi[4] * (mean(diag(Kad_wave2))-mean(Kad_wave2))
      VarCompi[5] <- VarCompi[5] * (mean(diag(Kdd_wave2))-mean(Kdd_wave2))
    }else{
      gaston_solni <- lmm.aireml(Y = y_wave,X = NULL,
                                 K = list(Khi_wave2,
                                          Kc_wave2),
                                 verbose = F)
      VarCompi <- c(gaston_solni$tau,gaston_solni$sigma2)
      VarCompi[1] <- VarCompi[1] * (mean(diag(Khi_wave2))-mean(Khi_wave2))
      VarCompi[2] <- VarCompi[2] * (mean(diag(Kc_wave2))-mean(Kc_wave2))
    }
    
    VarCompiProp <- VarCompi/sum(VarCompi)
    PVEi <- VarCompiProp[1]
    LRi <- -2*(gaston_soln0$logL-gaston_solni$logL)
	  pvali <- LRtoPvalue(LRi)
    
    ##--return the results
    LRres[i,"logLik_full"] <- gaston_solni$logL
    LRres[i,"LR"] <- LRi
	  LRres[i,"pval"] <- pvali
	  LRres[i,"PVE"] <- PVEi
    
    ##--point print the process
    # cat(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),   ## date()
    #     sprintf("...Done for SNPs %s", i),"\t",
    #     gaston_solni$logL,"\t",gaston_soln0$logL,"\t",LRi,"\t",PVEi,"\n",
    #     file = outlogfile,append = T)
    
    if(i %% 100 == 0){
      print(paste("##-----",i, "/", ncol(Adata),"-----##", sep = " "))
    }
    
  } ## i
  
  return(LRres)
  
}


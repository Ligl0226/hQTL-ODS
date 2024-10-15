##==================================================================##
##--R function: Powerful one-dimensional scan to detect 
##              the cumulative QTL effect of a loci, for example:
##--Highlight: (1) additive + cumulative additive-by-additive  (a+aa)
##--Highlight: (2) cumulative additive-by-additive + 
##                 cumulative additive-by-dominance +
##                 cumulative dominance-by-additive +
##                 cumulative dominance-by-dominance    (aa+ad+da+dd)
##--Highlight: (3) or other combinations based on your own research objects
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2024.10.14
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


QTL_ODS <- function(PHdata, GNdata, KinMatlist=NULL, P3D=TRUE, 
                    EffectsBeingTested="aa+ad+da+dd"){
  ## PHdata: should be formulated as a data frame with 4 columns
  ##         the first column "Genotype"
  ##         the second column "Female"
  ##         the third column "Male"
  ##         the fourth column "GrainYield" (target phenotype)
  
  ## GNdata should be formulated as a matrix with 
  ##        rows: genotypes including both parents and hybrids
  ##        colums: markers
  ##        the markers should be coded as 0,1,2 depending the number of minor alleles
  
  ## EffectsBeingTested could be "a+aa" or "aa+ad+da+dd"
  
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
  
  KinA <- KinMatlist$K_a
  KinD <- KinMatlist$K_d
  KinAA <- KinMatlist$K_aa
  KinAD <- KinMatlist$K_ad
  KinDD <- KinMatlist$K_dd
  
  Sys.time()
  ##==================================================================##
  ##--Generate Adata (-1,0,1) and Ddata (0,1)
  ##==================================================================##
  cat("\n##===Adata (-1,0,1) and Ddata (0,1)!===##\n")
  GenoData <- as.matrix(GNdata)
  Adata <- GenoData-1;
  Ddata <- matrix(1,nrow(Adata),ncol(Adata))-abs(Adata)

  
  Sys.time()
  ##==================================================================##
  ##--calculate AA' and DD' in advance
  ##==================================================================##
  cat("\n##===Calculate AA' and DD' in advance!===##\n")
  AtA <- Adata %*% t(Adata)
  DtD <- Ddata %*% t(Ddata)

  
  Sys.time()
  ##==================================================================##
  ##---common input variables for null model and full model 
  ##==================================================================##
  cat("\n##===Common input variables for null model and full model!===##\n")
  y=PHdata[,TraitName]
  y=matrix(y,length(y),1)

  
  Sys.time()
  ##==================================================================##
  ##---solve the null model by gaston
  ##==================================================================##
  cat("\n##===Fitting the null model using lmm.aireml() in gaston package!===##\n")
  if(EffectsBeingTested=="a+aa"){
    ## null model
    X1 <- matrix(1,length(y),1)
    gaston_soln0 <- lmm.aireml(Y = y,X = X1,K = list(KinA,KinAA),verbose = F,max_iter = 500)
    
    if(P3D==TRUE){
      cat("\n##===P3D (population parameters previously) methods!===##\n")
      lamdapar <- gaston_soln0$tau/gaston_soln0$sigma2
      KinMix <- lamdapar[1]*KinA + lamdapar[2]*KinAA
      KinMix <- KinMix/mean(diag(KinMix))
    }
    
  }else if(EffectsBeingTested=="aa+ad+da+dd"){
    ## null model
    X1 <- matrix(1,length(y),1)
    gaston_soln0 <- lmm.aireml(Y = y,X = X1,K = list(KinA,KinD,KinAA,KinAD,KinDD),verbose = F,max_iter = 500)
    
    if(P3D==TRUE){
      cat("\n##===P3D (population parameters previously) methods!===##\n")
      lamdapar <- gaston_soln0$tau/gaston_soln0$sigma2
      KinMix <- lamdapar[1]*KinA + lamdapar[2]*KinD + lamdapar[3]*KinAA + lamdapar[4]*KinAD + lamdapar[5]*KinDD
      KinMix <- KinMix/mean(diag(KinMix))
    }
  }
  # str(gaston_soln0)
  # round(c(gaston_soln0$tau,gaston_soln0$sigma2)/sum(gaston_soln0$tau,gaston_soln0$sigma2),4)
  # cat("gaston_soln0$logL:   ");gaston_soln0$logL


  Sys.time()
  ##==================================================================##
  ##--log-likelihood value evaluated from the full model 
  ##--(null model + cumulative effects of a loci)=full model
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
    
    if(EffectsBeingTested=="a+aa"){

      ## full model
      Ai <- as.matrix(Adata[,i])
      AitAi <- tcrossprod(Ai)
      AiA <- tcrossprod(Ai)*AtA-tcrossprod(Ai*Ai)
      MtM <- AitAi+AiA
      MtM <- MtM/mean(diag(MtM))

      if(P3D==FALSE){
        gaston_soln1 <- lmm.aireml(Y = y,X = X1,K = list(MtM,KinA,KinAA),verbose = F)
      }else{
        gaston_soln1 <- lmm.aireml(Y = y,X = X1,K = list(MtM,KinMix),verbose = F)
      }
      
      ## likelihood ratio test
      VarCompi <- c(gaston_soln1$tau,gaston_soln1$sigma2)
      VarCompiProp <- VarCompi/sum(VarCompi)
      PVEi <- VarCompiProp[1]
      LRi <- -2*(gaston_soln0$logL-gaston_soln1$logL)
      pvali <- LRtoPvalue(LRi)
      
      LRres[i,"logLik_full"] <- gaston_soln1$logL
      LRres[i,"LR"] <- LRi
      LRres[i,"pval"] <- pvali
      LRres[i,"PVE"] <- PVEi
      
      
    }else if(EffectsBeingTested=="aa+ad+da+dd"){

      ## full model
      Ai <- as.matrix(Adata[,i])
      Di <- matrix(1-abs(Ai))
      AitAi <- tcrossprod(Ai)
      DitDi <- tcrossprod(Di)
      AiA <- tcrossprod(Ai)*AtA-tcrossprod(Ai*Ai)
      AiD <- tcrossprod(Ai)*DtD-tcrossprod(Ai*Di)
      DiA <- tcrossprod(Di)*AtA-tcrossprod(Di*Ai)
      DiD <- tcrossprod(Di)*DtD-tcrossprod(Di*Di)
      MtM <- AiA+AiD+DiA+DiD
      MtM <- MtM/mean(diag(MtM))
      
      if(P3D==FALSE){
        gaston_soln1 <- lmm.aireml(Y = y,X = X1,K = list(MtM,KinA,KinD,KinAA,KinAD,KinDD),verbose = F)
      }else{
        gaston_soln1 <- lmm.aireml(Y = y,X = X1,K = list(MtM,KinMix),verbose = F)
      }
      
      ## likelihood ratio test
      VarCompi <- c(gaston_soln1$tau,gaston_soln1$sigma2)
      VarCompiProp <- VarCompi/sum(VarCompi)
      PVEi <- VarCompiProp[1]
      LRi <- -2*(gaston_soln0$logL-gaston_soln1$logL)
      pvali <- LRtoPvalue(LRi)
      
      LRres[i,"logLik_full"] <- gaston_soln1$logL
      LRres[i,"LR"] <- LRi
      LRres[i,"pval"] <- pvali
      LRres[i,"PVE"] <- PVEi
    }
    
    # ## The Adata and Ddata for marker i
    # Ai <- as.matrix(Adata[,i])
    # Di <- matrix(1-abs(Ai))
    # 
    # ##--Generate and calculate the MM' of the effect of marker i being tested
    # AitAi <- tcrossprod(Ai)
    # DitDi <- tcrossprod(Di)
    # AiA <- tcrossprod(Ai)*AtA-tcrossprod(Ai*Ai)
    # AiD <- tcrossprod(Ai)*DtD-tcrossprod(Ai*Di)
    # DiA <- tcrossprod(Di)*AtA-tcrossprod(Di*Ai)
    # DiD <- tcrossprod(Di)*DtD-tcrossprod(Di*Di)
    # 
    # MtM <- AitAi+DitDi+AiA+AiD+DiA+DiD
    # Khi_wave <- VinvTmat %*% MtM %*% tTmattVinv  ### Khi_wave <- Vinv %*% Tmat %*% tcrossprod(Mmati) %*% t(Tmat) %*% t(Vinv)
    # Khi_wave2 <- Khi_wave/mean(diag(Khi_wave))
    
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


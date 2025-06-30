##==================================================================##
##--R function: to investigate whether these component effects of hQTL are significant or not
##--Highlight: test individual heterotic effect component of detected hQTL 
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2025.06.30
##==================================================================##

# source("CalKinMat.R")
# library(stringr)
# library(gaston)

hQTL_CompEff <- function(PHdata, GNdata, KinMatlist=NULL, hSNPslist, CompEff="d"){
  ## PHdata: should be formulated as a data frame with 4 columns
  ##         the first column "Genotype", where the name of hybrids have to be composed of its two corresponding parents name linked by "_x_"
  ##         the second column "Female"
  ##         the third column "Male"
  ##         the fourth column "GrainYield" (target phenotype)
  
  ## GNdata should be formulated as a matrix with 
  ##        rows: genotypes including both parents and hybirds
  ##        colums: markers
  ##        the markers should be coded as 0,1,2 depending the number of minor alleles
  
  ## hSNPslist: the heterotic SNPs detected from hQTL-ODS scanning;
  
  ## CompEff: must be one of the following effects: d, aa, ad, da, or dd.
  
  
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
  P3D=TRUE
  if(P3D==TRUE){
    cat("\n##===P3D (population parameters previously) methods!===##\n")
    lamdapar <- gaston_soln0$tau/gaston_soln0$sigma2
    # cat("The lamda estimated in the null model:\n");print(round(lamdapar,10))
    
    Kc_wave <- lamdapar[1]*Kd_wave2 + lamdapar[2]*Kaa_wave2 + 
      lamdapar[3]*Kad_wave2 + lamdapar[4]*Kdd_wave2
    Kc_wave2 <- Kc_wave/mean(diag(Kc_wave))
    # cat("dim(Kc_wave2):\n");print(dim(Kc_wave2))
    # cat("Kc_wave2[1:5,1:5]:\n");print(Kc_wave2[1:5,1:5])
  }
  
  
  ##==================================================================##
  ## Using the combined kinship matrix controlling the genetic background
  ## transform the LMM to LM
  ##==================================================================##
  dim(Kc_wave2);Kc_wave2[1:5,1:5]
  eigK <- eigen(Kc_wave2)
  d <- eigK$values
  W <- eigK$vectors
  # y_star  <- as.vector(t(W)%*%y_wave)
  wt <- 1/(d+1)
  TRAN <- diag(sqrt(wt))%*%t(W)%*%VinvTmat
  cat("dim(TRAN):\n");print(dim(TRAN))
  
  Y <- TRAN%*%y
  cat("dim(Y):\n");print(dim(Y))
  
  trAdata <- TRAN%*%Adata
  # trAdata <- prod_XY_Cpp(TRAN,Adata)
  cat("dim(trAdata):\n");print(dim(trAdata))
  
  trDdata <- TRAN%*%Ddata
  # trDdata <- prod_XY_Cpp(TRAN,Ddata)
  cat("dim(trDdata):\n");print(dim(trDdata))
  
  
  ##==================================================================##
  ## Genome-wide scan for d/aa/ad/da/dd
  ##==================================================================##
  nmar <- ncol(Adata)
  total_ma_names <- colnames(Adata)
  
  if(CompEff=="d"){
    
    cat("\n##===Scanning for component effect: dominance effect for all markers.===##\n")
    DomRes <- data.frame(Marker=total_ma_names,Effect=NA,Pvalue=NA,PVE=NA)
    # head(DomRes);tail(DomRes);dim(DomRes)
    
    for (m in 1:nmar){
      subfinalDATA <- data.frame(y=Y,trDdata[,m])
      colnames(subfinalDATA)[2] <- total_ma_names[m]
      
      fit.lm <- NULL
      summytable <- NULL
      anovatable <- NULL
      effval <- NULL
      pvalue <- NULL
      pve <- NULL
      
      fit.lm <- lm(as.formula(paste("y ~",total_ma_names[m],"-1")),data=subfinalDATA)
      summytable <- summary(fit.lm)$coefficient
      anovatable <- anova(fit.lm)
      
      if(nrow(summytable)>=1 & nrow(anovatable)>=2){
        effval <- summary(fit.lm)$coefficient[1,1]
        pvalue <- summary(fit.lm)$coefficient[1,4]
        pve <- anovatable[total_ma_names[m],"Sum Sq"]/sum(anovatable[,"Sum Sq"])
      }else{
        effval <- 0
        pvalue <- 1
        pve <- 0
      }
      DomRes[m,"Effect"] <- effval
      DomRes[m,"Pvalue"] <- pvalue
      DomRes[m,"PVE"] <- pve

      if (m %% 1000 == 0){
        cat(date(),"Dominance_Effect_Scanning: Marker",m,"completed\n")
      }
    }
    return(DomRes)
    
    
  }else if(CompEff=="aa"){
    
    cat("\n##===Scanning for component effect: additive-by-additive effect between hSNPs and whole genomic background.===##\n")
    tempaa <- setdiff(hSNPslist,total_ma_names)
    if(length(tempaa)>0){
      stop(paste("The ",paste0(tempaa,collapse = ",")," are not in GenoData !!!",sep = ""))
    }
    ResAA <- data.frame(hSNP=rep(hSNPslist,each=nmar),
                        Marker=rep(total_ma_names,times=length(hSNPslist)),
                        Effect=NA,Pvalue=NA,PVE=NA)
    
    for (m in 1:length(hSNPslist)){
      for (n in 1:nmar){
        epidesign <- TRAN%*%(Adata[,m]*Adata[,n])
        subdata <- data.frame(y=Y,trAdata[,m],trAdata[,n],epi=epidesign)
        colnames(subdata)[2:3] <- c(hSNPslist[m],total_ma_names[n])
        
        fit.lm <- NULL
        infomat <- NULL
        infoaov <- NULL
        effval <- NULL
        pvalue <- NULL
        pve <- NULL
        fit.lm <- lm(as.formula(paste("y ~ -1 +",hSNPslist[m],"+",
                                      total_ma_names[n],"+ epi")),
                     data=subdata)
        infomat <- summary(fit.lm)$coefficient
        infoaov <- anova(fit.lm)

        if ("epi"%in%rownames(infomat)){
          effval <- summary(fit.lm)$coefficient["epi",1]
          pvalue <- summary(fit.lm)$coefficient["epi",4]
          pve <- infoaov["epi","Sum Sq"]/sum(infoaov[,"Sum Sq"])
        }else{
          effval <- 0
          pvalue <- 1
          pve <- 0
        }
        ResAA[(nmar*(m-1)+n),"Effect"] <- effval
        ResAA[(nmar*(m-1)+n),"Pvalue"] <- pvalue
        ResAA[(nmar*(m-1)+n),"PVE"] <- pve
      } 
      cat(date(),"additive-by-additive scan invloving hSNPs",m,hSNPslist[m],"completed.\n") 
    }
    return(ResAA)
    
    
  }else if(CompEff=="ad"){
    
    cat("\n##===Scanning for component effect: additive-by-dominance effect between hSNPs and whole genomic background.===##\n")
    tempad <- setdiff(hSNPslist,total_ma_names)
    if(length(tempad)>0){
      stop(paste("The ",paste0(tempad,collapse = ",")," are not in GenoData !!!",sep = ""))
    }
    ResAD <- data.frame(hSNP=rep(hSNPslist,each=nmar),
                        Marker=rep(total_ma_names,times=length(hSNPslist)),
                        Effect=NA,Pvalue=NA,PVE=NA)
    
    for (m in 1:length(hSNPslist)){
      for (n in 1:nmar){
        epidesign <- TRAN%*%(Adata[,m]*Ddata[,n])
        subdata <- data.frame(y=Y,trAdata[,m],trDdata[,n],epi=epidesign)
        colnames(subdata)[2:3] <- c(hSNPslist[m],total_ma_names[n])
        
        fit.lm <- NULL
        infomat <- NULL
        infoaov <- NULL
        effval <- NULL
        pvalue <- NULL
        pve <- NULL
        fit.lm <- lm(as.formula(paste("y ~ -1 +",hSNPslist[m],"+",
                                      total_ma_names[n],"+ epi")),
                     data=subdata)
        infomat <- summary(fit.lm)$coefficient
        infoaov <- anova(fit.lm)

        if ("epi"%in%rownames(infomat)){
          effval <- summary(fit.lm)$coefficient["epi",1]
          pvalue <- summary(fit.lm)$coefficient["epi",4]
          pve <- infoaov["epi","Sum Sq"]/sum(infoaov[,"Sum Sq"])
        }else{
          effval <- 0
          pvalue <- 1
          pve <- 0
        }
        ResAD[(nmar*(m-1)+n),"Effect"] <- effval
        ResAD[(nmar*(m-1)+n),"Pvalue"] <- pvalue
        ResAD[(nmar*(m-1)+n),"PVE"] <- pve
      } 
      cat(date(),"additive-by-dominance scan invloving hSNPs",m,hSNPslist[m],"completed.\n") 
    }
    return(ResAD)
    
    
  }else if(CompEff=="da"){
    
    cat("\n##===Scanning for component effect: dominance-by-additive effect between hSNPs and whole genomic background.===##\n")
    tempda <- setdiff(hSNPslist,total_ma_names)
    if(length(tempda)>0){
      stop(paste("The ",paste0(tempda,collapse = ",")," are not in GenoData !!!",sep = ""))
    }
    ResDA <- data.frame(hSNP=rep(hSNPslist,each=nmar),
                        Marker=rep(total_ma_names,times=length(hSNPslist)),
                        Effect=NA,Pvalue=NA,PVE=NA)
    
    for (m in 1:length(hSNPslist)){
      for (n in 1:nmar){
        epidesign <- TRAN%*%(Ddata[,m]*Adata[,n])
        subdata <- data.frame(y=Y,trDdata[,m],trAdata[,n],epi=epidesign)
        colnames(subdata)[2:3] <- c(hSNPslist[m],total_ma_names[n])
        
        fit.lm <- NULL
        infomat <- NULL
        infoaov <- NULL
        effval <- NULL
        pvalue <- NULL
        pve <- NULL
        fit.lm <- lm(as.formula(paste("y ~ -1 +",hSNPslist[m],"+",
                                      total_ma_names[n],"+ epi")),
                     data=subdata)
        infomat <- summary(fit.lm)$coefficient
        infoaov <- anova(fit.lm)
        
        if ("epi"%in%rownames(infomat)){
          effval <- summary(fit.lm)$coefficient["epi",1]
          pvalue <- summary(fit.lm)$coefficient["epi",4]
          pve <- infoaov["epi","Sum Sq"]/sum(infoaov[,"Sum Sq"])
        }else{
          effval <- 0
          pvalue <- 1
          pve <- 0
        }
        ResDA[(nmar*(m-1)+n),"Effect"] <- effval
        ResDA[(nmar*(m-1)+n),"Pvalue"] <- pvalue
        ResDA[(nmar*(m-1)+n),"PVE"] <- pve
      } 
      cat(date(),"dominance-by-additive scan invloving hSNPs",m,hSNPslist[m],"completed.\n") 
    }
    return(ResDA)
    
    
  }else if(CompEff=="dd"){
    
    cat("\n##===Scanning for component effect: dominance-by-dominance effect between hSNPs and whole genomic background.===##\n")
    tempdd <- setdiff(hSNPslist,total_ma_names)
    if(length(tempdd)>0){
      stop(paste("The ",paste0(tempdd,collapse = ",")," are not in GenoData !!!",sep = ""))
    }
    ResDD <- data.frame(hSNP=rep(hSNPslist,each=nmar),
                        Marker=rep(total_ma_names,times=length(hSNPslist)),
                        Effect=NA,Pvalue=NA,PVE=NA)
    
    for (m in 1:length(hSNPslist)){
      for (n in 1:nmar){
        epidesign <- TRAN%*%(Ddata[,m]*Ddata[,n])
        subdata <- data.frame(y=Y,trDdata[,m],trDdata[,n],epi=epidesign)
        colnames(subdata)[2:3] <- c(hSNPslist[m],total_ma_names[n])
        
        fit.lm <- NULL
        infomat <- NULL
        infoaov <- NULL
        effval <- NULL
        pvalue <- NULL
        pve <- NULL
        fit.lm <- lm(as.formula(paste("y ~ -1 +",hSNPslist[m],"+",
                                      total_ma_names[n],"+ epi")),
                     data=subdata)
        infomat <- summary(fit.lm)$coefficient
        infoaov <- anova(fit.lm)
        
        if ("epi"%in%rownames(infomat)){
          effval <- summary(fit.lm)$coefficient["epi",1]
          pvalue <- summary(fit.lm)$coefficient["epi",4]
          pve <- infoaov["epi","Sum Sq"]/sum(infoaov[,"Sum Sq"])
        }else{
          effval <- 0
          pvalue <- 1
          pve <- 0
        }
        ResDD[(nmar*(m-1)+n),"Effect"] <- effval
        ResDD[(nmar*(m-1)+n),"Pvalue"] <- pvalue
        ResDD[(nmar*(m-1)+n),"PVE"] <- pve
      } 
      cat(date(),"dominance-by-dominance scan invloving hSNPs",m,hSNPslist[m],"completed.\n") 
    }
    return(ResDD)

    
  }else{
    stop("CompEff must be one of the following effects: d, aa, ad, da, or dd.")
  }
  
} ## function





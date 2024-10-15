##==================================================================##
##--R function: Calculate the kinship matrix based on the genotype data
##--Highlight: coding as 0,1,2
##--Maintainer: GuoLiang Li (lig@ipk-gatersleben.de)
##--Date: 2023.11.10
##==================================================================##

CalKinMat <- function(GNdata){
  ## GNdata should be formulated as a matrix with 
  ## rows: genotypes including both parents and hybirds
  ## colums: markers
  ## the markers should be coded as 0,1,2 depending the number of minor alleles
  res <- list()
  
  GenoData <- as.matrix(GNdata)
  alf <- apply(GenoData,2,mean)/2
  W <- t(t(GenoData)-2*alf)
  c <- 2*sum(alf*(1-alf))
  
  K_a <- W %*% t(W)
  K_aa <- 0.5 * (K_a*K_a - (W*W) %*% (t(W)*t(W)))

  Ddata <- GenoData
  x1 <- -2*(1-alf)^2
  x2 <- 2*alf*(1-alf)
  x3 <- -2*alf^2
  for (i in 1:ncol(Ddata)){
    Ddata[which(GenoData[,i]==2),i] <- x1[i]
    Ddata[which(GenoData[,i]==1),i] <- x2[i]
    Ddata[which(GenoData[,i]==0),i] <- x3[i]
    Ddata[which(GenoData[,i]==1.5),i] <- (x1[i]+x2[i])/2
    Ddata[which(GenoData[,i]==0.5),i] <- (x3[i]+x2[i])/2
  }
  
  c2 <- sum(4*(alf*(1-alf))^2)
  W2 <- Ddata 
  
  K_d <- W2 %*% t(W2)
  K_dd <- 0.5 * (K_d*K_d - (W2*W2) %*% (t(W2)*t(W2)))
  K_ad <- K_a*K_d-(W*W2)%*%(t(W)*t(W2))

  res$K_a <- K_a/c
  res$K_aa <- K_aa/c^2
  res$K_d <- K_d/c2
  res$K_dd <- K_dd/c2^2
  res$K_ad <- K_ad/(c*c2)
  
  return(res)
}


#########################################################
##### SIMPLE LINEAR REGRESSION: "CORRECTED" DATASETS ####
##### SELECTION FREQUENCIES + PREDICTIVE PERFORMANCE ####
#########################################################

library(BayesVarSel)
library(HDInterval)
library(bestglm)


a <- Ozone35[,2:8]
a$c <- rep(1,178)

# Variance functions for ozone dataset
varhetx7 <- exp(1+0.052*a$x7)
varhetx7_10 <- 4.831*178/sum(varhetx7)*varhetx7
varhetx7_1 <- 48.31*178/sum(varhetx7)*varhetx7

varhetx6 <- exp(1+0.038*a$x6)
varhetx6_10 <- 4.831*178/sum(varhetx6)*varhetx6
varhetx6_1 <- 48.31*178/sum(varhetx6)*varhetx6



ozonefreqchol <- function(sdev, beta0, beta1, nsim){
  
  hpm <- matrix(nrow=nsim,ncol=8)
  mpm <- matrix(data=0,nrow=nsim,ncol=8)
  aic <- matrix(nrow=nsim,ncol=8)
  
  for (i in 1:nsim) {
    print(i)
    a$y <- beta0 + beta1 * a$x7 + rnorm(178, sd = sdev)
    
    ozonehomo <- a/sdev
    
    mod <- Bvs(data=ozonehomo,y~.-1,null.model=y~c-1)
    hpm[i,] <- c(mod$HPMbin,all(mod$HPMbin==c(0,0,0,1,0,0,0)))
    
    for (j in 1:7) {
      if (mod$inclprob[j] > 0.5) {
        mpm[i, j] <- 1
      }
    }
    
    mpm[i,8] <- all(mpm[i,]==c(0,0,0,1,0,0,0,0))

    bglmaic <- bestglm(ozonehomo,IC="AIC",RequireFullEnumerationQ=TRUE,intercept=FALSE)
    modaic <- unlist(bglmaic$BestModels[1,1:7])
    aic[i,] <- c(modaic,all(modaic==c(0,0,0,1,0,0,0)))
  }
  
  final <- data.frame(matrix(ncol=8,nrow=0))
  colnames(final) <- c("x4","x5","x6","x7","x8","x9","x10","mod certo")
  
  final[nrow(final)+1,] <- apply(hpm,2,mean)
  final[nrow(final)+1,] <- apply(mpm,2,mean)
  final[nrow(final)+1,] <- apply(aic,2,mean)

  return(final)
}

# SNR = 10

set.seed(4653)
freqocholhomo10 <- ozonefreqchol(sqrt(4.831),2,0.5,5000)
set.seed(56265)
freqocholhetx7_10 <- ozonefreqchol(sqrt(varhetx7_10),2,0.5,5000)
set.seed(9153)
freqocholhetx6_10  <- ozonefreqchol(sqrt(varhetx61),2,0.5,5000) 

# SNR = 1

set.seed(3500)
freqohomo1 <- ozonefreqchol(sqrt(48.31),2,0.5,5000) 
set.seed(43589)
freqocholhetx6_1 <- ozonefreqchol(sqrt(varhetx62),2,0.5,5000)
set.seed(16789)
freqocholhetx7_1 <- ozonefreqchol(sqrt(varhetx7_1),2,0.5,5000)

######################################################

# predictive performance

predozonechol <- function(sdev, beta0, beta1,nsim,nkeep){
  
  mse <- rep(NA,nsim)
  hdi <- rep(NA,nsim)
  coverage <- rep(NA,nsim)
  
  for (i in 1:nsim) {
    print(i)
    a$y <<- beta0 + beta1 * a$x7 + rnorm(178, sd = sdev)
    
    indextrain <- sample(seq_len(178), size = 160)
    
    ozonehomo <- a/sdev
    
    train <- ozonehomo[indextrain,]
    test <- ozonehomo[-indextrain,]
    
    mod <- Bvs(data=train,y~.-1,null.model=y~c-1,n.keep=nkeep)
    
    pred <- predict(mod, newdata=test[,1:8], n.sim=10000) 
    
    meanpred <- apply(pred,2,function(x) mean(x)) 
    
    msei <- sum((test$y - meanpred) ^ 2) / sum((test$y - mean(test$y))^2)
    
    mse[i] <- msei
    auxlength <- 0
    auxcoverage <- 0
    for(j in 1:ncol(pred)){
      dpred <- density(pred[,j])
      upper <- unname(hdi(dpred)[2])
      lower <- unname(hdi(dpred)[1])
      auxlength <- auxlength + (upper - lower)/(max(pred[,j])-min(pred[,j]))
      isinint <- ifelse(upper >= test$y[j] & lower <= test$y[j],1,0)
      auxcoverage <- auxcoverage+isinint/nrow(test)
    }
    hdi[i] <- auxlength / nrow(test)
    coverage[i] <- auxcoverage
  }
  return(c(mean(mse),mean(hdi),mean(coverage)))
}

# BMA

# SNR = 10

set.seed(214)
pred_ozone_homo_10_chol <- predozonechol(sqrt(4.831),2,0.5,10000,128)
set.seed(115)
pred_ozone_het_x7_10_chol <- predozonechol(sqrt(varhetx7_10),2,0.5,10000,128)
set.seed(932)
pred_ozone_het_x6_10_chol <- predozonechol(sqrt(varhetx6_10),2,0.5,10000,128)

# SNR = 1

set.seed(671)
pred_ozone_homo_1_chol <- predozonechol(sqrt(48.31),2,0.5,10000,128)
set.seed(873)
pred_ozone_het_x7_1_chol <- predozonechol(sqrt(varhetx7_1),2,0.5,10000,128)
set.seed(102)
pred_ozone_het_x6_1_chol <- predozonechol(sqrt(varhetx6_1),2,0.5,10000,128)

# HPM

# SNR = 10

set.seed(214)
pred_ozone_homo_10_chol_hpm <- predozonechol(sqrt(4.831),2,0.5,10000,1)
set.seed(115)
pred_ozone_het_x7_10_chol_hpm <- predozonechol(sqrt(varhetx7_10),2,0.5,10000,1)
set.seed(932)
pred_ozone_het_x6_10_chol_hpm <- predozonechol(sqrt(varhetx6_10),2,0.5,10000,1)

# SNR = 1

set.seed(671)
pred_ozone_homo_1_chol_hpm <- predozonechol(sqrt(48.31),2,0.5,10000,1)
set.seed(873)
pred_ozone_het_x7_1_chol_hpm <- predozonechol(sqrt(varhetx7_1),2,0.5,10000,1)
set.seed(102)
pred_ozone_het_x6_1_chol_hpm <- predozonechol(sqrt(varhetx6_1),2,0.5,10000,1)
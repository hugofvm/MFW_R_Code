#########################################################
###### SIMPLE LINEAR REGRESSION: ORIGINAL DATASET #######
##### SELECTION FREQUENCIES + PREDICTIVE PERFORMANCE ####
#########################################################

library(BayesVarSel)
library(HDInterval)

# dataset
a <- Ozone35[,2:8]

# variance functions
varhetx7 <- exp(1+0.052*a$x7)
varhetx7_10 <- 4.831*178/sum(varhetx7)*varhetx7
varhetx7_1 <- 48.31*178/sum(varhetx7)*varhetx7

varhetx6 <- exp(1+0.038*a$x6)
varhetx6_10 <- 4.831*178/sum(varhetx6)*varhetx6
varhetx6_1 <- 48.31*178/sum(varhetx6)*varhetx6

######################################################

# selection frequencies

ozonefreq <- function(sdev, beta0, beta1, nsim){
  
  hpm <- matrix(nrow=nsim,ncol=8)
  mpm <- matrix(data=0,nrow=nsim,ncol=8)
  lasso <- matrix(nrow=nsim,ncol=8)
  aic <- matrix(nrow=nsim,ncol=8)
  bic <- matrix(nrow=nsim,ncol=8)
  
  X <- scale(as.matrix(a[,-8]))
  
  for (i in 1:nsim) {
    print(i)
    a$y <<- beta0 + beta1 * a$x7 + rnorm(178, sd = sdev)
    
    mod <- Bvs(data = a, formula = y ~ .)
    hpm[i,] <- c(mod$HPMbin,all(mod$HPMbin==c(0,0,0,1,0,0,0)))
    
    for (j in 1:7) {
      if (mod$inclprob[j] > 0.5) {
        mpm[i, j] <- 1
      }
    }
    
    mpm[i,8] <- all(mpm[i,]==c(0,0,0,1,0,0,0,0)) # extra zero
    
    Y <- scale(as.numeric(a$y))
    modlasso <- cv.glmnet(X, Y)
    where.lambdamin <- which(modlasso$lambda == modlasso$lambda.min)
    gamma <- as.numeric(modlasso$glmnet.fit$beta[, where.lambdamin] != 0)
    
    lasso[i,] <- c(gamma,all(gamma==c(0,0,0,1,0,0,0)))
    
    bglmaic <- bestglm(a,IC="AIC",RequireFullEnumerationQ=TRUE)
    modaic <- unlist(bglmaic$BestModels[1,1:7])
    aic[i,] <- c(modaic,all(modaic==c(0,0,0,1,0,0,0)))
    
    bglmbic <- bestglm(a,IC="BIC",RequireFullEnumerationQ=TRUE)
    modbic <- unlist(bglmbic$BestModels[1,1:7])
    bic[i,] <- c(modbic,all(modbic==c(0,0,0,1,0,0,0)))
  }
  
  final <- data.frame(matrix(ncol=8,nrow=0))
  colnames(final) <- c("x4","x5","x6","x7","x8","x9","x10","mod certo")
  
  final[nrow(final)+1,] <- apply(hpm,2,mean)
  final[nrow(final)+1,] <- apply(mpm,2,mean)
  final[nrow(final)+1,] <- apply(aic,2,mean)
  final[nrow(final)+1,] <- apply(bic,2,mean)
  final[nrow(final)+1,] <- apply(lasso,2,mean)
  
  return(final)
}

# SNR = 10

set.seed(4653)
freqohomo10 <- ozonefreq(sqrt(4.831),2,0.5,5000)
set.seed(56265)
freqohetx7_10 <- ozonefreq(sqrt(varhetx7_10),2,0.5,5000)
set.seed(9153)
freqohetx6_10 <- ozonefreq(sqrt(varhetx6_10),2,0.5,5000)

# SNR = 1

set.seed(3500)
freqohomo1 <- ozonefreq(sqrt(48.31),2,0.5,5000) 
set.seed(16789)
freqohetx7_1 <- ozonefreq(sqrt(varhetx7_1),2,0.5,5000) 
set.seed(4358)
freqohetx6_1 <- ozonefreq(sqrt(varhetx6_1),2,0.5,5000)

######################################################

# predictive performance

predozone <- function(sdev, beta0, beta1, nsim, nkeep){
  
  mse <- rep(NA,nsim)
  hdi <- rep(NA,nsim)
  coverage <- rep(NA,nsim)
  
  for (i in 1:nsim) {
    print(i)
    
    a$y <<- beta0 + beta1 * a$x7 + rnorm(178, sd = sdev)
    
    indextrain <- sample(seq_len(178), size = 160) # 18 for prediction
    
    train <- a[indextrain,]
    test <- a[-indextrain,] 
    
    mod <- Bvs(data = train, formula = y ~ ., n.keep=nkeep) 
    
    pred <- predict(mod, newdata=test[,1:7], n.sim=10000)
 
    meanpred <- apply(pred,2,function(x) mean(x))
    
    msei <- sum((test$y - meanpred) ^ 2) / sum((test$y - mean(test$y))^2)

    mse[i] <- msei
    
    auxlength <- 0
    auxcoverage <- 0
    for(j in 1:ncol(pred)){
      dpred <- density(pred[,j])
      
      upper <- unname(hdi(dpred)[2])
      lower <- unname(hdi(dpred)[1])
      
      auxlength = auxlength + (upper - lower)/(max(pred[,j])-min(pred[,j]))
      
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
pred_ozone_homo_10 <- predozone(sqrt(4.831),2,0.5,10000,128)
set.seed(115)
pred_ozone_het_x7_10 <- predozone(sqrt(varhetx7_10),2,0.5,10000,128)
set.seed(932)
pred_ozone_het_x6_10 <- predozone(sqrt(varhetx6_10),2,0.5,10000,128)

# SNR = 1

set.seed(671)
pred_ozone_homo_1 <- predozone(sqrt(48.31),2,0.5,10000,128)
set.seed(873)
pred_ozone_het_x7_1 <- predozone(sqrt(varhetx7_1),2,0.5,10000,128)
set.seed(102)
pred_ozone_het_x6_1 <- predozone(sqrt(varhetx6_1),2,0.5,10000,128)

# HPM

# SNR = 10

set.seed(214)
pred_ozone_homo_10_hpm <- predozone(sqrt(4.831),2,0.5,10000,1)
set.seed(115)
pred_ozone_het_x7_10_hpm <- predozone(sqrt(varhetx7_10),2,0.5,10000,1)
set.seed(932)
pred_ozone_het_x6_10_hpm <- predozone(sqrt(varhetx6_10),2,0.5,10000,1)

# SNR = 1

set.seed(671)
pred_ozone_homo_1_hpm <- predozone(sqrt(48.31),2,0.5,10000,1)
set.seed(873)
pred_ozone_het_x7_1_hpm <- predozone(sqrt(varhetx7_1),2,0.5,10000,1)
set.seed(102)
pred_ozone_het_x6_1_hpm <- predozone(sqrt(varhetx6_1),2,0.5,10000,1)
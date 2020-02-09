#########################################################
## MULTIVARIATE LINEAR REGRESSION: "CORRECTED" DATASETS #
##### SELECTION FREQUENCIES + PREDICTIVE PERFORMANCE ####
#########################################################

library(BayesVarSel)
library(HDInterval)

# dataset
cr <- UScrime[,1:15]
cr$c <- rep(1,47)

# Variance functions
# 1 - exponential
varhet1 <- exp(seq(1,3.7,length.out=47))
varhet1_10 <- 312.4117*47/sum(varhet1)*varhet1
varhet1_1 <- 3124.117*47/sum(varhet1)*varhet1

# 2 - half observations with 12.25 times more variance than the other half
varhet2 <- c(rep(3124.117*3.5, times=23),rep(3124.117/3.5, times=24))
varhet2_10 <- 312.4117*47/sum(varhet2)*varhet2
varhet2_1 <- 3124.117*47/sum(varhet2)*varhet2

# 3 - exp(1+0.0186 Pop)
varhet3 <- exp(1+0.0186*cr$Pop)
varhet3_10 <- 312.4117*47/sum(varhet3)*varhet3
varhet3_1 <- 3124.117*47/sum(varhet3)*varhet3

# 4 - exp(1+0.0077 GDP)
varhet4 <- exp(1+0.0077 *cr$GDP)
varhet4_10 <- 312.4117*47/sum(varhet4)*varhet4
varhet4_1 <- 3124.117*47/sum(varhet4)*varhet4

######################################################

# selection frequencies

crimefreqchol <- function(sdev, beta0, beta, nsim){
  
  hpm <- matrix(nrow=nsim,ncol=16)
  mpm <- matrix(data=0,nrow=nsim,ncol=16)
  aic <- matrix(nrow=nsim,ncol=16)
  
  certo <- beta != 0
  
  for (i in 1:nsim) {
    print(i)
    cr$y <- beta0 + as.matrix(cr[,1:15]) %*% beta + rnorm(47, sd = sdev)    
    
    crhomo <- cr/sdev
    
    mod <- Bvs(data=crhomo,y~.-1,null.model=y~c-1)
    hpm[i,] <- c(mod$HPMbin,all(mod$HPMbin==certo))
    
    for (j in 1:15) {
      if (mod$inclprob[j] > 0.5) {
        mpm[i, j] <- 1
      }
    }
    
    mpm[i,16] <- all(mpm[i,]==c(certo,0))
    
    bglmaic <- bestglm(crhomo,IC="AIC",intercept=FALSE)
    modaic <- unlist(bglmaic$BestModels[1,1:15])
    aic[i,] <- c(modaic,all(modaic==certo))
  }
  
  final <- data.frame(matrix(ncol=16,nrow=0))
  colnames(final) <- c("M","So","Ed","Po1","Po2","LF","M.F","Pop","NW","U1","U2","GDP","Ineq","Prob","Time", "mod certo")
  
  final[nrow(final)+1,] <- apply(hpm,2,mean)
  final[nrow(final)+1,] <- apply(mpm,2,mean)
  final[nrow(final)+1,] <- apply(aic,2,mean)
  
  return(final)
}

# SNR = 10

set.seed(293)
freqcrime_chol_het1_10 <- crimefreqchol(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(104)
freqcrime_chol_het2_10 <- crimefreqchol(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(5634)
freqcrime_chol_het3_10 <- crimefreqchol(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(3843)
freqcrime_chol_het4_10 <- crimefreqchol(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)

# SNR = 1

set.seed(810)
freqcrime_chol_het1_1 <- crimefreqchol(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(435)
freqcrime_chol_het2_1 <- crimefreqchol(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(784)
freqcrime_chol_het3_1 <- crimefreqchol(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(183)
freqcrime_chol_het4_1 <- crimefreqchol(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)

######################################################

# predictive performance

predcrimechol <- function(sdev, beta0, beta,nsim,nkeep){
  
  mse <- rep(NA,nsim)
  hdi <- rep(NA,nsim)
  coverage <- rep(NA,nsim)
  
  for (i in 1:nsim) {
    print(i)
    cr$y <<- beta0 + as.matrix(cr[,1:15]) %*% beta + rnorm(47, sd = sdev)
    
    crhomo <- cr/sdev
    
    indextrain <- sample(seq_len(47), size = 40) # 7 for prediction
    
    train <- crhomo[indextrain,]
    test <- crhomo[-indextrain,]
    
    mod <- Bvs(data=train,y~.-1,null.model=y~c-1,n.keep=nkeep)
    
    pred <- predict(mod, newdata=test[,1:16], n.sim=10000)
    meanpred <- apply(pred,2,function(x) mean(x))
    msei <- sum((test$y - meanpred) ^ 2) / sum((test$y - mean(test$y))^2)
    
    mse[i] <- msei
    auxlength <- 0
    auxcoverage <- 0
    for(j in 1:ncol(pred)){
      dpred <- density(pred[,j])
      #dpred <- pred[,j]
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

set.seed(932)
pred_crime_homo_10_chol <- predcrimechol(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),500,256)
set.seed(185)
pred_crime_het1_10_chol <- predcrimechol(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(202)
pred_crime_het2_10_chol <- predcrimechol(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(640)
pred_crime_het3_10_chol <- predcrimechol(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(327)
pred_crime_het4_10_chol <- predcrimechol(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)

# SNR = 1

set.seed(730)
pred_crime_homo_1_chol <- predcrimechol(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),500,256)
set.seed(949)
pred_crime_het1_1_chol <- predcrimechol(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(893)
pred_crime_het2_1_chol <- predcrimechol(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(267)
pred_crime_het3_1_chol <- predcrimechol(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(620)
pred_crime_het4_1_chol <- predcrimechol(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)


# HPM

# SNR = 10

set.seed(932)
pred_crime_homo_10_chol_hpm <- predcrimechol(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(185)
pred_crime_het1_10_chol_hpm <- predcrimechol(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(202)
pred_crime_het2_10_chol_hpm <- predcrimechol(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(640)
pred_crime_het3_10_chol_hpm <- predcrimechol(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(327)
pred_crime_het4_10_chol_hpm <- predcrimechol(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)

# SNR = 1

set.seed(730)
pred_crime_homo_1_chol_hpm <- predcrimechol(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(949)
pred_crime_het1_1_chol_hpm <- predcrimechol(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(893)
pred_crime_het2_1_chol_hpm <- predcrimechol(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(267)
pred_crime_het3_1_chol_hpm <- predcrimechol(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(620)
pred_crime_het4_1_chol_hpm <- predcrimechol(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
#########################################################
### MULTIVARIATE LINEAR REGRESSION: ORIGINAL DATASET ####
##### SELECTION FREQUENCIES + PREDICTIVE PERFORMANCE ####
#########################################################

library(BayesVarSel)
library(HDInterval)

# dataset
cr <- UScrime[,1:15]

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

crimefreq <- function(sdev, beta0, beta, nsim){
  
  hpm <- matrix(nrow=nsim,ncol=16)
  mpm <- matrix(data=0,nrow=nsim,ncol=16)
  lasso <- matrix(nrow=nsim,ncol=16)
  aic <- matrix(nrow=nsim,ncol=16)
  bic <- matrix(nrow=nsim,ncol=16)
  
  X <- as.matrix(cr[,-16])
  certo <- beta != 0
  
  for (i in 1:nsim) {
    print(i)
    cr$y <- beta0 + as.matrix(cr[,1:15]) %*% beta + rnorm(47, sd = sdev)    
    
    mod <- Bvs(data = cr, formula = y ~ .)
    
    hpm[i,] <- c(mod$HPMbin,all(mod$HPMbin==certo))
    
    for (j in 1:15) {
      if (mod$inclprob[j] > 0.5) {
        mpm[i, j] <- 1
      }
    }
    
    mpm[i,16] <- all(mpm[i,]==c(certo,0)) # um zero a mais
    
    Y <- as.numeric(cr$y)
    modlasso <- cv.glmnet(X, Y)
    where.lambdamin <- which(modlasso$lambda == modlasso$lambda.min)
    gamma <- as.numeric(modlasso$glmnet.fit$beta[, where.lambdamin] != 0)
    
    lasso[i,] <- c(gamma,all(gamma==certo))
    
    bglmaic <- bestglm(cr,IC="AIC")
    modaic <- unlist(bglmaic$BestModels[1,1:15])
    aic[i,] <- c(modaic,all(modaic==certo))
    
    bglmbic <- bestglm(cr,IC="BIC")
    modbic <- unlist(bglmbic$BestModels[1,1:15])
    bic[i,] <- c(modbic,all(modbic==certo))
  }
  
  final <- data.frame(matrix(ncol=16,nrow=0))
  colnames(final) <- c("M","So","Ed","Po1","Po2","LF","M.F","Pop","NW","U1","U2","GDP","Ineq","Prob","Time", "mod certo")
  
  final[nrow(final)+1,] <- apply(hpm,2,mean)
  final[nrow(final)+1,] <- apply(mpm,2,mean)
  final[nrow(final)+1,] <- apply(aic,2,mean)
  final[nrow(final)+1,] <- apply(bic,2,mean)
  final[nrow(final)+1,] <- apply(lasso,2,mean)
  
  return(final)
}


# SNR = 10

set.seed(163)
freqcrimehomo10 <- crimefreq(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(1567)
freqcrimehet1_10 <- crimefreq(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(15648)
freqcrimehet2_10 <- crimefreq(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(564)
freqcrimehet3_10 <- crimefreq(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(384)
freqcrimehet4_10 <- crimefreq(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)

# SNR = 1

set.seed(1486)
freqcrimehomo1 <- crimefreq(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(1568)
freqcrimehet1_1 <- crimefreq(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(15612)
freqcrimehet2_1 <- crimefreq(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(375)
freqcrimehet3_1 <- crimefreq(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)
set.seed(476)
freqcrimehet4_1 <- crimefreq(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),5000)

######################################################

# predictive performance
predcrime <- function(sdev, beta0, beta, nsim,nkeep){
  
  mse <- rep(NA,nsim)
  hdi <- rep(NA,nsim)
  coverage <- rep(NA,nsim)
  
  for (i in 1:nsim) {
    print(i)
    cr$y <<- beta0 + as.matrix(cr[,1:15]) %*% beta + rnorm(47, sd = sdev)
    
    indextrain <- sample(seq_len(47), size = 40)
    
    train <- cr[indextrain,]
    test <- cr[-indextrain,] 
    
    mod <- Bvs(data = train, formula = y ~ ., n.keep=nkeep)
    
    pred <- predict(mod, newdata=test[,1:15], n.sim=10000)
    
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

set.seed(932)
pred_crime_homo_10 <- predcrime(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(185)
pred_crime_het1_10 <- predcrime(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(202)
pred_crime_het2_10 <- predcrime(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(640)
pred_crime_het3_10 <- predcrime(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(327)
pred_crime_het4_10 <- predcrime(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)

# SNR = 1

set.seed(730)
pred_crime_homo_1 <- predcrime(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(949)
pred_crime_het1_1 <- predcrime(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(893)
pred_crime_het2_1 <- predcrime(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(267)
pred_crime_het3_1 <- predcrime(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)
set.seed(620)
pred_crime_het4_1 <- predcrime(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,256)


# HPM

# SNR = 10

set.seed(932)
pred_crime_homo_10_hpm <- predcrime(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(185)
pred_crime_het1_10_hpm <- predcrime(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(202)
pred_crime_het2_10_hpm <- predcrime(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(640)
pred_crime_het3_10_hpm <- predcrime(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(327)
pred_crime_het4_10_hpm <- predcrime(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)

# SNR = 1

set.seed(730)
pred_crime_homo_1_hpm <- predcrime(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(949)
pred_crime_het1_1_hpm <- predcrime(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(893)
pred_crime_het2_1_hpm <- predcrime(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(267)
pred_crime_het3_1_hpm <- predcrime(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
set.seed(620)
pred_crime_het4_1_hpm <- predcrime(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)

#########################################################
### PRELIMINARY CALCULATIONS: SNR & SD RATIO CHECKS + ###
###            HETEROSCEDASTICITY TESTS               ###
#########################################################

library(BayesVarSel)
library(lmtest)

# SNR & SD Ratio functions

snr <- function(x, beta, Sigma){
  num <- var(as.matrix(x) %*% beta)
  return(num/(sum(Sigma)/nrow(x)))
}

sdratio <- function(var){
  sd <- sqrt(var)
  top <- mean(sd[order(sd,decreasing=TRUE)[1:floor(0.1*length(sd))]])
  bottom <- mean(sd[order(sd,decreasing=FALSE)[1:floor(0.1*length(sd))]])
  return(top/bottom)
}


# datasets
a <- Ozone35[,2:8]
cr <- UScrime[,1:15]

# Variance functions
# for ozone dataset
varhetx7 <- exp(1+0.052*a$x7)
varhetx7_10 <- 4.831*178/sum(varhetx7)*varhetx7
varhetx7_1 <- 48.31*178/sum(varhetx7)*varhetx7

varhetx6 <- exp(1+0.038*a$x6)
varhetx6_10 <- 4.831*178/sum(varhetx6)*varhetx6
varhetx6_1 <- 48.31*178/sum(varhetx6)*varhetx6

# for crime dataset
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

# Find out the variance that leads to SNR = 1
snr(a[,1:7],c(0,0,0,0.5,0,0,0),rep(1,times=178))
snr((cr[,1:15]),c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),rep(1,times=47))

# check if SD Ratios are in compliance (around 3.5)
sdratio(varhetx7)
sdratio(varhetx6)

sdratio(varhet1)
sdratio(varhet2)
sdratio(varhet3)
sdratio(varhet4)

############################
# heteroscedasticity tests #
############################

##### OZONE DATASET #####

hettestozone <- function(sdev, beta0, beta,nsim,chol){
  res <- data.frame(bpval = rep(NA,times=nsim), wpval = rep(NA,times=nsim))
  
  for (i in 1:nsim) {
    print(i)
    a$y <- beta0 + beta * a$x7 + rnorm(178, sd = sdev)
    
    if(chol==1){
      ozone <- a/sdev
    }
    else{ozone <- a}

    mod <- lm(data=ozone, y ~ .)
    
    res[i,1] <- unlist(bptest(mod)$p.value)
    res[i,2] <- unlist(bptest(mod,~mod$fitted.values+I(mod$fitted.values^2))$p.value)
  }
  return(res)
}

set.seed(869)
hettest_ozone_homo_1 <- hettestozone(sqrt(48.31),2,0.5,10000,0)
hettest_ozone_hetx7_1 <- hettestozone(sqrt(varhetx7_1),2,0.5,10000,0)
hettest_ozone_hetx6_1 <- hettestozone(sqrt(varhetx6_1),2,0.5,10000,0)

set.seed(538)
hettest_ozone_homo_10 <- hettestozone(sqrt(4.831),2,0.5,10000,0)
hettest_ozone_hetx7_10 <- hettestozone(sqrt(varhetx7_10),2,0.5,10000,0)
hettest_ozone_hetx6_10 <- hettestozone(sqrt(varhetx6_10),2,0.5,10000,0)

sum(hettest_ozone_homo_10$bpval < 0.05)/nrow(hettest_ozone_homo_10)
sum(hettest_ozone_homo_10$wpval < 0.05)/nrow(hettest_ozone_homo_10)
sum(hettest_ozone_hetx7_10$bpval < 0.05)/nrow(hettest_ozone_hetx7_10)
sum(hettest_ozone_hetx7_10$wpval < 0.05)/nrow(hettest_ozone_hetx7_10)
sum(hettest_ozone_hetx6_10$bpval < 0.05)/nrow(hettest_ozone_hetx6_10)
sum(hettest_ozone_hetx6_10$wpval < 0.05)/nrow(hettest_ozone_hetx6_10)

sum(hettest_ozone_homo_1$bpval < 0.05)/nrow(hettest_ozone_homo_1)
sum(hettest_ozone_homo_1$wpval < 0.05)/nrow(hettest_ozone_homo_1)
sum(hettest_ozone_hetx7_1$bpval < 0.05)/nrow(hettest_ozone_hetx7_1)
sum(hettest_ozone_hetx7_1$wpval < 0.05)/nrow(hettest_ozone_hetx7_1)
sum(hettest_ozone_hetx6_1$bpval < 0.05)/nrow(hettest_ozone_hetx6_1)
sum(hettest_ozone_hetx6_1$wpval < 0.05)/nrow(hettest_ozone_hetx6_1)

# for corrected dataset

set.seed(869)
hettest_ozone_homo_1_chol <- hettestozone(sqrt(48.31),2,0.5,10000,1)
hettest_ozone_hetx7_1_chol <- hettestozone(sqrt(varhetx7_1),2,0.5,10000,1)
hettest_ozone_hetx6_1_chol <- hettestozone(sqrt(varhetx6_1),2,0.5,10000,1)

set.seed(538)
hettest_ozone_homo_10_chol <- hettestozone(sqrt(4.831),2,0.5,10000,1)
hettest_ozone_hetx7_10_chol <- hettestozone(sqrt(varhetx7_10),2,0.5,10000,1)
hettest_ozone_hetx6_10_chol <- hettestozone(sqrt(varhetx6_10),2,0.5,10000,1)

sum(hettest_ozone_homo_10_chol$bpval < 0.05)/nrow(hettest_ozone_homo_10_chol)
sum(hettest_ozone_homo_10_chol$wpval < 0.05)/nrow(hettest_ozone_homo_10_chol)
sum(hettest_ozone_hetx7_10_chol$bpval < 0.05)/nrow(hettest_ozone_hetx7_10_chol)
sum(hettest_ozone_hetx7_10_chol$wpval < 0.05)/nrow(hettest_ozone_hetx7_10_chol)
sum(hettest_ozone_hetx6_10_chol$bpval < 0.05)/nrow(hettest_ozone_hetx6_10_chol)
sum(hettest_ozone_hetx6_10_chol$wpval < 0.05)/nrow(hettest_ozone_hetx6_10_chol)

sum(hettest_ozone_homo_1_chol$bpval < 0.05)/nrow(hettest_ozone_homo_1_chol)
sum(hettest_ozone_homo_1_chol$wpval < 0.05)/nrow(hettest_ozone_homo_1_chol)
sum(hettest_ozone_hetx7_1_chol$bpval < 0.05)/nrow(hettest_ozone_hetx7_1_chol)
sum(hettest_ozone_hetx7_1_chol$wpval < 0.05)/nrow(hettest_ozone_hetx7_1_chol)
sum(hettest_ozone_hetx6_1_chol$bpval < 0.05)/nrow(hettest_ozone_hetx6_1_chol)
sum(hettest_ozone_hetx6_1_chol$wpval < 0.05)/nrow(hettest_ozone_hetx6_1_chol)



##### CRIME DATASET #####

cr <- UScrime[,1:15]

hettestcrime <- function(sdev, beta0, beta,nsim){
  
  res <- data.frame(bpval = rep(NA, times=nsim),wpval = rep(NA, times=nsim))
  
  for (i in 1:nsim) {
    print(i)
    if(chol==1){
      crime <- cr/sdev
    }
    else{crime <- cr}
    
    mod <- lm(data=crime, y ~ .)
    
    res[i,1] <- unlist(bptest(mod)$p.value)
    res[i,2] <- unlist(bptest(mod,~mod$fitted.values+I(mod$fitted.values^2))$p.value)
    
  }
  
  return(res)
}

hettest_crime_homo_1 <- hettestcrime(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_homo_10 <- hettestcrime(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het1_1 <- hettestcrime(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het1_10 <- hettestcrime(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het2_1 <- hettestcrime(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het2_10 <- hettestcrime(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het3_1 <- hettestcrime(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het3_10 <- hettestcrime(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het4_1 <- hettestcrime(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)
hettest_crime_het4_10 <- hettestcrime(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,0)

sum(hettest_crime_homo_10$bpval < 0.05)/nrow(hettest_crime_homo_10)
sum(hettest_crime_homo_10$wpval < 0.05)/nrow(hettest_crime_homo_10)
sum(hettest_crime_het1_10$bpval < 0.05)/nrow(hettest_crime_het1_10)
sum(hettest_crime_het1_10$wpval < 0.05)/nrow(hettest_crime_het1_10)
sum(hettest_crime_het2_10$bpval < 0.05)/nrow(hettest_crime_het2_10)
sum(hettest_crime_het2_10$wpval < 0.05)/nrow(hettest_crime_het2_10)
sum(hettest_crime_het3_10$bpval < 0.05)/nrow(hettest_crime_het3_10)
sum(hettest_crime_het3_10$wpval < 0.05)/nrow(hettest_crime_het3_10)
sum(hettest_crime_het4_10$bpval < 0.05)/nrow(hettest_crime_het4_10)
sum(hettest_crime_het4_10$wpval < 0.05)/nrow(hettest_crime_het4_10)

sum(hettest_crime_homo_1$bpval < 0.05)/nrow(hettest_crime_homo_1)
sum(hettest_crime_homo_1$wpval < 0.05)/nrow(hettest_crime_homo_1)
sum(hettest_crime_het1_1$bpval < 0.05)/nrow(hettest_crime_het1_1)
sum(hettest_crime_het1_1$wpval < 0.05)/nrow(hettest_crime_het1_1)
sum(hettest_crime_het2_1$bpval < 0.05)/nrow(hettest_crime_het2_1)
sum(hettest_crime_het2_1$wpval < 0.05)/nrow(hettest_crime_het2_1) 
sum(hettest_crime_het3_1$bpval < 0.05)/nrow(hettest_crime_het3_1)
sum(hettest_crime_het3_1$wpval < 0.05)/nrow(hettest_crime_het3_1)
sum(hettest_crime_het4_1$bpval < 0.05)/nrow(hettest_crime_het4_1)
sum(hettest_crime_het4_1$wpval < 0.05)/nrow(hettest_crime_het4_1)


# for corrected dataset

hettest_crime_homo_1_chol <- hettestcrime(sqrt(3124.117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_homo_10_chol <- hettestcrime(sqrt(312.4117),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het1_1_chol <- hettestcrime(sqrt(varhet1_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het1_10_chol <- hettestcrime(sqrt(varhet1_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het2_1_chol <- hettestcrime(sqrt(varhet2_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het2_10_chol <- hettestcrime(sqrt(varhet2_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het3_1_chol <- hettestcrime(sqrt(varhet3_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het3_10_chol <- hettestcrime(sqrt(varhet3_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het4_1_chol <- hettestcrime(sqrt(varhet4_1),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)
hettest_crime_het4_10_chol <- hettestcrime(sqrt(varhet4_10),10,c(0,0,0.1,0,0.5,0,0,1.2,0,0,0,0,-0.2,100,-0.1),10000,1)

sum(hettest_crime_homo_10_chol$bpval < 0.05)/nrow(hettest_crime_homo_10_chol)
sum(hettest_crime_homo_10_chol$wpval < 0.05)/nrow(hettest_crime_homo_10_chol)
sum(hettest_crime_het1_10_chol$bpval < 0.05)/nrow(hettest_crime_het1_10_chol)
sum(hettest_crime_het1_10_chol$wpval < 0.05)/nrow(hettest_crime_het1_10_chol)
sum(hettest_crime_het2_10_chol$bpval < 0.05)/nrow(hettest_crime_het2_10_chol)
sum(hettest_crime_het2_10_chol$wpval < 0.05)/nrow(hettest_crime_het2_10_chol)
sum(hettest_crime_het3_10_chol$bpval < 0.05)/nrow(hettest_crime_het3_10_chol)
sum(hettest_crime_het3_10_chol$wpval < 0.05)/nrow(hettest_crime_het3_10_chol)
sum(hettest_crime_het4_10_chol$bpval < 0.05)/nrow(hettest_crime_het4_10_chol)
sum(hettest_crime_het4_10_chol$wpval < 0.05)/nrow(hettest_crime_het4_10_chol)

sum(hettest_crime_homo_1_chol$bpval < 0.05)/nrow(hettest_crime_homo_1_chol)
sum(hettest_crime_homo_1_chol$wpval < 0.05)/nrow(hettest_crime_homo_1_chol)
sum(hettest_crime_het1_1_chol$bpval < 0.05)/nrow(hettest_crime_het1_1_chol)
sum(hettest_crime_het1_1_chol$wpval < 0.05)/nrow(hettest_crime_het1_1_chol)
sum(hettest_crime_het2_1_chol$bpval < 0.05)/nrow(hettest_crime_het2_1_chol)
sum(hettest_crime_het2_1_chol$wpval < 0.05)/nrow(hettest_crime_het2_1_chol) 
sum(hettest_crime_het3_1_chol$bpval < 0.05)/nrow(hettest_crime_het3_1_chol)
sum(hettest_crime_het3_1_chol$wpval < 0.05)/nrow(hettest_crime_het3_1_chol)
sum(hettest_crime_het4_1_chol$bpval < 0.05)/nrow(hettest_crime_het4_1_chol)
sum(hettest_crime_het4_1_chol$wpval < 0.05)/nrow(hettest_crime_het4_1_chol)







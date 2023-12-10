### 100(1-alpha)% confidence intervals for cases ###

m <- 10000  # the number of times to calculate the confidence interval
alpha <- 0.05
n <- c(50,100,500,1000,5000,10000)  # sample size
coverage_prob1 <- numeric(length(n))  # a vector for the simple method
coverage_prob2 <- numeric(length(n))  # a vector for Fisher's z method
set.seed(2023)

# true value of TP,FP,FN,TN
tTP <- 0.45
tFP <- 0.05
tFN <- 0.05
tTN <- 1-tTP-tFP-tFN

# true value of MCC
tMCC <- (tTP*tTN-tFP*tFN)/sqrt((tTP+tFP)*(tTP+tFN)*(tFP+tTN)*(tFN+tTN))
tMCC

# the function that generates the limiting variance of p_hat
variance_matrix <- function(p) {
    mat <- -outer(p, p)
    diag(mat) <- p * (1-p)
    return(mat)
}

for(j in 1:length(n)){
    # data generation
    data <- rmultinom(m, size=n[j], prob = c(tTP,tFP,tFN,tTN))
    TP <- data[1,]/n[j]
    FP <- data[2,]/n[j]
    FN <- data[3,]/n[j]
    TN <- data[4,]/n[j]

    # a column of MCC estimates
    MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(FP+TN)*(FN+TN))

    # calculattion of the differentiation for the function phi
    denom1 <- 2 * ((TP+FP)*(TP+FN)*(FP+TN)*(FN+TN))^(3/2)   #denominator
    dphi_TP <- (FP+TN)*(FN+TN)*(2*TP*FP*FN + 2*FP*FN*TN + TP*FN*TN + TP*FP*TN + FP^2*FN + FP*FN^2)/denom1   # differentiation with respect to TP
    dphi_FP <- -(TP+FN)*(FN+TN)*(2*TP*FN*TN + 2*TP*FP*TN + TP*FP*FN + FP*FN*TN + TP^2*TN + TP*TN^2)/denom1   # differentiation with respect to FP
    dphi_FN <- -(TP+FP)*(FP+TN)*(2*TP*FP*TN + 2*TP*FN*TN + TP*FP*FN + FP*FN*TN + TP^2*TN + TP*TN^2)/denom1   # differentiation with respect to FN
    dphi_TN <- (TP+FP)*(TP+FN)*(2*FP*FN*TN + 2*TP*FP*FN + TP*FN*TN + TP*FP*TN + FP^2*FN + FP*FN^2)/denom1   # differentiation with respect to TN

    # differentiation of the composite function of f and phi
    denom2 <- 2 * sqrt((TP + FP) * (TP + FN) * (FP + TN) * (FN + TN)) * (FP * FN - TP * TN + sqrt((TP + FP) * (TP + FN) * (FP + TN) * (FN+ TN))) * (-FP * FN + TP * TN + sqrt((TP + FP) * (TP + FN) * (FP + TN) * (FN + TN)))   #denominator
    dfphi_TP <- ((FP + TN) * (FN + TN) * (2*TP*FP*FN + 2*FP*FN*TN + TP*FN*TN + TP*FP*TN + FP^2*FN + FP*FN^2)) / denom2   # differentiation with respect to TP
    dfphi_FP <- -((TP + FN) * (FN + TN) * (2*TP*FN*TN + 2*TP*FP*TN + TP*FP*FN + FP*FN*TN + TP^2*TN + TP*TN^2)) / denom2   # differentiation with respect to FP
    dfphi_FN <- -((TP + FP) * (FP + TN) * (2*TP*FP*TN + 2*TP*FN*TN + TP*FP*FN + FP*FN*TN + TP^2*TN + TP*TN^2)) / denom2   # differentiation with respect to FN
    dfphi_TN <- ((TP + FP) * (TP + FN) * (2*FP*FN*TN + 2*TP*FP*FN + TP*FN*TN + TP*FP*TN + FP^2*FN + FP*FN^2)) / denom2   # differentiation with respect to TN

    # the limiting variance of p_hat
    Sigma <-  lapply(1:m, function(i) variance_matrix(data[,i]/n[j]))

    ## the limiting variance of MCC calculated by simple method
    Var_simple <- lapply(1:m, function(i) c(dphi_TP[i], dphi_FP[i], dphi_FN[i], dphi_TN[i]) %*% Sigma[[i]] %*% c(dphi_TP[i], dphi_FP[i], dphi_FN[i],dphi_TN[i])/n[j])
    Var_simple <- unlist(Var_simple)
    
    ## Simple method
    k <- numeric(m)
    interval1 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(Var_simple[i]))
    interval1 <- unlist(interval1)
    k <- ifelse(MCC-interval1<tMCC & tMCC<MCC+interval1, 1, 0)
    k <- k[complete.cases(k)]
    coverage_prob1[j] <- mean(k)

    ## Fisher's z method
    z <- (1/2) * log((1+MCC)/(1-MCC))
    Var_fisher <- lapply(1:m, function(i) c(dfphi_TP[i], dfphi_FP[i], dfphi_FN[i], dfphi_TN[i]) %*% Sigma[[i]] %*% c(dfphi_TP[i], dfphi_FP[i], dfphi_FN[i], dfphi_TN[i]))
    Var_fisher <- unlist(Var_fisher)
    interval2 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(Var_fisher[i]/n[j]))
    interval2 <- unlist(interval2)
    zL <- z - interval2
    zU <- z + interval2
    L <- (exp(2*zL)-1)/((exp(2*zL)+1))
    U <- (exp(2*zU)-1)/((exp(2*zU)+1))
    l <- numeric(m)
    l <- ifelse(L<tMCC & tMCC<U, 1, 0)
    l <- l[complete.cases(l)]
    coverage_prob2[j] <- mean(l)
}
options(digits = 4)
result <- data.frame(row.names=c("n=50","n=100","n=500","n=1000","n=5000","n=10000"), Simple_method=coverage_prob1, Fishers_z_method=coverage_prob2)
print(paste("P(y=1)=", tTP + tFN, "MCC=", tMCC))
result




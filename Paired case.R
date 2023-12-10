### 100(1-alpha)% confidence intervals for paired design ###

m <- 10000  # the number of times to calculate the confidence interval
alpha <- 0.05
n <- c(50,100,500,1000,5000,100000) # sample size
coverage_prob1 <- numeric(length(n))  # a vector to contain the result for the simple method
coverage_prob2 <- numeric(length(n))  # a vector to contain the result for the Zou's method
coverage_prob3 <- numeric(length(n))  # a vector to contain the result for the simple method
k1 <- numeric(m)
k2 <- numeric(m)
k3 <- numeric(m)
set.seed(2023)

# true value of p_111 to p_000
p_111 <- 0.401
p_101 <- 0.049
p_011 <- 0.049
p_001 <- 0.001
p_110 <- 0.01
p_100 <- 0.04
p_010 <- 0.04
p_000 <- 1-p_111-p_101-p_011-p_001-p_110-p_100-p_010

# true value of MCC1
tMCC1 <- ((p_111+p_101)*(p_010+p_000)-(p_110+p_100)*(p_011+p_001))/sqrt(((p_111+p_101)+(p_110+p_100))*((p_111+p_101)+(p_011+p_001))*((p_110+p_100)+(p_010+p_000))*((p_011+p_001)+(p_010+p_000)))
tMCC1
# true value of MCC1
tMCC2 <- ((p_111+p_011)*(p_100+p_000)-(p_110+p_010)*(p_101+p_001))/sqrt(((p_111+p_011)+(p_110+p_010))*((p_111+p_011)+(p_101+p_001))*((p_110+p_010)+(p_100+p_000))*((p_101+p_001)+(p_100+p_000)))
tMCC2

# the function that generates the limiting variance of p_hat
variance_matrix <- function(p) {
    mat <- -outer(p, p)
    diag(mat) <- p * (1-p)
    return(mat)
}

for(j in 1:length(n)){
    data <- rmultinom(m, size=n[j], prob = c(p_111,p_101,p_011,p_001,p_110,p_100,p_010,p_000))
    p <- data/n[j]
    TP1 <- (p[1,]+p[2,])
    FP1 <- (p[5,]+p[6,])
    FN1 <- (p[3,]+p[4,])
    TN1 <- (p[7,]+p[8,])
    TP2 <- (p[1,]+p[3,])
    FP2 <- (p[5,]+p[7,])
    FN2 <- (p[2,]+p[4,])
    TN2 <- (p[6,]+p[8,])
    
    # a column of MCC estimates
    MCC1 <- (TP1*TN1-FP1*FN1)/sqrt((TP1+FP1)*(TP1+FN1)*(FP1+TN1)*(FN1+TN1))
    MCC2 <- (TP2*TN2-FP2*FN2)/sqrt((TP2+FP2)*(TP2+FN2)*(FP2+TN2)*(FN2+TN2))

    # the limiting variance of p_hat
    Sigma <-  lapply(1:m, function(i) variance_matrix(data[,i]/n[j]))

    ## perform substitutions to speed up the calculation
    term001 <- p[1,] + p[2,] + p[3,] + p[4,]
    term002 <- p[1,] + p[3,] + p[5,] + p[7,]
    term003 <- p[2,] + p[4,] + p[6,] + p[8,]
    term004 <- p[5,] + p[6,] + p[7,] + p[8,]
    term005 <- p[1,] + p[2,] + p[5,] + p[6,]
    term006 <- p[3,] + p[4,] + p[7,] + p[8,]

    term01 <- sqrt(term001 * term002 * term003 * term004)
    term02 <- sqrt(term001 * term005 * term006 * term004)
    term03 <- (term001 +term002)*term003*term004*(-p[2,]*(p[5,] + p[7,]) - p[4,]*(p[5,] + p[7,]) + (p[1,] + p[3,])*(p[6,] + p[8,]))
    term04 <- (term001+term005)*term006*term004*(-p[3,]*(p[5,] + p[6,]) - p[4,]*(p[5,] + p[6,]) + (p[1,] + p[2,])*(p[7,] + p[8,]))
    term05 <- term002*(term001+term003)*term004*(-p[2,]*(p[5,] + p[7,]) - p[4,]*(p[5,] + p[7,]) + (p[1,] + p[3,])*(p[6,] + p[8,]))
    term06 <- term005*(term001+term006)*term004*(-p[3,]*(p[5,] + p[6,]) - p[4,]*(p[5,] + p[6,]) + (p[1,] + p[2,])*(p[7,] + p[8,]))
    term07 <- term001*term003*(term002+term004)*(-p[2,]*(p[5,] + p[7,]) - p[4,]*(p[5,] + p[7,]) + (p[1,] + p[3,])*(p[6,] + p[8,]))
    term08 <- term001*term006*(term004+term005)*(-p[3,]*(p[5,] + p[6,]) - p[4,]*(p[5,] + p[6,]) + (p[1,] + p[2,])*(p[7,] + p[8,]))
    term09 <- term001*term002*(term003+term004)*(p[2,]*(p[5,] + p[7,]) + p[4,]*(p[5,] + p[7,]) - (p[1,] + p[3,])*(p[6,] + p[8,]))
    term10 <- term001*term005*(term006+term004)*(p[3,]*(p[5,] + p[6,]) + p[4,]*(p[5,] + p[6,]) - (p[1,] + p[2,])*(p[7,] + p[8,]))

    d1 <- -(p[6,] + p[8,]) / term01 + (p[7,] + p[8,]) / term02 + term03 / (2*term01^3) - term04 / (2*term02^3)
    d2 <- (p[5,] + p[7,]) / term01 + (p[7,] + p[8,]) / term02 + term05 / (2*term01^3) - term04 / (2*term02^3)
    d3 <- -(p[6,] + p[8,]) / term01 -(p[5,] + p[6,]) / term02 + term03 / (2*term01^3) - term06 / (2*term02^3)
    d4 <- (p[5,] + p[7,]) / term01 -(p[5,] + p[6,]) / term02 + term05 / (2*term01^3) - term06 / (2*term02^3)
    d5 <- (p[2,] + p[4,]) / term01 - (p[3,] + p[4,]) / term02 + term07 / (2*term01^3) - term08 / (2*term02^3)
    d6 <- -(p[1,] + p[3,]) / term01 - (p[3,] + p[4,]) / term02 - term09 / (2*term01^3) - term08 / (2*term02^3)
    d7 <- (p[2,] + p[4,]) / term01 + (p[1,] + p[2,]) / term02 + term07 / (2*term01^3) + term10 / (2*term02^3)
    d8 <- -(p[1,] + p[3,]) / term01 + (p[1,] + p[2,]) / term02 - term09 / (2*term01^3) + term10 / (2*term02^3)

    term1a <- (-((p[3,] + p[4,]) * (p[5,] + p[6,])) + (p[1,] + p[2,]) * (p[7,] + p[8,]))
    term1b <- (-((p[2,] + p[4,]) * (p[5,] + p[7,])) + (p[1,] + p[3,]) * (p[6,] + p[8,]))

    ### simple method  ###
    variance1 <-  lapply(1:m, function(i) c(d1[i],d2[i],d3[i],d4[i],d5[i],d6[i],d7[i],d8[i]) %*% matrix(unlist(Sigma[i]),8,8) %*% c(d1[i],d2[i],d3[i],d4[i],d5[i],d6[i],d7[i],d8[i]))
    variance1 <- unlist(variance1)
    interval1 <- qnorm(1-alpha/2) * sqrt(variance1/n[j])
    k1 <- ifelse(MCC1 - MCC2 - interval1 < tMCC1 - tMCC2 & tMCC1 - tMCC2 < MCC1 - MCC2 + interval1, 1, 0)

    
    ###  Zou's method  ###
    ## the result of differentiation for psi_tilde
    numerator_s11 <- (term006 * term004 * (2 * term001 * term005 * (p[7,] + p[8,]) - (term001 + term005) * term1a))
    denominator_s11 <- (2 * (term001 * term005 * term006 * term004)^(1.5))
    s11 <- numerator_s11 / denominator_s11
    numerator_s12 <- (term003 * term004 * (2 * term001 * term002 * (p[6,] + p[8,]) - (term001 + term002) * term1b))
    denominator_s12 <- (2 * (term001 * term002 * term003 * term004)^(1.5))
    s12 <- numerator_s12 / denominator_s12
    s21 <- s11
    numerator_s22 <- (term002 * term004 * (-2 * term001 * (p[5,] + p[7,]) * term003 - (term001 + term003) * term1b))
    denominator_s22 <- denominator_s12
    s22 <- numerator_s22 / denominator_s22
    s32 <- s12
    numerator_s31 <- (term005 * term004 * (-2 * term001 * (p[5,] + p[6,]) * term006 - (term001 + term006) * term1a))
    denominator_s31 <- (2 * (term001 * term005 * term006 * term004)^(1.5))
    s31 <- numerator_s31 / denominator_s31
    s41 <- s31
    s42 <- s22
    numerator_s51 <- (term001 * term006 * (-2 * (p[3,] + p[4,]) * term005 * term004 - (term004 + term005) * term1a))
    denominator_s51 <- denominator_s31
    s51 <- numerator_s51 / denominator_s51
    numerator_s52 <- (term001 * term003 * (-2 * (p[2,] + p[4,]) * term002 * term004 - (term002 + term004) *term1b))
    denominator_s52 <- denominator_s12
    s52 <- numerator_s52 / denominator_s52
    s61 <- s51
    numerator_s62 <- (term001 * term002 * (2 * (p[1,] + p[3,]) * term003 * term004 - (term003 + term004) * term1b))
    denominator_s62 <- denominator_s12
    s62 <- numerator_s62 / denominator_s62
    numerator_s71 <- (term001 * term005 * (2 * (p[1,] + p[2,]) * term006 * term004 - (term004 + term006) * term1a))
    denominator_s71 <- denominator_s31
    s71 <- numerator_s71 / denominator_s71
    s72 <- s52
    s81 <- s71
    s82 <- s62

    psi_tilde_derivative <- lapply(1:m, function(i) matrix(c(s11[i],s12[i],s21[i],s22[i],s31[i],s32[i],s41[i],s42[i],s51[i],s52[i],s61[i],s62[i],s71[i],s72[i],s81[i],s82[i]),nrow=2,ncol=8))
    covariance <-  lapply(1:m, function(i) psi_tilde_derivative[[i]] %*% Sigma[[i]]  %*% t(psi_tilde_derivative[[i]]))
    cor <- lapply(1:m, function(i) covariance[[i]][1,2]/sqrt(covariance[[i]][1,1] * covariance[[i]][2,2]))
    cor <- unlist(cor)

    # apply Fisher's z transformation to MCC1 and calculate the confidence interval for MCC1
    z1 <- (1/2) * log((1+MCC1)/(1-MCC1))
    interval1 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(1/(1-MCC1[i]^2)^2 * covariance[[i]][1,1]/n[j]))
    interval1 <- unlist(interval1)
    zl1 <- z1 - interval1
    zu1 <- z1 + interval1
    l1 <- (exp(2*zl1)-1)/((exp(2*zl1)+1))  # the confidence lower limit for MCC1
    u1 <- (exp(2*zu1)-1)/((exp(2*zu1)+1))  # the confidence upper limit for MCC1

    # apply Fisher's z transformation to MCC2 and calculate the confidence interval for MCC2
    z2 <- (1/2) * log((1+MCC2)/(1-MCC2))
    interval2 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(1/(1-MCC2[i]^2)^2 * covariance[[i]][2,2]/n[j]))
    interval2 <- unlist(interval2)
    zl2 <- z2 - interval2
    zu2 <- z2 + interval2
    l2 <- (exp(2*zl2)-1)/((exp(2*zl2)+1))  # the confidence lower limit for MCC2
    u2 <- (exp(2*zu2)-1)/((exp(2*zu2)+1))  # the confidence upper limit for MCC2

    # Based on Zou (2007), calculate the confidence lower and upper limits for MCC1-MCC2
    L <- MCC1-MCC2 - sqrt((MCC1-l1)^2 + (u2-MCC2)^2 - 2* cor * (MCC1-l1)*(u2-MCC2))
    U <- MCC1-MCC2 + sqrt((u1-MCC1)^2 + (MCC2-l2)^2 - 2* cor * (MCC1-l1)*(u2-MCC2))

    k2 <- ifelse(L < tMCC1 - tMCC2 & tMCC1 - tMCC2 < U, 1, 0)
    
    
    ###  MT method  ###
    kappa <- MCC1 - MCC2
    g_derivative <- function(x) {
    return(2 / (4 - x^2))
    }
    variance3 <-  lapply(1:m, function(i) g_derivative(kappa[i])^2 * c(d1[i],d2[i],d3[i],d4[i],d5[i],d6[i],d7[i],d8[i]) %*% matrix(unlist(Sigma[i]),8,8) %*% c(d1[i],d2[i],d3[i],d4[i],d5[i],d6[i],d7[i],d8[i]))
    variance3 <- unlist(variance3)
    z <- (1/2) * log((2 + kappa) / (2 - kappa))
    interval3 <- qnorm(1-alpha/2) * sqrt(variance3/n[j])
    zL <- z - interval3
    zU <- z + interval3
    L <- 2*(exp(2*zL)-1)/((exp(2*zL)+1))
    U <- 2*(exp(2*zU)-1)/((exp(2*zU)+1))
    k3 <- ifelse(L < tMCC1 - tMCC2 & tMCC1 - tMCC2 < U, 1, 0)
    
    
    ## store each result in coverage_prob
    
    k1 <- k1[complete.cases(k1)]
    k2 <- k2[complete.cases(k2)]
    k3 <- k3[complete.cases(k3)]
    coverage_prob1[j] <- mean(k1)
    coverage_prob2[j] <- mean(k2)
    coverage_prob3[j] <- mean(k3)
}

result <- data.frame(row.names=c("n=50","n=100","n=500","n=1000","n=5000","n=10000"), Method_1=coverage_prob1, Method_2=coverage_prob2, Method_3=coverage_prob3)
print(paste("P(y=1)=", p_111+p_101+p_011+p_001, "MCC_1=", tMCC1, "MCC_2=", tMCC2))
options(digits = 4)
result









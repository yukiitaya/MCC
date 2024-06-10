### 100(1-alpha)% confidence intervals for paired designs ###

m <- 1000000 # the number of times to calculate the confidence interval
alpha <- 0.05 # significance level
n <- c(50,100,500,1000,5000,10000) # sample sizes

# Initialize vectors to store coverage probabilities and counts of undefined cases
coverage_prob1 <- numeric(length(n))  # a vector for the Simple method
coverage_prob2 <- numeric(length(n))  # a vector for the Zou's method
coverage_prob3 <- numeric(length(n))  # a vector for the MT method
na_count1 <- numeric(length(n)) # count of NA values for the Simple method
na_count2 <- numeric(length(n)) # count of NA values for Zou's method
na_count3 <- numeric(length(n)) # count of NA values for the MT method

# Initialize vectors to store results for each simulation
k1 <- numeric(m)
k2 <- numeric(m)
k3 <- numeric(m)
set.seed(2023) # Set seed for reproducibility

# True values of TP, FP, FN, TN for two conditions
tTP1 <- 0.4
tFP1 <- 0.1
tFN1 <- 0.1
tTN1 <- 1 - tTP1 - tFP1 - tFN1
tTP2 <- 0.45
tFP2 <- 0.05
tFN2 <- 0.05
tTN2 <- 1 - tTP2 - tFP2 - tFN2

p_110 = 0.01
p_001 = 0.001
p_100 <- tFP1 - p_110
p_010 <- tFP2 - p_110
p_011 <- tFN1 - p_001
p_101 <- tFN2 - p_001
p_111 <- tTP1 - p_101
p_000 <- tTN1 - p_010


# Calculate true MCC values for two conditions
tMCC1 <- (tTP1*tTN1-tFP1*tFN1)/sqrt((tTP1+tFP1)*(tTP1+tFN1)*(tTN1+tFP1)*(tTN1+tFN1))
tMCC2 <- (tTP2*tTN2-tFP2*tFN2)/sqrt((tTP2+tFP2)*(tTP2+tFN2)*(tTN2+tFP2)*(tTN2+tFN2))

# Function to generate the limiting variance of p_hat
phi_matrix <- function(p) {
  mat <- -outer(p, p)
  diag(mat) <- p * (1-p)
  return(mat)
}

# Simulation loop for each sample size
for(j in 1:length(n)){
    # Data generation
    data <- rmultinom(m, size=n[j], prob = c(p_111,p_110,p_101,p_100,p_011,p_010,p_001,p_000))
    p <- data/n[j]
    
    # Recalculate probabilities from simulated data
    p111 <- p[1,]
    p110 <- p[2,]
    p101 <- p[3,]
    p100 <- p[4,]
    p011 <- p[5,]
    p010 <- p[6,]
    p001 <- p[7,]
    p000 <- p[8,]
    
    # Compute TP, FP, FN, TN based on recalculated probabilities
    TP1 <- p111 + p101  # p1.1
    FP1 <- p110 + p100  # p1.0
    FN1 <- p011 + p001  # p0.1
    TN1 <- p010 + p000  # p0.0
    TP2 <- p111 + p011  # p.11
    FP2 <- p110 + p010  # p.10
    FN2 <- p101 + p001  # p.01
    TN2 <- p100 + p000  # p.00
    
    # Calculate MCC from simulated data
    MCC1 <- (TP1*TN1-FP1*FN1)/sqrt((TP1+FP1)*(TP1+FN1)*(FP1+TN1)*(FN1+TN1))
    MCC2 <- (TP2*TN2-FP2*FN2)/sqrt((TP2+FP2)*(TP2+FN2)*(FP2+TN2)*(FN2+TN2))

    phi_p <-  lapply(1:m, function(i) phi_matrix(p[,i]))  # phi(p)

    p1.. <- p111 + p110 + p101 + p100
    p0.. <- p011 + p010 + p001 + p000
    p.1. <- p111 + p110 + p011 + p010
    p.0. <- p101 + p100 + p001 + p000
    p..1 <- p111 + p101 + p011 + p001
    p..0 <- p110 + p100 + p010 + p000
    
    D1 <- sqrt(p1.. * p..1 * p0.. * p..0)
    D2 <- sqrt(p.1. * p..1 * p.0. * p..0)
    
    # Partial derivatives for MCC1
    dMCC1_p111 <- dMCC1_p101 <- TN1/D1 - (p..1+p1..)/(2*p..1*p1..)*MCC1
    dMCC1_p110 <- dMCC1_p100 <- -FN1/D1 - (p1..+p..0)/(2*p1..*p..0)*MCC1
    dMCC1_p011 <- dMCC1_p001 <- -FP1/D1 - (p..1+p0..)/(2*p..1*p0..)*MCC1
    dMCC1_p010 <- dMCC1_p000 <- TP1/D1 - (p0..+p..0)/(2*p0..*p..0)*MCC1
    
    # Partial derivatives for MCC2
    dMCC2_p111 <- dMCC2_p011 <- TN2/D2 - (p..1+p.1.)/(2*p..1*p.1.)*MCC2
    dMCC2_p110 <- dMCC2_p010 <- -FN2/D2 - (p.1.+p..0)/(2*p.1.*p..0)*MCC2
    dMCC2_p101 <- dMCC2_p001 <- -FP2/D2 - (p..1+p.0.)/(2*p..1*p.0.)*MCC2
    dMCC2_p100 <- dMCC2_p000 <- TP2/D2 - (p.0.+p..0)/(2*p.0.*p..0)*MCC2
    
    # Partial derivatives for the difference between MCC1 and MCC2
    dpsi_p111 <- dMCC1_p111 - dMCC2_p111
    dpsi_p110 <- dMCC1_p110 - dMCC2_p110
    dpsi_p101 <- dMCC1_p101 - dMCC2_p101
    dpsi_p100 <- dMCC1_p100 - dMCC2_p100
    dpsi_p011 <- dMCC1_p011 - dMCC2_p011
    dpsi_p010 <- dMCC1_p010 - dMCC2_p010
    dpsi_p001 <- dMCC1_p001 - dMCC2_p001
    dpsi_p000 <- dMCC1_p000 - dMCC2_p000
    
    ## Simple method ##
    # derivative for psi
    dpsi <- lapply(1:m, function(i) c(dpsi_p111[i],dpsi_p110[i],dpsi_p101[i],dpsi_p100[i],dpsi_p011[i],dpsi_p010[i],dpsi_p001[i],dpsi_p000[i]))
    # the limiting variance
    variance_simple <-  lapply(1:m, function(i) {
        dpsi[[i]] %*% phi_p[[i]] %*% dpsi[[i]]
    })
    variance_simple <- unlist(variance_simple)
    interval_simple <- qnorm(1-alpha/2) * sqrt(variance_simple/n[j])
    k1 <- ifelse(MCC1 - MCC2 - interval_simple < tMCC1 - tMCC2 & tMCC1 - tMCC2 < MCC1 - MCC2 + interval_simple, 1, 0)

    ## Zou's method ##
    # derivative for psi_tilde
    dpsi_tilde_p111 <- c(dMCC1_p111, dMCC2_p111)
    dpsi_tilde_p110 <- c(dMCC1_p110, dMCC2_p110)
    dpsi_tilde_p101 <- c(dMCC1_p101, dMCC2_p101)
    dpsi_tilde_p100 <- c(dMCC1_p100, dMCC2_p100)
    dpsi_tilde_p011 <- c(dMCC1_p011, dMCC2_p011)
    dpsi_tilde_p010 <- c(dMCC1_p010, dMCC2_p010)
    dpsi_tilde_p001 <- c(dMCC1_p001, dMCC2_p001)
    dpsi_tilde_p000 <- c(dMCC1_p000, dMCC2_p000)

    dpsi_tilde <- lapply(1:m, function(i)
                            matrix(c(dMCC1_p111[i], dMCC2_p111[i],
                                     dMCC1_p110[i], dMCC2_p110[i],
                                     dMCC1_p101[i], dMCC2_p101[i],
                                     dMCC1_p100[i], dMCC2_p100[i],
                                     dMCC1_p011[i], dMCC2_p011[i],
                                     dMCC1_p010[i], dMCC2_p010[i],
                                     dMCC1_p001[i], dMCC2_p001[i],
                                     dMCC1_p000[i], dMCC2_p000[i]),
                                    nrow=2,ncol=8))
    covariance <-  lapply(1:m, function(i) dpsi_tilde[[i]] %*% phi_p[[i]]  %*% t(dpsi_tilde[[i]]))
    cor <- lapply(1:m, function(i) covariance[[i]][1,2]/sqrt(covariance[[i]][1,1] * covariance[[i]][2,2]))
    cor <- unlist(cor)

    # Fisher's z method for MCC1
    xi1 <- (1/2) * log((1+MCC1)/(1-MCC1))
    interval1 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(1/(1-MCC1[i]^2)^2 * covariance[[i]][1,1]/n[j]))
    interval1 <- unlist(interval1)
    zl1 <- xi1 - interval1
    zu1 <- xi1 + interval1
    l1 <- (exp(2*zl1)-1)/((exp(2*zl1)+1))
    u1 <- (exp(2*zu1)-1)/((exp(2*zu1)+1))

    # Fisher's z method for MCC2
    xi2 <- (1/2) * log((1+MCC2)/(1-MCC2))
    interval2 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(1/(1-MCC2[i]^2)^2 * covariance[[i]][2,2]/n[j]))
    interval2 <- unlist(interval2)
    zl2 <- xi2 - interval2
    zu2 <- xi2 + interval2
    l2 <- (exp(2*zl2)-1)/((exp(2*zl2)+1))
    u2 <- (exp(2*zu2)-1)/((exp(2*zu2)+1))

    # Zou's approach for the confidence interval
    L <- MCC1-MCC2 - sqrt((MCC1-l1)^2 + (u2-MCC2)^2 - 2* cor * (MCC1-l1)*(u2-MCC2))
    U <- MCC1-MCC2 + sqrt((u1-MCC1)^2 + (MCC2-l2)^2 - 2* cor * (MCC1-l1)*(u2-MCC2))

    k2 <- ifelse(L < tMCC1 - tMCC2 & tMCC1 - tMCC2 < U, 1, 0)



    ## MT method ##
    difference <- MCC1 - MCC2
    g <- function(x) {
    return(2 / (4 - x^2))
    }
    variance_MT <-  lapply(1:m, function(i) g(difference[i])^2 * dpsi[[i]] %*% phi_p[[i]] %*% dpsi[[i]])
    variance_MT <- unlist(variance_MT)
    xi <- (1/2) * log((2 + difference) / (2 - difference))
    interval_MT <- qnorm(1-alpha/2) * sqrt(variance_MT/n[j])
    zL <- xi - interval_MT
    zU <- xi + interval_MT
    L <- 2*(exp(2*zL)-1)/((exp(2*zL)+1))
    U <- 2*(exp(2*zU)-1)/((exp(2*zU)+1))
    k3 <- ifelse(L < tMCC1 - tMCC2 & tMCC1 - tMCC2 < U, 1, 0)
    

    na_count1[j] <- sum(is.na(k1))
    na_count2[j] <- sum(is.na(k2))
    na_count3[j] <- sum(is.na(k3))

    k1 <- k1[complete.cases(k1)]
    k2 <- k2[complete.cases(k2)]
    k3 <- k3[complete.cases(k3)]
    coverage_prob1[j] <- mean(k1)
    coverage_prob2[j] <- mean(k2)
    coverage_prob3[j] <- mean(k3)
}

# Display result                           
result <- data.frame(row.names=c("n=50","n=100","n=500","n=1000","n=5000","n=10000"), Simple=coverage_prob1, na1=na_count1, Zou=coverage_prob2, na2=na_count2, MT=coverage_prob3, na3=na_count3)
options(digits = 2)
print(paste("P(y=1)=", p_111+p_110+p_101+p_100, "MCC_1=", tMCC1, "MCC_2=", tMCC2))
options(digits = 4)
result
                           

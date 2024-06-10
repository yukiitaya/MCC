### 100(1-alpha)% confidence intervals for single cases ###

m <- 1000000  # the number of times to calculate the confidence interval
alpha <- 0.05
n <- c(50,100,500,1000,5000,10000)  # sample sizes
coverage_prob1 <- numeric(length(n))  # a vector for the simple method
coverage_prob2 <- numeric(length(n))  # a vector for Fisher's z method
coverage_prob3 <- numeric(length(n))  # a vector for Fisher's z transformation with the asymptotic variance 1/(n-3)

# Vectors to count the number of undefined cases
na_count1 <- numeric(length(n))
na_count2 <- numeric(length(n))
na_count3 <- numeric(length(n))
set.seed(2023)

# True values of TP, FP, FN, TN
tTP <- 0.45  # p11
tFP <- 0.05  # p01
tFN <- 0.05  # p10
tTN <- 1 - tTP - tFP - tFN  # p00

# True value of MCC
tMCC <- (tTP * tTN - tFP * tFN) / sqrt((tTP + tFP) * (tTP + tFN) * (tFP + tTN) * (tFN + tTN))
tMCC

# Function to generate the limiting variance of p_hat
variance_matrix <- function(p) {
  mat <- -outer(p, p)
  diag(mat) <- p * (1 - p)
  return(mat)
}

for (j in 1:length(n)) {
  # Data generation
  data <- rmultinom(m, size = n[j], prob = c(tTP, tFP, tFN, tTN))
  p11 <- data[1, ] / n[j]
  p01 <- data[2, ] / n[j]
  p10 <- data[3, ] / n[j]
  p00 <- data[4, ] / n[j]

  p1. <- p11 + p10
  p0. <- p01 + p00
  p.1 <- p11 + p01
  p.0 <- p10 + p00
  
  D <- sqrt(p1.*p0.*p.1*p.0)
  MCC <- (p11 * p00 - p10 * p01) / D
  
  dphi_p11 <- p00/D - (p1.+p.1)/(2*p1.*p.1)*MCC   # differentiation with respect to p11
  dphi_p01 <- -p10/D - (p.1+p0.)/(2*p.1*p0.)*MCC  # differentiation with respect to p01
  dphi_p10 <- -p01/D - (p1.+p.0)/(2*p1.*p.0)*MCC  # differentiation with respect to p10
  dphi_p00 <- p11/D - (p0.+p.0)/(2*p0.*p.0)*MCC   # differentiation with respect to p00

  # differentiation of phi and the composition of f and phi
  dphi <- lapply(1:m, function(i) c(dphi_p11[i], dphi_p01[i], dphi_p10[i], dphi_p00[i]))
  dfphi <- lapply(1:m, function(i) 1/(1-MCC[i]^2)*dphi[[i]])

  # Calculate the limiting variance of p_hat
  sigma <- lapply(1:m, function(i) variance_matrix(data[,i] / n[j]))

  ## Limiting variance of MCC calculated by the simple method
  var_simple <- lapply(1:m, function(i) dphi[[i]] %*% sigma[[i]] %*% dphi[[i]] / n[j])
  var_simple <- unlist(var_simple)

  ## Simple method
  k <- numeric(m)
  interval1 <- lapply(1:m, function(i) qnorm(1 - alpha / 2) * sqrt(var_simple[i]))
  interval1 <- unlist(interval1)
  k <- ifelse(MCC - interval1 < tMCC & tMCC < MCC + interval1, 1, 0)
  na_count1[j] <- sum(is.na(k))
  k <- k[complete.cases(k)]
  coverage_prob1[j] <- mean(k)

  ## Fisher's z method
  z <- (1 / 2) * log((1 + MCC) / (1 - MCC))
  var_fisher1 <- lapply(1:m, function(i) dfphi[[i]] %*% sigma[[i]] %*% dfphi[[i]])
  var_fisher1 <- unlist(var_fisher1)
  interval2 <- qnorm(1 - alpha / 2) * sqrt(var_fisher1 / n[j])
  zL1 <- z - interval2
  zU1 <- z + interval2
  L1 <- (exp(2 * zL1) - 1) / (exp(2 * zL1) + 1)
  U1 <- (exp(2 * zU1) - 1) / (exp(2 * zU1) + 1)
  l1 <- numeric(m)
  l1 <- ifelse(L1 < tMCC & tMCC < U1, 1, 0)
  na_count2[j] <- sum(is.na(l1))
  l1 <- l1[complete.cases(l1)]
  coverage_prob2[j] <- mean(l1)

  ## Fisher's z transformation with the asymptotic variance 1/(n-3)
  var_fisher2 <- 1 / (n[j]-3)
  interval3 <- qnorm(1 - alpha / 2) / sqrt(n[j]-3)
  zL2 <- z - interval3
  zU2 <- z + interval3
  L2 <- (exp(2 * zL2) - 1) / (exp(2 * zL2) + 1)
  U2 <- (exp(2 * zU2) - 1) / (exp(2 * zU2) + 1)
  l2 <- numeric(m)
  l2 <- ifelse(L2 < tMCC & tMCC < U2, 1, 0)
  na_count3[j] <- sum(is.na(l2))
  l2 <- l2[complete.cases(l2)]
  coverage_prob3[j] <- mean(l2)
}

# Display results
options(digits = 4)
result <- data.frame(row.names = c("n=50", "n=100", "n=500", "n=1000", "n=5000", "n=10000"),
                     Simple_method = coverage_prob1,
                     na1 = na_count1,
                     Fishers_z_method = coverage_prob2,
                     na2 = na_count2,
                     Fishers_z_trans = coverage_prob3,
                     na3 = na_count3)

print(paste("P(y=1)=", tTP + tFN, "MCC=", tMCC))
result

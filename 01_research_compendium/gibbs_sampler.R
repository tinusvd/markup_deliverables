source("./project_files/packages.R")

# Simulated data
set.seed(1234)
X <- rnorm(300, 10, 3)
M <- rnorm(300, 0.21*X + 20, 5)
Y <- rnorm(300, 0.13*M + 0.11*X + 10, 2) 
data <- data.frame(X = X, M = M, Y = Y, N = length(X))

# Hyperparameters
## M = B0 + B1*X + e1
mu0 <-0 
sigma0 <- 10 
mu1 <- 0 
sigma1 <- 10 
alpha00 <- 1 
beta00 <- 1 

## Y = B2 + B3*m + B4*x + e2
mu2 <- 0 
sigma2 <- 10 
mu3 <- 0
sigma3 <- 10
mu4 <- 0
sigma4 <- 10
alpha01 <- 2
beta01 <- 2

combined_gibbs <- function(mu0, sigma0, mu1, sigma1, alpha00, beta00,
                           mu2, sigma2, mu3, sigma3, mu4, sigma4, alpha01, beta01,
                           iter = 215000, burnin = 15000, thin = 40,
                           X, M, Y) {
  # Initial values for parameters
  b0 <- 1  
  b1 <- 1  
  s1 <- 1  
  b2 <- 1 
  b3 <- 1 
  b4 <- 1  
  s2 <- 1 
  # Initialize empty data sets
  N_M <- length(M) 
  N_Y <- length(Y) 
  # Desired size of posterior samples
  total_samples <- 5000  
  # Matrix for storing the parameter estimates
  vals <- matrix(NA, nrow = total_samples, ncol = 7)
  
  sample_index <- 1
  
  for (i in 1:iter) {
    # Update for the first model: M = b0 + b1*X + e1
    # Intercept of M
    # Sample auxiliary g for robust estimates
    g0 <- 1 / rgamma(1, shape = 1 / 2, rate = N_M / 2)
    post.mu0 <- ((sum(M - b1 * X)) / s1 + mu0 / (g0 * sigma0)) / (N_M / s1 + 1 / (g0 * sigma0))
    post.sigma0 <- 1 / (N_M / s1 + 1 / (g0 * sigma0))
    b0 <- rnorm(1, mean = post.mu0, sd = sqrt(post.sigma0))
    
    # Regression coefficient of X: alpha
    # Auxiliary variable g for robuster estimates
    g1 <- 1 / rgamma(1, shape = 1 / 2, rate = N_M / 2)
    post.mu1 <- ((sum(X * (M - b0))) / s1) / (sum(X^2) / s1 + 1 / (g1 * sigma1))
    post.sigma1 <- 1 / (sum(X^2) / s1 + 1 / (g1 * sigma1))
    b1 <- rnorm(1, mean = post.mu1, sd = sqrt(post.sigma1))
    
    # Error variance of M
    post.alpha0 <- (N_M / 2) + alpha00
    post.beta0 <- ((sum((M - (b0 + b1 * X))^2)) / 2) + beta00
    s1 <- 1 / rgamma(1, shape = post.alpha0, rate = post.beta0)
    
    # Update for the second model: Y = b2 + b3*M + b4*X + e2
    # Intercept for Y
    # Auxiliary variable g for robuster estimates
    g2 <- 1 / rgamma(1, shape = 1 / 2, rate = N_Y / 2)
    post.mu2 <- ((sum(Y - b3 * M - b4 * X)) / s2 + mu2 / (g2 * sigma2)) / (N_Y / s2 + 1 / (g2 * sigma2))
    post.sigma2 <- 1 / (N_Y / s2 + 1 / (g2 * sigma2))
    b2 <- rnorm(1, mean = post.mu2, sd = sqrt(post.sigma2))
    
    # Regression coefficient for M: beta
    # Auxiliary variable g for robuster estimates
    g3 <- 1 / rgamma(1, shape = 1 / 2, rate = N_Y / 2)
    post.mu3 <- ((sum(M * (Y - b2 - b4 * X))) / s2) / (sum(M^2) / s2 + 1 / (g3 * sigma3))
    post.sigma3 <- 1 / (sum(M^2) / s2 + 1 / (g3 * sigma3))
    b3 <- rnorm(1, mean = post.mu3, sd = sqrt(post.sigma3))
    
    # Regression coefficient for X: tau'
    # Auxiliary variable g for robuster estimates
    g4 <- 1 / rgamma(1, shape = 1 / 2, rate = N_Y / 2)
    post.mu4 <- ((sum(X * (Y - b2 - b3 * M))) / s2) / (sum(X^2) / s2 + 1 / (g4 * sigma4))
    post.sigma4 <- 1 / (sum(X^2) / s2 + 1 / (g4 * sigma4))
    b4 <- rnorm(1, mean = post.mu4, sd = sqrt(post.sigma4))
    
    # Error variance for Y
    post.alpha01 <- (N_Y / 2) + alpha01
    post.beta01 <- ((sum((Y - (b2 + b3 * M + b4 * X))^2)) / 2) + beta01
    s2 <- 1 / rgamma(1, shape = post.alpha01, rate = post.beta01)
    
    # Thinning interval of 40 for minimizing of autocorrelation at lags 1-40
    if (i > burnin && (i - burnin) %% thin == 0) {
      vals[sample_index, ] <- c(b0, b1, s1, b2, b3, b4, s2)
      sample_index <- sample_index + 1
      if (sample_index > total_samples) break
    }
  }
  
  return(vals)
}

# Run the sampler three times for three chains
set.seed(1234)
chain1 <- combined_gibbs(
  mu0,
  sigma0,
  mu1,
  sigma1,
  alpha00,
  beta00,
  mu2,
  sigma2,
  mu3,
  sigma3,
  mu4,
  sigma4,
  alpha01,
  beta01,
  X = X,
  M = M,
  Y = Y
)
set.seed(1243)
chain2 <- combined_gibbs(
  mu0,
  sigma0,
  mu1,
  sigma1,
  alpha00,
  beta00,
  mu2,
  sigma2,
  mu3,
  sigma3,
  mu4,
  sigma4,
  alpha01,
  beta01,
  X = X,
  M = M,
  Y = Y
)
set.seed(2134)
chain3 <- combined_gibbs(
  mu0,
  sigma0,
  mu1,
  sigma1,
  alpha00,
  beta00,
  mu2,
  sigma2,
  mu3,
  sigma3,
  mu4,
  sigma4,
  alpha01,
  beta01,
  X = X,
  M = M,
  Y = Y
)

# Multiplying the alpha and beta for the indirect effect for each chain
abchain1 <- chain1[, 2] * chain1[, 5]
abchain2 <- chain2[, 2] * chain2[, 5]
abchain3 <- chain3[, 2] * chain3[, 5]

# bind the new columns to their chains
chain_1 <- cbind(chain1, abchain1)
chain_2 <- cbind(chain2, abchain2)
chain_3 <- cbind(chain3, abchain3)

# Giving the columns names
column_names_chains <- c("Intercept M",
                         "Alpha",
                         "Sigma_M",
                         "Intercept Y",
                         "Beta",
                         "C'",
                         "Sigma_Y",
                         "ab")

colnames(chain_1) <- column_names_chains
colnames(chain_2) <- column_names_chains
colnames(chain_3) <- column_names_chains

# Create one data frame of the three chains
finalchain <- as.data.frame(rbind(chain_1, chain_2, chain_3))

source("./project_files/AAFBF_estimator.R")

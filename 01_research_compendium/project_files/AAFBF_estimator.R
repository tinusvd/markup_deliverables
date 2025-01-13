########## SAMPLES ##########################################################
# Small data set with only parameters of first model
finalchainM <- finalchain[, 1:3]
# Small data set with only parameters of second model
finalchainY <- finalchain[, 4:7]

################## SPECIFICATION FOR BF ESTIMATION #############################

# Creating mean structures from the posterior
posterior_mean_alpha <- colMeans(finalchainM)
posterior_mean_beta <- colMeans(finalchainY)

# Creating covariance matrices from the posterior
posterior_cov_alpha <- cov(finalchainM)
posterior_cov_beta <- cov(finalchainY)

# Focal points for the parameters of interest
prior_meanh_alpha <- c(0, 0, 0)
prior_meanh_beta <- c(0, 0, 0, 0)

# Specifying the prior covariance matrices
prior_cov_alpha <- posterior_cov_alpha
prior_cov_beta <- posterior_cov_beta

# Setting lower limits for integration. 0 for alpha, because h2: \alpha > 0
lower_limits_alpha_g0 <- c(-Inf, 0, -Inf)
# Setting upper limit to infinity for alpha
upper_limits_alpha_g0 <- c(Inf, Inf, Inf)

# Setting lower limits for integration. 0 for beta, because h2: \beta > 0
lower_limits_beta_g0 <- c(-Inf, 0, -Inf, -Inf)
# Setting upper limit to infinity for beta
upper_limits_beta_g0 <- c(Inf, Inf, Inf, Inf)


# Define the mean vector and covariance matrix for the fractional prior distribution
mean_prior_alpha <- rep(0, length(posterior_mean_alpha))
mean_prior_beta <- rep(0, length(posterior_mean_beta))

# J/N = 2/length(X)
J <- 2
fraction <- J / length(X)

# Multiplying the covariance matrices by the fraction
prior_cov_alpha <- fraction * posterior_cov_alpha
prior_cov_beta <- fraction * posterior_cov_beta


################################### ESTIMATION #################################


# Calculate the proportion of probability under the bivariate normal distribution for fit
fit_alpha_g0 <- pmvnorm(
  lower = lower_limits_alpha_g0,
  upper = upper_limits_alpha_g0,
  mean = posterior_mean_alpha,
  sigma = posterior_cov_alpha
)
fit_beta_g0 <- pmvnorm(
  lower = lower_limits_beta_g0,
  upper = upper_limits_beta_g0,
  mean = posterior_mean_beta,
  sigma = posterior_cov_beta
)

# Calculate the proportion of probability under the bivariate normal distribution for complexity
complexity_alpha_g0 <- pmvnorm(
  lower = lower_limits_alpha_g0,
  upper = upper_limits_alpha_g0,
  mean = mean_prior_alpha,
  sigma = prior_cov_alpha
)
complexity_beta_g0 <- pmvnorm(
  lower = lower_limits_beta_g0,
  upper = upper_limits_beta_g0,
  mean = mean_prior_beta,
  sigma = prior_cov_beta
)

# Calculate the Bayes Factor
BF_alpha <- fit_alpha_g0[1] / complexity_alpha_g0[1]
BF_beta <- fit_beta_g0[1] / complexity_beta_g0[1]

########################### CROSS-VALIDATION ###################################
set.seed(1234)
testlmM <- lm(M ~ X, data = data)
bain_alpha <- bain(testlmM, "X > 0")
BF_bain_alpha <- bain_alpha$fit["BF.u"][1,]
PMP_bain_alpha <- bain_alpha$fit["PMPb"][1,]


testlmYXM <- lm(Y ~ M + X, data = data)
bain_beta <- bain(testlmYXM, "M > 0")
BF_bain_beta <- bain_beta$fit["BF.u"][1,]
PMP_bain_beta <- bain_beta$fit["PMPb"][1,]

######################### BFmed CALCULATION ####################################

# Specify prior odds (equal belief in both hypotheses, reflecting uncertainty)
priorodds_alpha <- 1
priorodds_beta <- 1
priorodds_med <- priorodds_alpha * priorodds_beta / (1 + priorodds_alpha + priorodds_beta)

# Calculate posterior odds
postodds_alpha <- priorodds_alpha * BF_alpha
postodds_beta <- priorodds_beta * BF_beta

# Compute posterior probabilities
posteriorprob_alpha <- postodds_alpha / (1 + postodds_alpha)
posteriorprob_beta <- postodds_beta / (1 + postodds_beta)

# Calculate mediation effect BF
BF_med <- BF_alpha * BF_beta

# Compute the posterior probability for the mediation effect
priorodds_med <- 1
posteriorodds_med <- priorodds_med * BF_med
posteriorprob_med <- posteriorodds_med / (1 + posteriorodds_med)

cat("Gibbs sampler results")
cat("The Bayes Factor for Alpha is:", round(BF_alpha, 2))
cat("The PMP for Alpha is:", round(posteriorprob_alpha, 3))
cat("The Bayes Factor for Beta is:", round(BF_beta, 2))
cat("The PMP for Beta is:", round(posteriorprob_beta, 3))
cat("The Bayes Factor for Alpha*Beta is:", round(BF_med, 2))
cat("The PMP for Alpha*Beta is:", round(posteriorprob_med, 3))

cat("Bain results")
cat("The Bayes Factor for Alpha with bain is:", round(BF_bain_alpha, 2))
cat("The PMP for Alpha with bain is:", round(PMP_bain_alpha, 3))
cat("The Bayes Factor for Beta with bain is:", round(BF_bain_beta, 2))
cat("The PMP for Beta with bain is:", round(PMP_bain_beta, 3))
cat("The Bayes Factor for alpha*beta with bain is:", round(BF_bain_alpha*BF_bain_beta, 2))




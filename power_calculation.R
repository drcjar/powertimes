# Required libraries
library(survival)
library(MASS)

# Model 
# h(t) The hazard function, represents the instantaneous rate of occurrence of the event of interest (e.g., lung cancer death) at a specific time t, given that the individual has survived up to that time
# h(t) = h0(t) * exp(β1 * Occupation + β2 * Smoking + β3 * (Occupation * Smoking) + Σ (from i=4 to k) βi * OtherCovariates_i) where k is total number of covariates


# Parameters
n_sim <- 1000 # N of simulations
beta_occupation <- log(2) # 'doubling of risk' accepted criterion for attribution
beta_smoking <- log(5) # Odds ratio for 'ever smoking' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3505152/
beta_interaction <- log(2)  # Interaction effect
baseline_hazard <- 0.0008  # Baseline hazard (lung cancer rate in population)
censor_time <- 10  # Follow-up time in years
alpha <- 0.05  # Significance level
power <- 0.8  # Desired power
prop_smokers <- 0.3 # Proportion of smokers
prop_occupation <- 0.2  # Proportion of high-risk occupations
R2_other <- 0.1  # Proportion of variance explained by other covariates for lung cancer (estimated, e.g. https://onlinelibrary.wiley.com/doi/full/10.1002/ijc.30084)

# Initial sample size and step size
start_n <- 10000
step_size <- 1000


# Function to simulate data and calculate p-value for interaction term
simulate_power_cox <- function(n, beta_occupation, beta_smoking, beta_interaction, 
                                 baseline_hazard, censor_time, prop_smokers, 
                                 prop_occupation, R2_other) {
  # Simulate covariates
  occupation <- rbinom(n, 1, prop_occupation)
  smoking <- rbinom(n, 1, prop_smokers)

  # Simulate other covariates based on R2_other
  X_other <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, R2_other, R2_other, 1), ncol = 2))
  
  # Calculate combined effects
  hazard <- baseline_hazard * exp(beta_occupation * occupation + 
                                   beta_smoking * smoking + 
                                   beta_interaction * occupation * smoking +
                                   X_other[,1] + X_other[,2])  # Include other covariates
  
  # Simulate survival times and censoring
  survival_time <- rexp(n, hazard)
  censor_time <- runif(n, 0, censor_time)
  event_time <- pmin(survival_time, censor_time)
  event_status <- as.numeric(survival_time <= censor_time)
  
  # Fit Cox model and get p-value for interaction term
  cox_model <- coxph(Surv(event_time, event_status) ~ occupation * smoking + X_other[,1] + X_other[,2])
  summary(cox_model)$coefficients["occupation:smoking", "Pr(>|z|)"]
}


# Function to estimate power
estimate_power <- function(n) {
  p_values <- replicate(n_sim, simulate_power_cox(n, beta_occupation, beta_smoking, 
                                                    beta_interaction, baseline_hazard, 
                                                    censor_time, prop_smokers, 
                                                    prop_occupation, R2_other))
  mean(p_values < alpha)
}

# Find required sample size
find_sample_size <- function(start_n, step_size, target_power) {
  n <- start_n
  while (TRUE) {
    estimated_power <- estimate_power(n)
    cat("Current sample size:", n, "Power:", estimated_power, "\n")
    if (estimated_power >= target_power) {
      return(n)
    }
    n <- n + step_size
  }
}


# Estimate required sample size
required_sample_size <- find_sample_size(start_n, step_size, power)
cat("Required sample size:", required_sample_size, "\n")


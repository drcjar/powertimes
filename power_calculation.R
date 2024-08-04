# Required libraries
library(survival)
library(MASS)

# Parameters
n_sim <- 100  # Reduced number of simulations for testing
beta_occupation <- log(2)  # 'Doubling of risk' accepted criterion for attribution
beta_smoking <- log(5)  # Odds ratio for 'ever smoking' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3505152/
beta_interaction <- log(2)  # Interaction effect
baseline_hazard <- 0.0008  # Baseline hazard (lung cancer rate in population)
censor_time <- 10  # Follow-up time in years
alpha <- 0.05  # Significance level
power <- 0.8  # Desired power
prop_smokers <- 0.3  # Proportion of smokers
n_occupations <- 9  # Number of occupational categories
prop_occupation <- rep(1 / n_occupations, n_occupations)  # Equal proportion for each occupational category

# Initial sample size and step size
start_n <- 10000
step_size <- 10000

# Function to simulate data and calculate p-value for interaction term
simulate_power_cox <- function(n, beta_occupation, beta_smoking, beta_interaction, 
                               baseline_hazard, censor_time, prop_smokers, 
                               prop_occupation, n_occupations) {
  # Simulate covariates
  occupation <- sample(1:n_occupations, n, replace = TRUE, prob = prop_occupation)
  smoking <- rbinom(n, 1, prop_smokers)
  age <- rnorm(n, mean = 50, sd = 10)  # Assuming a normal distribution for age
  sex <- rbinom(n, 1, 0.5)  # Assuming equal proportion of males and females
  
  # Simulate random error term
  error_term <- rnorm(n, mean = 0, sd = 1)  # Mean 0 and some standard deviation
  
  # Calculate hazard with added error term
  hazard <- baseline_hazard * exp(beta_occupation * as.numeric(occupation == 1) + 
                                beta_smoking * smoking + 
                                beta_interaction * as.numeric(occupation == 1) * smoking + 
                                error_term)
  
  # Simulate survival times and censoring
  survival_time <- rexp(n, rate = hazard)
  censor_time <- runif(n, 0, censor_time)
  event_time <- pmin(survival_time, censor_time)
  event_status <- as.numeric(survival_time <= censor_time)
  
  # Fit Cox model and get p-value for interaction term
  cox_model <- coxph(Surv(event_time, event_status) ~ factor(occupation) * smoking + age + sex)
  
  # Extract p-value for interaction term
  term_names <- rownames(summary(cox_model)$coefficients)
  interaction_terms <- grep("factor\\(occupation\\)[0-9]+:smoking", term_names, value = TRUE)
  if (length(interaction_terms) == 0) {
    return(NA)
  } else {
    p_values <- summary(cox_model)$coefficients[interaction_terms, "Pr(>|z|)"]
    if (length(p_values) == 0) {
      return(NA)
    } else {
      return(mean(p_values, na.rm = TRUE))
    }
  }
}

# Function to estimate power
estimate_power <- function(n) {
  p_values <- replicate(n_sim, simulate_power_cox(n, beta_occupation, beta_smoking, 
                                                  beta_interaction, baseline_hazard, 
                                                  censor_time, prop_smokers, 
                                                  prop_occupation, n_occupations))
  mean(p_values < alpha, na.rm = TRUE)
}

# Find required sample size
find_sample_size <- function(start_n, step_size, target_power) {
  n <- start_n
  while (TRUE) {
    estimated_power <- estimate_power(n)
    cat("Current sample size:", n, "Power:", estimated_power, "\n")
    if (!is.na(estimated_power) && estimated_power >= target_power) {
      return(n)
    }
    n <- n + step_size
  }
}

# Estimate required sample size
required_sample_size <- find_sample_size(start_n, step_size, power)
cat("Required sample size:", required_sample_size, "\n")

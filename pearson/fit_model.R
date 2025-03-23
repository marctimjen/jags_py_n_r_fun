#!/usr/bin/env Rscript

# Load required libraries
library(rjags)
library(coda)
library(ggplot2)
library(parallel)

# Set up parallel processing
n_cores <- detectCores() - 1  # Use all cores but one
if(n_cores < 1) n_cores <- 1
cat("Using", n_cores, "cores for parallel processing\n")

# Set MCMC parameters
n_chains <- 40  # Use either 4 chains or number of cores, whichever is smaller
n_adapt <- 2000      # Adaptation iterations
n_burnin <- 5000     # Burn-in iterations
n_iter <- 10000      # Post burn-in iterations
n_thin <- 5          # Thinning rate

# Read the generated data
cat("Loading generated data...\n")
data <- read.csv("pearson/generated_data.csv")
true_params <- read.csv("pearson/true_parameters.csv")

# Extract the observed data (x variables)
x_obs <- cbind(data$x1, data$x2)
n_obs <- nrow(x_obs)

# Known measurement error standard deviations (in practice, these might be estimated or come from calibration data)
error_sd1 <- 0.5  
error_sd2 <- 0.3  

# Handle locale warnings
Sys.setlocale("LC_ALL", "C")

# Initialize the JAGS model using corr_2.txt
cat("Preparing data for model fitting...\n")
jags_data <- list(
  x = x_obs,
  n = n_obs,
  lambdaerror = c(1/(error_sd1^2), 1/(error_sd2^2))  # Precision of measurement error
)

# Initialize parameter values
init_values <- function() {
  list(
    mu = c(mean(x_obs[,1]), mean(x_obs[,2])),
    lambda = c(1/var(x_obs[,1]), 1/var(x_obs[,2])),
    r = cor(x_obs[,1], x_obs[,2]),
    y = matrix(0, nrow=n_obs, ncol=2)  # Initialize y matrix with zeros
  )
}

# Check if model file exists
if(!file.exists("pearson/corr_2.txt")) {
  cat("ERROR: Model file not found at 'pearson/corr_2.txt'\n")
  quit(status = 1)
}

cat("Setting up parallel processing with", n_chains, "chains across", n_cores, "cores...\n")

# Define a function to run a single chain
run_chain <- function(chain_id) {
  # Load required libraries within the function for parallel execution
  library(rjags)
  library(coda)
  
  # Set seed based on chain ID for reproducibility
  set.seed(123 + chain_id)
  
  # Create initialization for this chain
  init <- init_values()
  
  # Compile model
  cat("Chain", chain_id, ": Compiling model...\n")
  jags_model <- jags.model(
    file = "pearson/corr_2.txt",
    data = jags_data,
    inits = init,
    n.chains = 1,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  
  # Burn-in phase
  cat("Chain", chain_id, ": Burn-in phase...\n")
  update(jags_model, n.iter = n_burnin, progress.bar = "none")
  
  # Sampling phase
  cat("Chain", chain_id, ": Sampling phase...\n")
  samples <- coda.samples(
    model = jags_model,
    variable.names = c("mu", "sigma", "r"),
    n.iter = n_iter,
    thin = n_thin,
    progress.bar = "none"
  )
  
  cat("Chain", chain_id, ": Completed\n")
  return(samples[[1]])  # Return the MCMC samples
}

# Execute chains in parallel
cat("Starting parallel execution of", n_chains, "chains...\n")
start_time <- Sys.time()

# Use mclapply for parallel execution (better on Linux)
chain_results <- mclapply(1:n_chains, run_chain, mc.cores = n_cores)
cat("All chains completed\n")

# Calculate runtime
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("Total runtime:", round(runtime, 2), "minutes\n")

# Convert results to mcmc.list object
jags_mcmc <- as.mcmc.list(lapply(chain_results, function(x) as.mcmc(x)))

# Calculate summary statistics
cat("Calculating summary statistics...\n")
jags_summary <- summary(jags_mcmc)
print(jags_summary$statistics)
print(jags_summary$quantiles)

# Combine MCMC chains
all_samples <- do.call(rbind, jags_mcmc)

# Extract parameter samples
r_samples <- all_samples[, "r"]
mu1_samples <- all_samples[, "mu[1]"]
mu2_samples <- all_samples[, "mu[2]"]
sigma1_samples <- all_samples[, "sigma[1]"]
sigma2_samples <- all_samples[, "sigma[2]"]

# Compare with true values
comparison <- data.frame(
  Parameter = c("r", "mu1", "mu2", "sigma1", "sigma2"),
  True = c(
    true_params$true_value[true_params$param == "r"],
    true_params$true_value[true_params$param == "mu1"],
    true_params$true_value[true_params$param == "mu2"],
    true_params$true_value[true_params$param == "sigma1"],
    true_params$true_value[true_params$param == "sigma2"]
  ),
  Estimated = c(
    mean(r_samples),
    mean(mu1_samples),
    mean(mu2_samples),
    mean(sigma1_samples),
    mean(sigma2_samples)
  )
)
comparison$Difference <- comparison$Estimated - comparison$True
comparison$StdErr <- c(
  sd(r_samples),
  sd(mu1_samples),
  sd(mu2_samples),
  sd(sigma1_samples),
  sd(sigma2_samples)
)

# Print the comparison
cat("\nComparison of true values with posterior estimates:\n")
print(comparison)

# Create posterior distribution plots
# 1. Correlation
p_r <- ggplot(data.frame(r = r_samples), aes(x = r)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = comparison$True[1], color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = comparison$Estimated[1], color = "blue", linetype = "solid", size = 1) +
  labs(
    title = "Posterior Distribution of Correlation (r)",
    subtitle = paste("True:", round(comparison$True[1], 3), "  Estimated:", round(comparison$Estimated[1], 3)),
    x = "Correlation (r)",
    y = "Frequency"
  ) +
  theme_minimal()

# 2. Means
p_mu <- ggplot() +
  geom_histogram(data = data.frame(mu = mu1_samples), aes(x = mu, fill = "mu1"), bins = 30, alpha = 0.5) +
  geom_histogram(data = data.frame(mu = mu2_samples), aes(x = mu, fill = "mu2"), bins = 30, alpha = 0.5) +
  geom_vline(xintercept = comparison$True[2], color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = comparison$True[3], color = "darkred", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("mu1" = "blue", "mu2" = "green"), 
                     labels = c("mu1" = "Mean 1", "mu2" = "Mean 2"),
                     name = "Parameter") +
  labs(
    title = "Posterior Distributions of Means",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal()

# 3. Sigmas
p_sigma <- ggplot() +
  geom_histogram(data = data.frame(sigma = sigma1_samples), aes(x = sigma, fill = "sigma1"), bins = 30, alpha = 0.5) +
  geom_histogram(data = data.frame(sigma = sigma2_samples), aes(x = sigma, fill = "sigma2"), bins = 30, alpha = 0.5) +
  geom_vline(xintercept = comparison$True[4], color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = comparison$True[5], color = "darkred", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("sigma1" = "purple", "sigma2" = "orange"), 
                     labels = c("sigma1" = "Std Dev 1", "sigma2" = "Std Dev 2"),
                     name = "Parameter") +
  labs(
    title = "Posterior Distributions of Standard Deviations",
    x = "Value",
    y = "Frequency"
  ) +
  theme_minimal()

# Save the plots
png("pearson/posterior_r.png", width = 800, height = 600, res = 100)
print(p_r)
dev.off()

png("pearson/posterior_means.png", width = 800, height = 600, res = 100)
print(p_mu)
dev.off()

png("pearson/posterior_sigmas.png", width = 800, height = 600, res = 100)
print(p_sigma)
dev.off()

# Save the MCMC samples for future reference
saveRDS(jags_mcmc, file = "pearson/mcmc_samples.rds")

# Save the parameter comparison
write.csv(comparison, "pearson/parameter_comparison.csv", row.names = FALSE)

# Create trace plots to assess convergence
png("pearson/trace_plots.png", width = 1000, height = 800, res = 100)
plot(jags_mcmc)
dev.off()

# Gelman-Rubin diagnostic
gelman_diag <- gelman.diag(jags_mcmc)
cat("\nGelman-Rubin diagnostics for convergence:\n")
print(gelman_diag)

# Effective sample size
eff_size <- effectiveSize(jags_mcmc)
cat("\nEffective sample sizes:\n")
print(eff_size)

# Performance metrics
performance <- data.frame(
  Metric = c("Runtime_mins", "Min_ESS", "Max_PSRF"),
  Value = c(
    round(as.numeric(runtime), 2),
    min(eff_size),
    max(gelman_diag$psrf[,1])
  )
)
write.csv(performance, "pearson/performance_metrics.csv", row.names = FALSE)

cat("\nModel fitting completed in", round(runtime, 2), "minutes\n")
cat("Results saved to the pearson directory.\n") 
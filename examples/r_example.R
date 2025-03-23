library(rjags)
library(coda)
library(ggplot2)
library(vscDebugger)  # Load the debugger package

# Define the JAGS model
model_string <- "
model {
    # Priors
    mu ~ dnorm(0, 0.001)
    tau ~ dgamma(0.001, 0.001)
    sigma <- 1/sqrt(tau)
    
    # Likelihood
    for (i in 1:n) {
        y[i] ~ dnorm(mu, tau)
    }
}
"

# Generate synthetic data
set.seed(42)
true_mu <- 5
true_sigma <- 2

n_obs <- 50
data <- rnorm(n_obs, true_mu, true_sigma)

# DEBUG: Print data summary
print("Data Summary:")
print(summary(data))
print(paste("Data length:", length(data)))

# Set up the model
jags_data <- list(
  y = data,
  n = n_obs
)

# DEBUG: Check jags_data
print("JAGS data structure:")
str(jags_data)

# Write the model to a file
model_file <- "model.jags"
cat(model_string, file = model_file)

# DEBUG: Check if file was created
if(file.exists(model_file)) {
  print(paste("Model file created at:", getwd(), "/", model_file))
} else {
  stop("Failed to create model file")
}

# Compile the model
print("Compiling JAGS model...")
tryCatch({
  jags_model <- jags.model(
    file = model_file,
    data = jags_data,
    n.chains = 4,
    n.adapt = 1000
  )
  print("Model compiled successfully")
}, error = function(e) {
  print(paste("Error in model compilation:", e))
  stop(e)
})

# Perform sampling
print("Beginning MCMC sampling...")
tryCatch({
  jags_samples <- coda.samples(
    model = jags_model,
    variable.names = c("mu", "sigma"),
    n.iter = 2000
  )
  print("Sampling completed")
}, error = function(e) {
  print(paste("Error in sampling:", e))
  stop(e)
})

# Examine the results
summary_stats <- summary(jags_samples)
print("MCMC Summary Statistics:")
print(summary_stats)

# Combine chains for plotting
all_samples <- do.call(rbind, jags_samples)
samples_df <- data.frame(
  mu = all_samples[, "mu"],
  sigma = all_samples[, "sigma"]
)

# Print summary statistics
cat("True mu:", true_mu, "\n")
cat("Estimated mu:", mean(samples_df$mu), "std:", sd(samples_df$mu), "\n")
cat("True sigma:", true_sigma, "\n")
cat("Estimated sigma:", mean(samples_df$sigma), "std:", sd(samples_df$sigma), "\n")

# DEBUG: Check convergence diagnostics
print("Gelman-Rubin diagnostics:")
print(gelman.diag(jags_samples))
print("Effective sample sizes:")
print(effectiveSize(jags_samples))

# Plot the posterior distributions
p1 <- ggplot(samples_df, aes(x = mu)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "skyblue") +
  geom_vline(xintercept = true_mu, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Posterior distribution of mu", x = "mu", y = "Count") +
  theme_minimal()

p2 <- ggplot(samples_df, aes(x = sigma)) +
  geom_histogram(bins = 30, alpha = 0.7, fill = "skyblue") +
  geom_vline(xintercept = true_sigma, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Posterior distribution of sigma", x = "sigma", y = "Count") +
  theme_minimal()

# Save plots
png("r_jags_posterior_mu.png", width = 800, height = 600)
print(p1)
dev.off()

png("r_jags_posterior_sigma.png", width = 800, height = 600)
print(p2)
dev.off()

# Clean up
file.remove(model_file) 
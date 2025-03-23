#!/usr/bin/env Rscript

# Load required libraries
library(rjags)
library(coda)
library(ggplot2)

set.seed(123) # For reproducibility

# Set true parameters for data generation
true_r <- 0.7          # True correlation
true_mu1 <- 10         # True mean of variable 1
true_mu2 <- 5          # True mean of variable 2
true_sigma1 <- 2       # True standard deviation of variable 1
true_sigma2 <- 1.5     # True standard deviation of variable 2
error_sd1 <- 0.5       # Measurement error for variable 1
error_sd2 <- 0.3       # Measurement error for variable 2
n_obs <- 20           # Number of observations

# Create directory for output
data_dir <- "pearson/generated_data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory created:", data_dir, "\n")

# Handle locale issues
Sys.setlocale("LC_ALL", "C")

# Prepare data for JAGS
jags_data <- list(
  n = n_obs,
  lambdaerror = c(1/(error_sd1^2), 1/(error_sd2^2)),
  mu = c(true_mu1, true_mu2),
  lambda = c(1/(true_sigma1^2), 1/(true_sigma2^2)),
  r = true_r
)


# Compile the model
cat("Compiling model for data generation...\n")
jags_model <- jags.model(
  file = "pearson/corr_2.txt", # Use the data generation model
  data = jags_data,
  n.chains = 1,
  n.adapt = 1000
)

# Generate variables to track
cat("Generating simulated data...\n")
samples <- coda.samples(
  model = jags_model,
  variable.names = c("x", "y"),
  n.iter = 1
)

# Extract the data from JAGS output
x_obs <- matrix(0, nrow = n_obs, ncol = 2)
y_true <- matrix(0, nrow = n_obs, ncol = 2)

for (i in 1:n_obs) {
  for (j in 1:2) {
    x_obs[i, j] <- samples[[1]][1, paste0("x[", i, ",", j, "]")]
    y_true[i, j] <- samples[[1]][1, paste0("y[", i, ",", j, "]")]
  }
}

# Convert to data frame
data_df <- data.frame(
  x1 = x_obs[, 1],
  x2 = x_obs[, 2],
  y1 = y_true[, 1],
  y2 = y_true[, 2]
)

# Calculate the observed correlation
obs_cor <- cor(x_obs[, 1], x_obs[, 2])
true_cor <- cor(y_true[, 1], y_true[, 2])

# Create and save a scatter plot
p <- ggplot(data_df, aes(x = x1, y = x2)) +
  geom_point(alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", color = "red", fill = "#FF000080") +
  labs(
    title = "Generated Data with Measurement Error",
    subtitle = paste("True r =", round(true_r, 3), 
                     ", Empirical true r =", round(true_cor, 3),
                     ", Observed r =", round(obs_cor, 3)),
    x = "Variable 1",
    y = "Variable 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the plot
png(file.path("pearson/generated_data", "generated_data_plot.png"), width = 800, height = 600, res = 100)
print(p)
dev.off()

# Save the data
write.csv(data_df, file.path("pearson/generated_data", "generated_data.csv"), row.names = FALSE)

# Save the true parameters for later comparison
true_params <- data.frame(
  param = c("mu1", "mu2", "sigma1", "sigma2", "r"),
  true_value = c(true_mu1, true_mu2, true_sigma1, true_sigma2, true_r)
)
write.csv(true_params, file.path("pearson/generated_data", "true_parameters.csv"), row.names = FALSE)

cat("Data generation completed.\n")
cat("Files saved in pearson/generated_data:\n")
cat("- generated_data.csv (the data for analysis)\n")
cat("- true_parameters.csv (true parameter values)\n")
cat("- generated_data_plot.png (scatter plot of generated data)\n") 
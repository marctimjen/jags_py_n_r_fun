#!/usr/bin/env Rscript

# Load required libraries
library(R2jags)
library(coda)
library(ggplot2)
library(parallel)

# Set MCMC parameters
chains_total <<- 40        # Number of MCMC chains
burn_baby <<- 5000     # Burn-in iterations
iterations <<- 10000      # Post burn-in iterations
thinning <<- 5          # Thinning rate

# Set up parallel processing
n_cores <- detectCores() - 1  # Use all cores but one
if(n_cores < 1) n_cores <- 1
cat("Running JAGS with", chains_total, "chains on", n_cores, "cores\n")

# Handle locale warnings
Sys.setlocale("LC_ALL", "C")

# Read the generated data
cat("Loading generated data...\n")
data <- read.csv("pearson/generated_data/generated_data.csv")
true_params <- read.csv("pearson/generated_data/true_parameters.csv")

# Extract the observed data (x variables)
x_obs <- cbind(data$x1, data$x2)
n_obs <- nrow(x_obs)

# Known measurement error standard deviations
error_sd1 <- 0.5  
error_sd2 <- 0.3  

# Prepare data for JAGS
cat("Preparing data for model fitting...\n")
jags_data <- list(
  x = x_obs,
  n = n_obs
)

# Parameters to monitor
params <- c("mu", "sigma", "r")

# Check if model file exists
if(!file.exists("pearson/corr_2.txt")) {
  cat("ERROR: Model file not found at 'pearson/corr_2.txt'\n")
  quit(status = 1)
}

# Set global seed for reproducibility
set.seed(123)

# Start timing
start_time <- Sys.time()
cat("Starting JAGS with parallel processing at", format(start_time), "\n")

# Run JAGS using jags.parallel from R2jags
cat("Running JAGS with parallel approach...\n")

# Run the model with parallel processing
jags_fit <- jags.parallel(
  data = jags_data,
  parameters.to.save = params,
  model.file = "pearson/corr_1.txt",
  n.chains = chains_total,
  n.iter = 10000,  # Hard-coded: n_burnin + n_iter (5000 + 10000)
  n.burnin = 5000, # Hard-coded
  n.thin = 5,
  DIC = TRUE,
  jags.seed = 123,
  n.cluster = n_cores
)

cat("Parallel JAGS completed successfully\n")

# Calculate runtime
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("JAGS completed in", round(runtime, 2), "minutes\n")

# Extract MCMC samples
cat("Processing results...\n")

# Fix for jags.parallel output format
if (inherits(jags_fit, "jags")) {
  # Standard JAGS output
  jags_mcmc <- as.mcmc.list(jags_fit$BUGSoutput$sims.array)
  cat("Standard JAGS format detected\n")
} else {
  # jags.parallel has a different structure
  # Convert the results into a proper mcmc.list object
  cat("Parallel JAGS format detected\n")
  
  # Debug the structure of jags_fit
  cat("Structure of jags_fit:\n")
  print(names(jags_fit))
  cat("Structure of BUGSoutput:\n")
  print(names(jags_fit$BUGSoutput))
  
  # Check if sims.array exists in the output
  if (!is.null(jags_fit$BUGSoutput$sims.array)) {
    sims <- jags_fit$BUGSoutput$sims.array
    cat("sims.array dimensions:", dim(sims), "\n")
    
    # If sims.array has proper dimensions, create mcmc.list
    if (length(dim(sims)) == 3) {
      chains <- list()
      for (i in 1:dim(sims)[2]) {
        # Create an mcmc object for each chain
        chains[[i]] <- coda::mcmc(sims[, i, ])
      }
      jags_mcmc <- coda::mcmc.list(chains)
      cat("Converted parallel JAGS output to mcmc.list format with", length(chains), "chains\n")
    } else {
      # If dimensions are wrong, try alternate approach
      cat("Unexpected sims.array dimensions, trying alternate approach\n")
      
      # Check if sims.list exists as an alternative
      if (!is.null(jags_fit$BUGSoutput$sims.list)) {
        # Convert sims.list to mcmc.list
        sims_list <- jags_fit$BUGSoutput$sims.list
        cat("Using sims.list instead. Available parameters:", names(sims_list), "\n")
        
        # Create a data frame from sims.list
        sims_df <- as.data.frame(sims_list)
        
        # Create a single mcmc object (no chain information)
        jags_mcmc <- coda::mcmc.list(coda::mcmc(sims_df))
        cat("Created mcmc.list from sims.list with 1 chain\n")
      } else {
        cat("ERROR: Could not find suitable MCMC samples in the output\n")
        # Create empty mcmc list as fallback
        empty_df <- data.frame(r = NA)
        jags_mcmc <- coda::mcmc.list(coda::mcmc(empty_df))
      }
    }
  } else if (!is.null(jags_fit$BUGSoutput$sims.list)) {
    # If sims.array doesn't exist but sims.list does
    sims_list <- jags_fit$BUGSoutput$sims.list
    cat("Using sims.list. Available parameters:", names(sims_list), "\n")
    
    # Create a data frame from sims.list
    sims_df <- as.data.frame(sims_list)
    
    # Create a single mcmc object (no chain information)
    jags_mcmc <- coda::mcmc.list(coda::mcmc(sims_df))
    cat("Created mcmc.list from sims.list with 1 chain\n")
  } else {
    cat("ERROR: Could not find either sims.array or sims.list in the output\n")
    # Create empty mcmc list as fallback
    empty_df <- data.frame(r = NA)
    jags_mcmc <- coda::mcmc.list(coda::mcmc(empty_df))
  }
}

# Get DIC
dic_value <- jags_fit$BUGSoutput$DIC
cat("DIC:", dic_value, "\n")

# Print summary statistics
cat("Summary statistics:\n")
jags_summary <- summary(jags_mcmc)
print(jags_summary$statistics)
print(jags_summary$quantiles)

# Combine MCMC chains
all_samples <- do.call(rbind, jags_mcmc)

# Debug parameter names
cat("\nParameter names in MCMC samples:\n")
param_names <- colnames(all_samples)
print(param_names)

# Extract parameter samples - with safer approach to handle different naming patterns
# Function to find the right parameter name
find_param <- function(pattern, param_names) {
  matches <- grep(pattern, param_names, fixed = FALSE)
  if (length(matches) > 0) {
    return(param_names[matches[1]])
  } else {
    warning(paste("Parameter matching", pattern, "not found"))
    return(NULL)
  }
}

# Find parameter names that match expected patterns
r_param <- find_param("^r$", param_names)
mu1_param <- find_param("mu.1.", param_names)  # Try alternative format
if (is.null(mu1_param)) mu1_param <- find_param("mu.*1", param_names)
mu2_param <- find_param("mu.2.", param_names)  # Try alternative format
if (is.null(mu2_param)) mu2_param <- find_param("mu.*2", param_names)
sigma1_param <- find_param("sigma.1.", param_names)  # Try alternative format
if (is.null(sigma1_param)) sigma1_param <- find_param("sigma.*1", param_names)
sigma2_param <- find_param("sigma.2.", param_names)  # Try alternative format
if (is.null(sigma2_param)) sigma2_param <- find_param("sigma.*2", param_names)

cat("\nMatched parameter names:\n")
cat("r parameter:", r_param, "\n")
cat("mu1 parameter:", mu1_param, "\n")
cat("mu2 parameter:", mu2_param, "\n")
cat("sigma1 parameter:", sigma1_param, "\n")
cat("sigma2 parameter:", sigma2_param, "\n")

# Safely extract samples
r_samples <- if (!is.null(r_param)) all_samples[, r_param] else rep(NA, nrow(all_samples))
mu1_samples <- if (!is.null(mu1_param)) all_samples[, mu1_param] else rep(NA, nrow(all_samples))
mu2_samples <- if (!is.null(mu2_param)) all_samples[, mu2_param] else rep(NA, nrow(all_samples))
sigma1_samples <- if (!is.null(sigma1_param)) all_samples[, sigma1_param] else rep(NA, nrow(all_samples))
sigma2_samples <- if (!is.null(sigma2_param)) all_samples[, sigma2_param] else rep(NA, nrow(all_samples))

# Compare with true values
comparison <- data.frame(
  Parameter = c("r", "mu1", "mu2", "sigma1", "sigma2"),
  True = c(
    true_params$true_value[true_params$param == "r"],
    true_params$true_value[true_params$param == "mu1"],
    true_params$true_value[true_params$param == "mu2"],
    true_params$true_value[true_params$param == "sigma1"],
    true_params$true_value[true_params$param == "sigma2"]
  )
)
comparison$Estimated <- c(
  if (all(is.na(r_samples))) NA else mean(r_samples, na.rm = TRUE),
  if (all(is.na(mu1_samples))) NA else mean(mu1_samples, na.rm = TRUE),
  if (all(is.na(mu2_samples))) NA else mean(mu2_samples, na.rm = TRUE),
  if (all(is.na(sigma1_samples))) NA else mean(sigma1_samples, na.rm = TRUE),
  if (all(is.na(sigma2_samples))) NA else mean(sigma2_samples, na.rm = TRUE)
)
comparison$Difference <- comparison$Estimated - comparison$True
comparison$StdErr <- c(
  if (all(is.na(r_samples))) NA else sd(r_samples, na.rm = TRUE),
  if (all(is.na(mu1_samples))) NA else sd(mu1_samples, na.rm = TRUE),
  if (all(is.na(mu2_samples))) NA else sd(mu2_samples, na.rm = TRUE),
  if (all(is.na(sigma1_samples))) NA else sd(sigma1_samples, na.rm = TRUE),
  if (all(is.na(sigma2_samples))) NA else sd(sigma2_samples, na.rm = TRUE)
)

# Print the comparison
cat("\nComparison of true values with posterior estimates:\n")
print(comparison)

# Create posterior distribution plots - only for available parameters
# 1. Correlation (if available)
if (!all(is.na(r_samples))) {
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
    
  # Save the plot
  png(file.path("pearson/posterior_plots/posterior_r.png"), width = 800, height = 600, res = 100)
  print(p_r)
  dev.off()
} else {
  cat("Skipping r plot - parameter not found\n")
}

# 2. Means (if available)
mu_plots_available <- !all(is.na(mu1_samples)) || !all(is.na(mu2_samples))
if (mu_plots_available) {
  p_data <- data.frame()
  if (!all(is.na(mu1_samples))) {
    p_data <- rbind(p_data, data.frame(mu = mu1_samples, type = "mu1"))
  }
  if (!all(is.na(mu2_samples))) {
    p_data <- rbind(p_data, data.frame(mu = mu2_samples, type = "mu2"))
  }
  
  if (nrow(p_data) > 0) {
    p_mu <- ggplot(p_data, aes(x = mu, fill = type)) +
      geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
      scale_fill_manual(values = c("mu1" = "blue", "mu2" = "green"), 
                         labels = c("mu1" = "Mean 1", "mu2" = "Mean 2"),
                         name = "Parameter") +
      labs(
        title = "Posterior Distributions of Means",
        x = "Value",
        y = "Frequency"
      ) +
      theme_minimal()
      
    # Add vertical lines for true values if they exist
    if (!all(is.na(mu1_samples)) && !is.na(comparison$True[2])) {
      p_mu <- p_mu + geom_vline(xintercept = comparison$True[2], color = "red", linetype = "dashed", size = 1)
    }
    if (!all(is.na(mu2_samples)) && !is.na(comparison$True[3])) {
      p_mu <- p_mu + geom_vline(xintercept = comparison$True[3], color = "darkred", linetype = "dashed", size = 1)
    }
    
    # Save the plot
    png(file.path("pearson/posterior_plots/posterior_means.png"), width = 800, height = 600, res = 100)
    print(p_mu)
    dev.off()
  }
} else {
  cat("Skipping means plot - parameters not found\n")
}

# 3. Sigmas (if available)
sigma_plots_available <- !all(is.na(sigma1_samples)) || !all(is.na(sigma2_samples))
if (sigma_plots_available) {
  p_data <- data.frame()
  if (!all(is.na(sigma1_samples))) {
    p_data <- rbind(p_data, data.frame(sigma = sigma1_samples, type = "sigma1"))
  }
  if (!all(is.na(sigma2_samples))) {
    p_data <- rbind(p_data, data.frame(sigma = sigma2_samples, type = "sigma2"))
  }
  
  if (nrow(p_data) > 0) {
    p_sigma <- ggplot(p_data, aes(x = sigma, fill = type)) +
      geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
      scale_fill_manual(values = c("sigma1" = "purple", "sigma2" = "orange"), 
                         labels = c("sigma1" = "Std Dev 1", "sigma2" = "Std Dev 2"),
                         name = "Parameter") +
      labs(
        title = "Posterior Distributions of Standard Deviations",
        x = "Value",
        y = "Frequency"
      ) +
      theme_minimal()
      
    # Add vertical lines for true values if they exist
    if (!all(is.na(sigma1_samples)) && !is.na(comparison$True[4])) {
      p_sigma <- p_sigma + geom_vline(xintercept = comparison$True[4], color = "red", linetype = "dashed", size = 1)
    }
    if (!all(is.na(sigma2_samples)) && !is.na(comparison$True[5])) {
      p_sigma <- p_sigma + geom_vline(xintercept = comparison$True[5], color = "darkred", linetype = "dashed", size = 1)
    }
    
    # Save the plot
    png(file.path("pearson/posterior_plots/posterior_sigmas.png"), width = 800, height = 600, res = 100)
    print(p_sigma)
    dev.off()
  }
} else {
  cat("Skipping sigmas plot - parameters not found\n")
}

# Save the MCMC samples
saveRDS(jags_mcmc, file = file.path("pearson/posterior_plots/mcmc_samples.rds"))

# Save the parameter comparison
write.csv(comparison, file.path("pearson/posterior_plots/parameter_comparison.csv"), row.names = FALSE)

# Create trace plots to assess convergence
png(file.path("pearson/posterior_plots/trace_plots.png"), width = 1000, height = 800, res = 100)
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
  Metric = c("Runtime_mins", "DIC", "Min_ESS", "Max_PSRF"),
  Value = c(
    round(as.numeric(runtime), 2),
    round(dic_value, 2),
    min(eff_size),
    max(gelman_diag$psrf[,1])
  )
)
write.csv(performance, file.path("pearson/posterior_plots/performance_metrics.csv"), row.names = FALSE)

# Create additional trace plots for specific parameters
# Trace plot for mu parameters
if (!is.null(mu1_param) || !is.null(mu2_param)) {
  # Filter for mu parameters only
  mu_params <- grep("mu", param_names, value = TRUE)
  if (length(mu_params) > 0) {
    cat("Creating trace plots for mu parameters:", paste(mu_params, collapse=", "), "\n")
    
    # Save mu trace plots - improved version
    png(file.path("pearson/posterior_plots/trace_plots_mu.png"), width = 1200, height = 800, res = 100)
    
    # Extract mu samples from all chains
    if (length(jags_mcmc) > 0 && nrow(all_samples) > 0) {
      # Setup plot layout based on number of parameters
      par(mfrow = c(length(mu_params), 1))
      
      # For each mu parameter
      for (param in mu_params) {
        # Plot all chains for this parameter with different colors
        if (param %in% colnames(all_samples)) {
          # Extract data for each chain
          chain_data <- list()
          chain_colors <- rainbow(length(jags_mcmc))
          
          # Prepare data from each chain
          for (i in 1:length(jags_mcmc)) {
            if (param %in% colnames(jags_mcmc[[i]])) {
              chain_data[[i]] <- jags_mcmc[[i]][, param]
            }
          }
          
          # Find y-axis limits
          y_min <- min(sapply(chain_data, function(x) if(length(x) > 0) min(x) else NA), na.rm=TRUE)
          y_max <- max(sapply(chain_data, function(x) if(length(x) > 0) max(x) else NA), na.rm=TRUE)
          
          # Create empty plot first
          plot(NA, xlim=c(1, max(sapply(chain_data, length))), 
               ylim=c(y_min, y_max),
               main=paste("Trace plot for", param),
               xlab="Iteration", ylab="Value")
          
          # Add each chain's data
          for (i in 1:length(chain_data)) {
            if (length(chain_data[[i]]) > 0) {
              lines(chain_data[[i]], col=chain_colors[i], lwd=1)
            }
          }
          
          # Add legend
          if (length(chain_data) > 1) {
            legend("topright", legend=paste("Chain", 1:length(chain_data)), 
                   col=chain_colors, lty=1, cex=0.8)
          }
        }
      }
    } else {
      # Fallback if we can't extract parameter data
      plot(1:10, 1:10, type="n", xlab="", ylab="", 
           main="Cannot create mu trace plots with current MCMC format")
      text(5, 5, "Parameter extraction failed")
    }
    
    dev.off()
    cat("Trace plots for mu parameters saved to pearson/posterior_plots/trace_plots_mu.png\n")
  } else {
    cat("No mu parameters found for trace plots\n")
  }
}

cat("\nR2jags model fitting completed in", round(runtime, 2), "minutes\n")
cat("Results saved to pearson/posterior_plots directory\n") 
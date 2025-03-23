#!/usr/bin/env Rscript

# Load required libraries
library(R2jags)
library(coda)
library(ggplot2)
library(parallel)

# Set MCMC parameters
n_chains <- 4            # Number of MCMC chains
n_burnin <- 5000         # Burn-in iterations
n_iter <- 10000          # Post burn-in iterations
n_thin <- 5              # Thinning rate

# Set up parallel processing
n_cores <- detectCores() - 1  # Use all cores but one
if(n_cores < 1) n_cores <- 1
cat("Running JAGS with", n_chains, "chains on", n_cores, "cores\n")

# Enable multi-core processing for MCMC
options(mc.cores = n_cores)

# Handle locale warnings
Sys.setlocale("LC_ALL", "C")

# Set paths
model_file <- "kappa/kappa.txt"
data_dir <- "kappa/generated_data"
output_dir <- "kappa/posterior_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Explicit check for kappa.txt
cat("Using model file:", model_file, "\n")
if(!file.exists(model_file)) {
  cat("ERROR: Kappa model file not found at", model_file, "\n")
  cat("Please make sure kappa.txt exists in the kappa directory.\n")
  quit(status = 1)
}

# Load the data
cat("Loading generated data...\n")
if(file.exists(file.path(data_dir, "jags_data.rds"))) {
  jags_data <- readRDS(file.path(data_dir, "jags_data.rds"))
  cat("Loaded data from jags_data.rds\n")
  cat("Data summary: n =", jags_data$n, "with counts:", jags_data$y, "\n")
} else if(file.exists(file.path(data_dir, "generated_data.csv"))) {
  # Load from CSV and convert to JAGS format
  data_df <- read.csv(file.path(data_dir, "generated_data.csv"))
  jags_data <- list(
    y = data_df$count,
    n = sum(data_df$count)
  )
  cat("Loaded and converted data from generated_data.csv\n")
  cat("Data summary: n =", jags_data$n, "with counts:", jags_data$y, "\n")
} else {
  cat("ERROR: No data found in", data_dir, "\n")
  quit(status = 1)
}

# Load true parameters
true_params <- NULL
if(file.exists(file.path(data_dir, "true_parameters.csv"))) {
  true_params <- read.csv(file.path(data_dir, "true_parameters.csv"))
  cat("Loaded true parameters for comparison\n")
} else {
  cat("Note: No true parameters found for comparison\n")
}

# Parameters to monitor
params <- c("alpha", "beta", "gamma", "xi", "psi", "kappa")

# Set global seed for reproducibility
set.seed(123)

# Start timing
start_time <- Sys.time()
cat("Starting JAGS with parallel processing at", format(start_time), "\n")

# Run JAGS using jags.parallel from R2jags
cat("Running JAGS with R2jags parallel approach...\n")

# Run the model with parallel processing
jags_fit <- jags.parallel(
  data = jags_data,
  parameters.to.save = params,
  model.file = model_file,
  n.chains = 40,
  n.iter = 10000,  # Total iterations (burn-in + sampling)
  n.burnin = 5000,
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
        empty_df <- data.frame(alpha = NA)
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
    empty_df <- data.frame(alpha = NA)
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

# Extract individual parameter samples
extract_param <- function(param_name) {
  if(param_name %in% param_names) {
    return(all_samples[, param_name])
  } else {
    warning(paste("Parameter", param_name, "not found in samples"))
    return(rep(NA, nrow(all_samples)))
  }
}

alpha_samples <- extract_param("alpha")
beta_samples <- extract_param("beta")
gamma_samples <- extract_param("gamma")
kappa_samples <- extract_param("kappa")

# Compare with true values if available
if(!is.null(true_params)) {
  get_true_value <- function(param_name) {
    idx <- which(true_params$param == param_name)
    if(length(idx) > 0) {
      return(true_params$true_value[idx])
    } else {
      return(NA)
    }
  }
  
  comparison <- data.frame(
    Parameter = c("alpha", "beta", "gamma", "kappa"),
    True = c(
      get_true_value("alpha"),
      get_true_value("beta"),
      get_true_value("gamma"),
      get_true_value("kappa")
    )
  )
  
  comparison$Estimated <- c(
    mean(alpha_samples, na.rm = TRUE),
    mean(beta_samples, na.rm = TRUE),
    mean(gamma_samples, na.rm = TRUE),
    mean(kappa_samples, na.rm = TRUE)
  )
  
  comparison$Difference <- comparison$Estimated - comparison$True
  comparison$StdErr <- c(
    sd(alpha_samples, na.rm = TRUE),
    sd(beta_samples, na.rm = TRUE),
    sd(gamma_samples, na.rm = TRUE),
    sd(kappa_samples, na.rm = TRUE)
  )
  
  # Print the comparison
  cat("\nComparison of true values with posterior estimates:\n")
  print(comparison)
  
  # Save the comparison
  write.csv(comparison, file.path(output_dir, "parameter_comparison.csv"), row.names = FALSE)
}

# Create posterior plots
create_posterior_plot <- function(samples, param_name, true_value = NA) {
  if(all(is.na(samples))) {
    cat("Skipping", param_name, "plot - parameter not found\n")
    return(NULL)
  }
  
  df <- data.frame(value = samples)
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(
      title = paste("Posterior Distribution of", param_name),
      x = param_name,
      y = "Frequency"
    ) +
    theme_minimal()
  
  # Add true value line if available
  if(!is.na(true_value)) {
    p <- p + 
      geom_vline(xintercept = true_value, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = mean(samples, na.rm = TRUE), color = "blue", linetype = "solid", size = 1)
  }
  
  # Save the plot
  plot_file <- file.path(output_dir, paste0("posterior_", tolower(param_name), ".png"))
  png(plot_file, width = 800, height = 600, res = 100)
  print(p)
  dev.off()
  cat("Saved posterior plot for", param_name, "to", plot_file, "\n")
  
  return(p)
}

# Create individual parameter plots
if(!is.null(true_params)) {
  create_posterior_plot(alpha_samples, "Alpha", get_true_value("alpha"))
  create_posterior_plot(beta_samples, "Beta", get_true_value("beta"))
  create_posterior_plot(gamma_samples, "Gamma", get_true_value("gamma"))
  create_posterior_plot(kappa_samples, "Kappa", get_true_value("kappa"))
} else {
  create_posterior_plot(alpha_samples, "Alpha")
  create_posterior_plot(beta_samples, "Beta")
  create_posterior_plot(gamma_samples, "Gamma")
  create_posterior_plot(kappa_samples, "Kappa")
}

# Create trace plots for convergence assessment
png(file.path(output_dir, "trace_plots.png"), width = 1000, height = 800, res = 100)
plot(jags_mcmc)
dev.off()
cat("Saved trace plots to", file.path(output_dir, "trace_plots.png"), "\n")

# Create trace plots for specific parameters
# Trace plot for key parameters
key_params <- c("alpha", "beta", "gamma", "kappa")
for (param in key_params) {
  if (param %in% param_names) {
    # Save parameter-specific trace plots
    png(file.path(output_dir, paste0("trace_plot_", param, ".png")), width = 1200, height = 600, res = 100)
    
    # Extract this parameter from all chains
    if (length(jags_mcmc) > 0) {
      # Setup plot for just this parameter
      par(mfrow=c(1,1))
      
      # Plot all chains for this parameter with different colors
      chain_colors <- rainbow(length(jags_mcmc))
      
      # Find y-axis limits across all chains
      y_limits <- range(sapply(jags_mcmc, function(chain) {
        if (param %in% colnames(chain)) {
          return(range(chain[, param]))
        } else {
          return(c(NA, NA))
        }
      }), na.rm = TRUE)
      
      # Create empty plot first
      plot(NA, xlim=c(1, nrow(jags_mcmc[[1]])), 
           ylim=y_limits,
           main=paste("Trace plot for", param),
           xlab="Iteration", ylab="Value")
      
      # Add each chain's data
      for (i in 1:length(jags_mcmc)) {
        if (param %in% colnames(jags_mcmc[[i]])) {
          lines(jags_mcmc[[i]][, param], col=chain_colors[i], lwd=1)
        }
      }
      
      # Add legend
      legend("topright", legend=paste("Chain", 1:length(jags_mcmc)), 
             col=chain_colors, lty=1, cex=0.8)
    } else {
      # Fallback if we can't extract parameter data
      plot(1:10, 1:10, type="n", xlab="", ylab="", 
          main=paste("Cannot create trace plot for", param))
      text(5, 5, "Parameter extraction failed")
    }
    
    dev.off()
    cat("Saved trace plot for", param, "to", file.path(output_dir, paste0("trace_plot_", param, ".png")), "\n")
  }
}

# Save MCMC samples
saveRDS(jags_mcmc, file.path(output_dir, "mcmc_samples.rds"))
cat("Saved MCMC samples to", file.path(output_dir, "mcmc_samples.rds"), "\n")

# Calculate Gelman-Rubin diagnostics
if(length(jags_mcmc) > 1) {
  gelman_diag <- gelman.diag(jags_mcmc)
  cat("\nGelman-Rubin diagnostics for convergence:\n")
  print(gelman_diag)
  
  # Check for convergence issues
  psrf_max <- max(gelman_diag$psrf[, 1])
  if(psrf_max > 1.1) {
    cat("\nWARNING: Potential convergence issues detected (max PSRF =", psrf_max, ")\n")
    cat("Consider running more iterations or checking model specification\n")
  } else {
    cat("\nNo convergence issues detected (max PSRF =", psrf_max, ")\n")
  }
}

# Calculate effective sample sizes
eff_size <- effectiveSize(jags_mcmc)
cat("\nEffective sample sizes:\n")
print(eff_size)

# Performance metrics
performance <- data.frame(
  Metric = c("Runtime_mins", "DIC", "Min_ESS", "Max_PSRF"),
  Value = c(
    round(as.numeric(runtime), 2),
    round(dic_value, 2),
    round(min(eff_size), 2),
    ifelse(length(jags_mcmc) > 1, round(max(gelman_diag$psrf[, 1]), 4), NA)
  )
)
write.csv(performance, file.path(output_dir, "performance_metrics.csv"), row.names = FALSE)

# Create joint posterior plot for alpha, beta, and gamma
if(!all(is.na(alpha_samples)) && !all(is.na(beta_samples)) && !all(is.na(gamma_samples))) {
  params_df <- data.frame(
    alpha = alpha_samples,
    beta = beta_samples,
    gamma = gamma_samples
  )
  
  # Create pairs plot
  png(file.path(output_dir, "joint_posterior.png"), width = 1000, height = 1000, res = 100)
  pairs(params_df, main = "Joint Posterior Distributions", pch = ".")
  dev.off()
  cat("Saved joint posterior plot to", file.path(output_dir, "joint_posterior.png"), "\n")
}

cat("\nR2jags model fitting completed in", round(runtime, 2), "minutes\n")
cat("Results saved to", output_dir, "directory\n") 
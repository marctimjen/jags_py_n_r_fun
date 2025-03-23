#!/usr/bin/env Rscript

# Load required libraries
library(rjags)
library(coda)
library(ggplot2)

set.seed(123) # For reproducibility

# Set true parameters for data generation
true_alpha <- 0.5  # Rate objective method decides 'one'
true_beta <- 0.99   # Rate surrogate method decides 'one' when objective method decides 'one'
true_gamma <- 0.99  # Rate surrogate method decides 'zero' when objective method decides 'zero'
n_total <- 10000    # Total number of observations/trials

# Create directory for output if it doesn't exist
data_dir <- "kappa/generated_data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory created:", data_dir, "\n")

# Handle locale issues
Sys.setlocale("LC_ALL", "C")

# Prepare data for JAGS
jags_data <- list(
  n = n_total,
  alpha = true_alpha,
  beta = true_beta,
  gamma = true_gamma
)

# Compile the model for data generation
cat("Compiling JAGS model for data generation...\n")
jags_model <- jags.model(
  file = "kappa/kappa.txt",
  data = jags_data,
  n.chains = 1,
  n.adapt = 1000,
  quiet = FALSE
)

# Generate the data using the model
cat("Generating simulated data from the model...\n")
variables_to_monitor <- c("y", "pi", "xi", "psi", "kappa")
simulation <- coda.samples(
  model = jags_model,
  variable.names = variables_to_monitor,
  n.iter = 1,
  thin = 1
)

# Extract the data from JAGS output
counts <- as.numeric(simulation[[1]][1, grep("^y\\[", colnames(simulation[[1]]))])
pi_values <- as.numeric(simulation[[1]][1, grep("^pi\\[", colnames(simulation[[1]]))])
true_xi <- as.numeric(simulation[[1]][1, "xi"])
true_psi <- as.numeric(simulation[[1]][1, "psi"])
true_kappa <- as.numeric(simulation[[1]][1, "kappa"])

cat("Generated counts:", counts, "\n")
cat("Generated pi values:", pi_values, "\n")
cat("True xi:", true_xi, "\n")
cat("True psi:", true_psi, "\n")
cat("True kappa:", true_kappa, "\n")

# Create data frame for the counts
data_df <- data.frame(
  category = c("Both One", "Obj One, Surr Zero", "Obj Zero, Surr One", "Both Zero"),
  count = counts,
  prob = pi_values
)

# Create a bar plot
p1 <- ggplot(data_df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Generated Data for Kappa Agreement Model",
    subtitle = paste("Alpha =", round(true_alpha, 3),
                    ", Beta =", round(true_beta, 3),
                    ", Gamma =", round(true_gamma, 3),
                    ", Kappa =", round(true_kappa, 3)),
    x = "Category",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Create a pie chart for proportions
p2 <- ggplot(data_df, aes(x = "", y = count, fill = category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  labs(
    title = "Proportion of Categories",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

# Save the plots
png(file.path(data_dir, "counts_bar_plot.png"), width = 800, height = 600, res = 100)
print(p1)
dev.off()

png(file.path(data_dir, "counts_pie_chart.png"), width = 800, height = 600, res = 100)
print(p2)
dev.off()

# Save the data
write.csv(data_df, file.path(data_dir, "generated_data.csv"), row.names = FALSE)

# Prepare data for JAGS format (for the inference model)
jags_data_for_inference <- list(
  y = counts,
  n = n_total
)
saveRDS(jags_data_for_inference, file.path(data_dir, "jags_data.rds"))

# Save the true parameters for later comparison
true_params <- data.frame(
  param = c("alpha", "beta", "gamma", "kappa"),
  true_value = c(true_alpha, true_beta, true_gamma, true_kappa)
)
write.csv(true_params, file.path(data_dir, "true_parameters.csv"), row.names = FALSE)

cat("Data generation completed.\n")
cat("Files saved in kappa/generated_data:\n")
cat("- generated_data.csv (the count data for analysis)\n")
cat("- jags_data.rds (data in JAGS format for inference)\n")
cat("- true_parameters.csv (true parameter values)\n")
cat("- counts_bar_plot.png (bar plot of the counts)\n")
cat("- counts_pie_chart.png (pie chart of the proportions)\n")
cat("True kappa value:", round(true_kappa, 4), "\n") 
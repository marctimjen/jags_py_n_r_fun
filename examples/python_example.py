import numpy as np
import pyjags
import matplotlib.pyplot as plt

# Define the JAGS model
model_code = """
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
"""

# Generate some synthetic data
np.random.seed(42)
true_mu = 5
true_sigma = 2
n_obs = 50
data = np.random.normal(true_mu, true_sigma, n_obs)

# Set up the model
model_data = {
    'y': data,
    'n': n_obs
}

# Compile the model
model = pyjags.Model(code=model_code, data=model_data,
                    init=[{'mu': 0, 'tau': 1}],
                    chains=4, adapt=1000)

# Perform sampling
samples = model.sample(2000, vars=['mu', 'sigma'])

# Examine the results
mu_samples = samples['mu']
sigma_samples = samples['sigma']

# Print summary statistics
print("True mu:", true_mu)
print("Estimated mu: mean =", np.mean(mu_samples), "std =", np.std(mu_samples))
print("True sigma:", true_sigma)
print("Estimated sigma: mean =", np.mean(sigma_samples), "std =", np.std(sigma_samples))

# Plot the posterior distributions
plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
plt.hist(mu_samples, bins=30, alpha=0.7)
plt.axvline(true_mu, color='red', linestyle='dashed', linewidth=2)
plt.title('Posterior distribution of mu')
plt.xlabel('mu')

plt.subplot(1, 2, 2)
plt.hist(sigma_samples, bins=30, alpha=0.7)
plt.axvline(true_sigma, color='red', linestyle='dashed', linewidth=2)
plt.title('Posterior distribution of sigma')
plt.xlabel('sigma')

plt.tight_layout()
plt.savefig('python_jags_posterior.png')
plt.show() 
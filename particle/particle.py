import numpy as np
import plotly.express as px
from filterpy.kalman import KalmanFilter

a = 0.8
sig = np.sqrt(1 - a**2)
omega = 3

X1 = np.array([0.0])
X2 = np.array([0.0])
Y = np.array([X2[-1] + omega * np.random.normal(0, 1)])

for _ in range(128 - 1):
    X1 = np.append(X1, a * X1[-1] + sig * np.random.normal(0, 1))
    X2 = np.append(X2, X1 + X2[-1])
    Y = np.append(Y, X2[-1] + omega * np.random.normal(0, 1))

# Kalman Filter setup
kf = KalmanFilter(dim_x=2, dim_z=1)
kf.x = np.array([0., 0.])  # Initial state (position and velocity)
kf.F = np.array([[1, 1], [0, 1]])  # State transition matrix
kf.H = np.array([[1, 0]])  # Measurement function
kf.P *= 1000  # Covariance matrix
kf.R = omega**2  # Measurement noise
kf.Q = np.array([[0.1, 0], [0, 0.1]])  # Process noise

# Apply Kalman filter to the Y time series
filtered_Y = []
for y in Y:
    kf.predict()
    kf.update(y)
    filtered_Y.append(kf.x[0])  # Store the filtered position

# Particle Filter parameters
n_particles = 10000
resample_period = 8
particles = np.zeros((n_particles, len(Y)))
weights = np.ones(n_particles) / n_particles

particles_X1 = np.zeros((n_particles, len(Y)))
particles_X2 = np.zeros((n_particles, len(Y)))

particle_adjustment = np.zeros(n_particles)

# Initialize particles
for i in range(n_particles):
    particles[i, 0] = Y[0] + np.random.normal(0, omega)

# Run particle filter
for t in range(1, len(Y)):
    # Predict step: propagate particles
    for i in range(n_particles):
        particles_X1[i, t] = a * particles_X1[i, t-1] + sig * np.random.normal(0, 1)
        particles_X2[i, t] = particles_X1[i, t] + particles_X2[i, t-1]
        particles[i, t] = particles_X2[i, t] + omega * np.random.normal(0, 1) + particle_adjustment[i]

    # Update step
    if t % resample_period == 0:
        # Calculate weights based on likelihood
        weights = np.exp(-0.5 * np.sum((Y[t-resample_period:t] - particles[:, t-resample_period:t])**2, axis=1) / (omega**2))
        weights = weights / np.sum(weights)

        weights[weights < 1e-3] = 0  # Avoid numerical issues
        weights = weights / np.sum(weights)

        particle_adjustment = particle_adjustment + Y[t] - particles[:, t]

        # Resample particles
        # indices = np.random.choice(n_particles, size=n_particles, p=weights)
        # particles[:, t] = particles[indices, t]
        # particles[:, t] = Y[t]
        # weights = np.ones(n_particles) / n_particles

# Calculate mean of particles for visualization
particle_mean = np.sum(particles * weights[:, None], axis=0)

# Calculate weighted standard deviation of particles
particle_std = np.sqrt(np.sum(weights[:, None] * (particles - particle_mean)**2, axis=0))

# Update plotting
fig = px.scatter(x=range(len(Y)), y=Y, labels={"x": "Index", "y": "Y values"}, 
                 title="Scatter Plot of Y with Kalman and Particle Filters")
fig.add_scatter(x=list(range(len(filtered_Y))), y=np.array(filtered_Y), 
                mode='lines', name='Kalman Filter')
fig.add_scatter(x=list(range(len(particle_mean))), y=particle_mean, 
                mode='lines', name='Particle Filter')

# Add shaded region for weighted standard deviation
fig.add_scatter(x=list(range(len(particle_mean))), 
                y=particle_mean + particle_std, 
                mode='lines', line=dict(dash='dot', color='gray'), 
                name='Particle Filter +1 Std Dev')
fig.add_scatter(x=list(range(len(particle_mean))), 
                y=particle_mean - particle_std, 
                mode='lines', line=dict(dash='dot', color='gray'), 
                name='Particle Filter -1 Std Dev')

# Add 5 particle trajectories with dashed lines
for i in range(min(5, n_particles)):  # Show first 5 particles
    fig.add_scatter(x=list(range(len(Y))), y=particles[i, :], 
                    mode='lines', line=dict(dash='dash'), 
                    opacity=0.5, name=f'Particle {i}')

fig.show()

print()

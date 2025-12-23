import numpy as np


# environment/orbit parameters
mu_earth = 3.986e14             # gravitational parameter [m^3/s^2]
r_earth  = 6371                 # planet radius [km]

# satellite parameters
J = np.diag([1.0, 1.0, 1.0])    # inertia tensor [kg*m^2]
h = 7000e3                      # altitude [m]
s = 7.5e3                       # speed [m/s]
r0 = np.array([h, 0, 0])        # initial position vector [m]
v0 = np.array([0, s, 0])        # initial velocity vector [m/s]

# attitude controller parameters
k_p = 0.025
k_i = k_p*0.01
k_d = 0.1

# numerical integration
method = "RK4"
dt     = 1.0                    # time step [s]
N      = 1*3600                 # number of steps
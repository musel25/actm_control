import numpy as np
import matplotlib.pyplot as plt

# Model Parameters
N = 20               # Number of sections
K = 360 * 2          # Time horizon (2 hours in simulation steps)
T = 10 / 3600        # Sampling time interval (10 s in hours)

# On-ramp indicator: Ir[i]=1 if there is an on-ramp in section i
Ir = np.zeros(N)
Ir[4] = 1   # Section 5 (MATLAB index 5 -> Python index 4)
Ir[9] = 1   # Section 10
Ir[11] = 1  # Section 12

# Section lengths (in km; 0.5 km = 500 m)
L = np.full(N, 0.5)

# Maximum densities
rho_max = np.full(N, 300)

# Free-flow speeds and congested wave speeds
vf = np.full(N, 100)
w = np.full(N, 27)

# Initial Conditions
rho_iniz = np.full(N, 60)
fi_iniz = rho_iniz * 100  # initial interface flows

# Boundary conditions for the first section (all time steps)
rho_sez0 = np.full(K, 60)
fi_sez0 = rho_sez0 * 100

# Boundary conditions for the last section (all time steps)
rho_sezfin = np.full(K, 80)
fi_sezfin = rho_sezfin * 100

# ON-ramp Demand Profile (same for both cases)
dem = np.zeros((N, K))
# High demand for first 100 time steps in on-ramp sections:
for k in range(100):
    dem[4, k] = 2000   # Section 5
    dem[9, k] = 2000   # Section 10
    dem[11, k] = 2000  # Section 12
# Lower demand afterwards:
for k in range(100, K):
    dem[4, k] = 100    # Section 5
    dem[9, k] = 100    # Section 10
    dem[11, k] = 50    # Section 12

# Simulation: No Control
# In this case, the ramp flow is simply equal to the demand.
rho_nc = np.zeros((N, K + 1))  # densities (no control)
fi_nc = np.zeros((N, K + 1))   # interface flows (no control)
r_nc = dem.copy()             # ramp flows: unchanged from demand

# Initialize densities at time step 0 (using given initial conditions)
for i in range(N - 1):
    rho_nc[i, 0] = rho_iniz[i] + T / L[i] * (fi_iniz[i] - fi_iniz[i + 1])
rho_nc[N - 1, 0] = rho_iniz[N - 1] + T / L[N - 1] * (fi_iniz[N - 1] - vf[N - 1] * rho_iniz[N - 1])

# Time stepping loop for No Control
for k in range(K):
    # Compute interface flow at upstream boundary
    fi_nc[0, k] = min(fi_sez0[k], w[0] * (rho_max[0] - rho_nc[0, k]) - r_nc[0, k])
    # Compute interface flows for sections 2 to N
    for i in range(1, N):
        fi_nc[i, k] = min(vf[i - 1] * rho_nc[i - 1, k], w[i] * (rho_max[i] - rho_nc[i, k]) - r_nc[i, k])
    # Update densities for sections 1 to N-1
    for i in range(N - 1):
        rho_nc[i, k + 1] = rho_nc[i, k] + T / L[i] * (fi_nc[i, k] + r_nc[i, k] - fi_nc[i + 1, k])
    # Downstream boundary update for the last section
    fi_sezfin_nc = min(vf[N - 1] * rho_nc[N - 1, k], fi_sezfin[k])
    rho_nc[N - 1, k + 1] = rho_nc[N - 1, k] + T / L[N - 1] * (fi_nc[N - 1, k] + r_nc[N - 1, k] - fi_sezfin_nc)

# Simulation: ALINEA Control
# ALINEA control parameters
K_alinea = 50         # Gain for ALINEA (vehicles/h per (veh/km))
rho_ref = 60          # Target mainline density (veh/km)
r_max_control = 2000  # Maximum allowed ramp metering rate (vehicles/h)

rho_al = np.zeros((N, K + 1))  # densities (with ALINEA)
fi_al = np.zeros((N, K + 1))   # interface flows (with ALINEA)
r_al = dem.copy()             # initial ramp flows (will be updated for on-ramps)

# Initialize densities at time step 0
for i in range(N - 1):
    rho_al[i, 0] = rho_iniz[i] + T / L[i] * (fi_iniz[i] - fi_iniz[i + 1])
rho_al[N - 1, 0] = rho_iniz[N - 1] + T / L[N - 1] * (fi_iniz[N - 1] - vf[N - 1] * rho_iniz[N - 1])

# Time stepping loop with ALINEA Control
for k in range(K):
    # --- ALINEA Control Update ---
    for i in range(N):
        if Ir[i] == 1 and k > 0:
            r_al[i, k] = r_al[i, k - 1] + K_alinea * (rho_ref - rho_al[i, k])
            r_al[i, k] = np.clip(r_al[i, k], 0, r_max_control)
    # --- Compute Interface Flows ---
    fi_al[0, k] = min(fi_sez0[k], w[0] * (rho_max[0] - rho_al[0, k]) - r_al[0, k])
    for i in range(1, N):
        fi_al[i, k] = min(vf[i - 1] * rho_al[i - 1, k], w[i] * (rho_max[i] - rho_al[i, k]) - r_al[i, k])
    # --- Update Densities ---
    for i in range(N - 1):
        rho_al[i, k + 1] = rho_al[i, k] + T / L[i] * (fi_al[i, k] + r_al[i, k] - fi_al[i + 1, k])
    fi_sezfin_al = min(vf[N - 1] * rho_al[N - 1, k], fi_sezfin[k])
    rho_al[N - 1, k + 1] = rho_al[N - 1, k] + T / L[N - 1] * (fi_al[N - 1, k] + r_al[N - 1, k] - fi_sezfin_al)

# Plotting the Results Side by Side with a Common Color Scale
# Determine the common min and max for both simulations
global_min = min(np.min(rho_nc), np.min(rho_al))
global_max = max(np.max(rho_nc), np.max(rho_al))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# No Control plot with common color limits
im0 = axes[0].imshow(rho_nc, aspect='auto', origin='lower', vmin=global_min, vmax=global_max)
axes[0].set_title('No Control')
axes[0].set_xlabel('Time Step')
axes[0].set_ylabel('Section')
fig.colorbar(im0, ax=axes[0], label='Density (veh/km)')

# ALINEA Control plot with common color limits
im1 = axes[1].imshow(rho_al, aspect='auto', origin='lower', vmin=global_min, vmax=global_max)
axes[1].set_title('ALINEA Control')
axes[1].set_xlabel('Time Step')
axes[1].set_ylabel('Section')
fig.colorbar(im1, ax=axes[1], label='Density (veh/km)')

plt.tight_layout()
plt.show()

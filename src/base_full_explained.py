"""
ACTM freeway simulation – “no-control” scenario
Implements equations (13) & (14) of the slides, using the same
parameter values that were hard-coded in the original MATLAB file.
"""

# ---------------------------------------------------------------------
# 1. Imports
# ---------------------------------------------------------------------
import numpy as np              # Vector/matrix arithmetic
import matplotlib.pyplot as plt # 3-D surface plot for the density map

# ---------------------------------------------------------------------
# 2. Model-wide constants (slide titles: *Modelling freeway traffic systems*)
# ---------------------------------------------------------------------
N = 20          # Number of cells/sections along the mainline
K = 360 * 2     # Simulation length = 720 steps  (≈ 2 h because every step is 10 s)
T = 10 / 3600   # Sampling period T = 10 s expressed in hours  (10⁄3600 h)

# -------------------- geometric layout --------------------------------
Ir = np.zeros(N, dtype=int) # Ir[i] = 1  → cell i has an on-ramp
Ir[[4, 9, 11]] = 1          # on-ramps at sections 5, 10, 12  (MATLAB is 1-based)

L         = np.full(N, 0.5) # Cell length Lᵢ = 0.5 km (500 m)                (slide 1)
rho_max   = np.full(N, 300) # Jam density ρ̄ᵢ = 300 veh/km                    (slide 6)
vf        = np.full(N, 100) # Free-flow speed vᵢ = 100 km/h                  (slide 6)
w         = np.full(N, 27)  # Congestion wave speed wᵢ = 27 km/h             (slide 6)
# (No off-ramps are considered here → βᵢ = 0, so the extra terms in eqs. 10–12 vanish.)

# ---------------------------------------------------------------------
# 3. Initial / boundary conditions  (slide 1 & MATLAB code)
# ---------------------------------------------------------------------
rho  = np.zeros((N, K + 1))   # ρᵢ(k)    – mainline density trajectories
ql   = np.zeros((N, K + 1))   # lᵢ(k)    – on-ramp queue trajectories
fi   = np.zeros((N, K + 1))   # φᵢ(k)    – mainline interface flows
ramp = np.zeros((N, K))       # rᵢ(k)    – actual flow entering from the ramp

# uniform initial state (ρ = 60 veh/km ; φ = ρ · v = 60·100 = 6000 veh/h)
rho[:, 0] = 60
fi[:, 0]  = rho[:, 0] * 100   # only used for the very first density update

# Up-stream boundary (cell 0 of the slides) held constant at 60 veh/km
rho0 = np.full(K, 60)
fi0  = rho0 * 100             # sending flow from the virtual cell 0

# Down-stream boundary (after cell N) held at 80 veh/km
rhof = np.full(K, 80)
fif  = rhof * 100

# ---------------------------------------------------------------------
# 4. On-ramp *demand* profile dᵢ(k)  (slide 1 variable list)
#    – here we hard-code exactly what the MATLAB file did.
# ---------------------------------------------------------------------
dem = np.zeros((N, K))              # dᵢ(k)   – desired ramp inflow
dem[[4, 9, 11], :100]  = 2000       # heavy burst (first 100 steps ≈ 17 min)
dem[4,  100:] = 100                 # afterwards moderate flow
dem[9,  100:] = 100
dem[11, 100:] =  50

# ---------------------------------------------------------------------
# 5. First density update (k = 0 → 1) uses the *initial* interface flows
#    ρᵢ(1) = ρᵢ(0) + (T/Lᵢ) ( φᵢ(0) − φᵢ₊₁(0) )
# ---------------------------------------------------------------------
for i in range(N - 1):
    rho[i, 1] = rho[i, 0] + (T / L[i]) * (fi[i, 0] - fi[i + 1, 0])

# last cell uses downstream boundary instead of φ_{N+1}
rho[N - 1, 1] = rho[N - 1, 0] + (T / L[N - 1]) * \
                (fi[N - 1, 0] - vf[N - 1] * rho[N - 1, 0])

# ---------------------------------------------------------------------
# 6. MAIN TIME-STEPPING LOOP – implements eqs. 13 & 14 exactly
# ---------------------------------------------------------------------
for k in range(K):

    # 6-A. **Compute mainline interface flow φ₁(k) at the upstream boundary**
    #      φ₁ = min{ sending capacity of the virtual cell 0 ,
    #                receiving capacity of cell 1 minus current ramp inflow }
    fi[0, k] = min(fi0[k], w[0] * (rho_max[0] - rho[0, k]) - ramp[0, k])

    # 6-B. **Compute interface flows φᵢ(k)   (eq. 14)**
    for i in range(1, N):
        sending    = vf[i - 1] * rho[i - 1, k]               # vᵢ₋₁ ρᵢ₋₁(k)
        receiving  = w[i] * (rho_max[i] - rho[i, k])         # wᵢ(ρ̄ᵢ − ρᵢ(k))
        fi[i, k]   = min(sending, receiving - ramp[i, k])    # minus ramp term ⇒ ACTM

    # 6-C. **Update densities (eq. 13 with βᵢ = 0)** -------------------------
    for i in range(N - 1):
        rho[i, k + 1] = rho[i, k] + (T / L[i]) * \
                        (fi[i, k] + ramp[i, k] - fi[i + 1, k])

    # last cell uses the downstream boundary again
    fif[k] = min(vf[N - 1] * rho[N - 1, k], fif[k])          # (not strictly needed)
    rho[N - 1, k + 1] = rho[N - 1, k] + (T / L[N - 1]) * \
                        (fi[N - 1, k] + ramp[N - 1, k] - fif[k])

    # 6-D. **Update on-ramp queues (eq. 7 simplified)**
    ql[:, k + 1] = ql[:, k] + T * (dem[:, k] - ramp[:, k])

    # 6-E. **With no control, actual ramp flow equals demand** -------------
    ramp[:, k] = dem[:, k]   # (ALINEA or any metering logic would modify this line)

# ---------------------------------------------------------------------
# 7. Visualise results (simple 3-D mesh like the MATLAB `mesh` command)
# ---------------------------------------------------------------------
fig = plt.figure(figsize=(9, 6))
ax  = fig.add_subplot(projection='3d')
X, Y = np.meshgrid(np.arange(K + 1), np.arange(N))  # time, space grids
ax.plot_surface(X, Y, rho, rstride=1, cstride=1, linewidth=0, antialiased=False)
ax.set_xlabel('Time step  k')
ax.set_ylabel('Cell index  i')
ax.set_zlabel('Density  ρᵢ(k)  [veh/km]')
plt.show()
"""
===============================================================================
ACTM + ALINEA DEMO
-------------------------------------------------------------------------------
• Implements the Asymmetric Cell-Transmission Model (CTM slides 13-14).
• Adds an ALINEA ramp-metering controller (Feedback-controller slide (3)).
• Generates – for ONE scenario – the full set of plots you requested:
      ├─ For *each* on-ramp:   demand / density / downstream flow
      ├─ For both regimes:     (blue = no-control, red = ALINEA)
      ├─ Global average-density comparison
      └─ Density heat-maps


# --------------------------------------------------------------------------- #
# 0. STANDARD IMPORTS                                                         #
# --------------------------------------------------------------------------- #
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --------------------------------------------------------------------------- #
# 1.  SCENARIO-WIDE CONSTANTS  (slides “Modelling freeway traffic systems”)   #
# --------------------------------------------------------------------------- #
N = 20                  # number of cells                                   
K = 360 * 2             # horizon = 720 steps  (≈2 h at 10 s/step)           
T = 10 / 3600           # sampling time  T =10 s expressed in hours          
L = np.full(N, 0.5)     # cell length Lᵢ =0.5 km (500 m)                     
rho_max = np.full(N,300)# jam density ρ̄ᵢ                                     
vf       = np.full(N,100)# free-flow speed vᵢ                                
w        = np.full(N, 27)# congestion wave speed wᵢ                          
# (No off-ramps in this example ⇒ split ratios βᵢ =0 everywhere.)

# Indices of the cells that have on-ramps  (MATLAB → Python: 5→4,10→9,12→11)
ramp_idx = np.array([4, 9, 11], dtype=int)
Ir       = np.zeros(N, dtype=int)
Ir[ramp_idx] = 1        # Indicator (1 if cell has an on-ramp)               

# --------------------------------------------------------------------------- #
# 2.  INITIAL AND BOUNDARY CONDITIONS                                         #
# --------------------------------------------------------------------------- #
# Uniform “light congestion” initial state (60 veh/km)
rho0_init = 60
rho       = np.full((N, K + 1), rho0_init)   # mainline densities ρᵢ(k)
fi        = np.zeros((N, K + 1))             # interface flows φᵢ(k)
fi[:, 0]  = rho0_init * vf                   # φᵢ(0) = ρ · v  (slide 2, qᵢ=ρᵢvᵢ)

# Upstream boundary (virtual cell 0) is kept at 60 veh/km throughout.
rho_inlet = np.full(K, 60)
fi_inlet  = rho_inlet * vf[0]

# Downstream boundary (after cell N) is fixed at 80 veh/km.
rho_outlet = np.full(K, 80)
fi_outlet  = rho_outlet * vf[-1]

# --------------------------------------------------------------------------- #
# 3.  ON-RAMP *DEMAND* dᵢ(k)  – same burst profile as the MATLAB code         #
# --------------------------------------------------------------------------- #
dem = np.zeros((N, K))              # demand matrix dᵢ(k)
dem[ramp_idx[:, None], :100] = 2000 # heavy burst (first 100 steps ≈17 min)
dem[4, 100:]  = 100                 # moderate flow afterwards
dem[9, 100:]  = 100
dem[11,100:]  =  50

# --------------------------------------------------------------------------- #
# 4.  ALINEA PARAMETERS  (Feedback-controller slide (3))                      #
# --------------------------------------------------------------------------- #
rho_ref        = 40             # set-point for *down-stream* density ρᵢ⁺₁*(k)
K_alinea       = 70             # integral gain  K_R         (tune if needed)
r_min_alinea   = 0              # minimum metered flow (veh/h)
r_max_alinea   = 2200           # maximum metered flow (veh/h)  ~2 veh/3 s
# Array that will store the *controlled* ramp flows
r_alinea = np.zeros((N, K))     # rᵢᶜ(k) in the slide notation

# --------------------------------------------------------------------------- #
# 5.  HELPER FUNCTION – COMPUTE φᵢ(k) (slide 14)                              #
# --------------------------------------------------------------------------- #
def interface_flow(i, k, rho_mat, ramp_flow):
    """
    φᵢ(k) = min{ v_{i-1} ρ_{i-1}(k) , wᵢ (ρ̄ᵢ − ρᵢ(k)) − rᵢ(k) }
    Implements ACTM eq. (14).   • i ∈ [0, N-1]
    """
    if i == 0:
        # Up-stream boundary uses the inlet sending flow instead of v_{-1}ρ_{-1}
        sending = fi_inlet[k]
    else:
        sending = vf[i - 1] * rho_mat[i - 1, k]

    receiving = w[i] * (rho_max[i] - rho_mat[i, k]) - ramp_flow[i, k]
    return min(sending, receiving)

# --------------------------------------------------------------------------- #
# 6.  MAIN LOOP – SIMULATE **two** regimes in parallel                        #
# --------------------------------------------------------------------------- #
# 6-A. Containers for “no-control” scenario  (ramp = demand)
r_nc   = dem.copy()                      # ramp flow == demand (no metering)
rho_nc = rho.copy()                      # deep-copy of initial densities
fi_nc  = fi.copy()

# 6-B. Containers for **ALINEA** scenario (ramp flow set by control law)
rho_al = rho.copy()
fi_al  = fi.copy()
# r_alinea (created above) already zero-initialised.

for k in range(K):

    # ------------------------------------------------------------------ #
    # 6-C-1. NO-CONTROL – compute interface flows φᵢ(k) and advance ρ.   #
    # ------------------------------------------------------------------ #
    for i in range(N):
        fi_nc[i, k] = interface_flow(i, k, rho_nc, r_nc)  # ACTM eq. (14)

    # density update (ACTM eq. (13) with βᵢ=0)
    for i in range(N - 1):
        rho_nc[i, k + 1] = rho_nc[i, k] + (T / L[i]) * \
                           (fi_nc[i, k] + r_nc[i, k] - fi_nc[i + 1, k])

    # last cell uses downstream boundary
    fi_outlet[k]       = min(vf[-1] * rho_nc[-1, k], fi_outlet[k])
    rho_nc[-1, k + 1]  = rho_nc[-1, k] + (T / L[-1]) * \
                         (fi_nc[-1, k] + r_nc[-1, k] - fi_outlet[k])

    # ------------------------------------------------------------------ #
    # 6-C-2. ALINEA – *first* compute φᵢ(k) using current r_alinea,      #
    #         *then* update densities, *then* apply the controller to    #
    #         produce rᵢ(k+1) for the *next* step.                       #
    # ------------------------------------------------------------------ #
    for i in range(N):
        fi_al[i, k] = interface_flow(i, k, rho_al, r_alinea)

    for i in range(N - 1):
        rho_al[i, k + 1] = rho_al[i, k] + (T / L[i]) * \
                           (fi_al[i, k] + r_alinea[i, k] - fi_al[i + 1, k])

    fi_outlet[k]       = min(vf[-1] * rho_al[-1, k], fi_outlet[k])
    rho_al[-1, k + 1]  = rho_al[-1, k] + (T / L[-1]) * \
                         (fi_al[-1, k] + r_alinea[-1, k] - fi_outlet[k])

    # --------------- ALINEA CONTROL LAW  (slide (3), eq. (3)) ----------- #
    #  rᵢᶜ(k+1) = rᵢᶜ(k) + K_R * (ρ*  − ρ_down(k))
    for i in ramp_idx:
        # density immediately *down-stream* of ramp = cell i+1
        ds_idx       = min(i + 1, N - 1)
        error        = rho_ref - rho_al[ds_idx, k]  # (ρ* − ρ_down)
        r_next       = r_alinea[i, k] + K_alinea * error
        # enforce physical bounds
        r_next       = np.clip(r_next, r_min_alinea, r_max_alinea)
        if k < K - 1:
            r_alinea[i, k + 1] = r_next

    # Non-ramp cells keep r=0 for the next step automatically.

# --------------------------------------------------------------------------- #
# 7.  VISUALISATION – exactly the layout in your multi-row code snippet       #
# --------------------------------------------------------------------------- #
num_ramps  = len(ramp_idx)
total_rows = 3 * num_ramps + 1 + 1        # per-ramp rows + overall + heat-maps

fig = plt.figure(figsize=(12, 3 * total_rows))
gs  = gridspec.GridSpec(total_rows, 2, figure=fig)

time_dem   = np.arange(K)       # length-K  vectors   (flows / demand)
time_state = np.arange(K + 1)   # length K+1 (densities)

for j, ramp in enumerate(ramp_idx):
    ds = min(ramp + 1, N - 1)   # index of downstream density section

    # --- COLUMN 0 :  NO-CONTROL ----------------------------------------
    ax0 = fig.add_subplot(gs[3 * j, 0])
    ax0.plot(time_dem, r_nc[ramp, :], 'b-', lw=2)
    ax0.set(title=f'No-Ctrl – Ramp {ramp} Flow', ylabel='veh/h')

    ax1 = fig.add_subplot(gs[3 * j + 1, 0])
    ax1.plot(time_state, rho_nc[ds, :], 'b-', lw=2)
    ax1.set(title=f'No-Ctrl – ρ downstream sec {ds}', ylabel='veh/km')

    ax2 = fig.add_subplot(gs[3 * j + 2, 0])
    flow_ds_nc = fi_nc[ds, :K] if ramp < N - 1 else \
                 np.minimum(vf[-1] * rho_nc[-1, :K], fi_outlet)
    ax2.plot(time_dem, flow_ds_nc, 'b-', lw=2)
    ax2.set(title=f'No-Ctrl – φ downstream sec {ds}', ylabel='veh/h',
            xlabel='time step k')

    # --- COLUMN 1 :  ALINEA -------------------------------------------
    bx0 = fig.add_subplot(gs[3 * j, 1])
    bx0.plot(time_dem, r_alinea[ramp, :], 'r-', lw=2)
    bx0.set(title=f'ALINEA – Ramp {ramp} Flow')

    bx1 = fig.add_subplot(gs[3 * j + 1, 1])
    bx1.plot(time_state, rho_al[ds, :], 'r-', lw=2)
    bx1.set(title=f'ALINEA – ρ downstream sec {ds}')

    bx2 = fig.add_subplot(gs[3 * j + 2, 1])
    flow_ds_al = fi_al[ds, :K] if ramp < N - 1 else \
                 np.minimum(vf[-1] * rho_al[-1, :K], fi_outlet)
    bx2.plot(time_dem, flow_ds_al, 'r-', lw=2)
    bx2.set(title=f'ALINEA – φ downstream sec {ds}', xlabel='time step k')

# -------- Overall average-density comparison -------------------------------
row_ov = 3 * num_ramps
ax_ov  = fig.add_subplot(gs[row_ov, :])
ax_ov.plot(time_state, rho_nc.mean(0), 'b-', lw=2, label='No-control')
ax_ov.plot(time_state, rho_al.mean(0), 'r-', lw=2, label='ALINEA')
ax_ov.set(title='Network-average density', ylabel='veh/km', xlabel='time step k')
ax_ov.legend()

# -------- Density heat-maps -------------------------------------------------
row_hm = row_ov + 1
vmin   = min(rho_nc.min(), rho_al.min())
vmax   = max(rho_nc.max(), rho_al.max())

ax_hm0 = fig.add_subplot(gs[row_hm, 0])
im0 = ax_hm0.imshow(rho_nc, aspect='auto', origin='lower',
                    vmin=vmin, vmax=vmax)
ax_hm0.set(title='No-Ctrl – Density Heat-map', ylabel='cell i')
fig.colorbar(im0, ax=ax_hm0, label='veh/km')

ax_hm1 = fig.add_subplot(gs[row_hm, 1])
im1 = ax_hm1.imshow(rho_al, aspect='auto', origin='lower',
                    vmin=vmin, vmax=vmax)
ax_hm1.set(title='ALINEA – Density Heat-map', ylabel='cell i', xlabel='time step k')
fig.colorbar(im1, ax=ax_hm1, label='veh/km')

fig.suptitle('ACTM + ALINEA – Single-scenario visualisation', fontsize=16, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

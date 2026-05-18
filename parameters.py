import math

# ---- Optical parameters ----
lam_m  = 1550e-9       # wavelength (m)
lam_um = lam_m * 1e6   # wavelength (um) - used in Mie coefficient equations
lam_nm = lam_m * 1e9   # wavelength (nm) - used in geometrical scattering

# ---- TX / RX geometry ----
theta_T = 15e-6   # TX beam divergence angle (rad)
D_R     = 0.08    # receiver aperture diameter (m)
err_T   = 1e-6    # TX pointing error (rad)
err_R   = 1e-6    # RX pointing error (rad)

# ---- Orbit / geometry ----
h_E  = 1.0      # ground station altitude (km)
h_S  = 550.0    # satellite altitude (km)
h_A  = 20.0     # troposphere height (km)
R_E  = 6378.1   # Earth radius (km)
d_SS = 1000.0   # ISL distance (km)

theta_E     = 40.0                  # elevation angle (degrees)
theta_E_rad = math.radians(theta_E) # elevation angle (radians)

# ---- Atmospheric parameters ----
L_W = 0.5       # liquid water content (g/m3)
N   = 3.128e-4  # cloud droplet concentration (cm-3)
phi = 1.6       # particle size coefficient

# ---- Link budget ----
P_T_inter  = 15.32  # transmitted power ISL (dBm)
P_T_updown = 13.98  # transmitted power uplink/downlink (dBm)
eta_T = -0.97       # transmitter efficiency (dB)
eta_R = -0.97       # receiver efficiency (dB)
P_req = -35.5       # receiver sensitivity at 10 Gbps (dBm)

# ---- APD / SNR parameters - from Table A1, Giggenbach 2022 ----
q   = 1.602e-19  # elementary charge (C)
h   = 6.626e-34  # Planck constant (Js)
c   = 2.998e8    # speed of light (m/s)
eta = 0.8        # quantum efficiency
M   = 10         # APD multiplication factor
kA  = 0.2        # ionisation ratio
Idm = 2.5e-9     # dark current - multiplied portion (A)
it  = 6.6e-12    # TIA thermal noise density (A/sqrt(Hz))
r   = 10e9       # data rate (10 Gbps)
B   = r / 2      # receiver bandwidth (Hz)

PBGL_inter  = 0.0    # background light - zero for intersatellite (no atmosphere)
PBGL_updown = 50e-9  # background light for up/downlink - 50nW typical daytime (Table 2)

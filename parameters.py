import math

# References
# [Liang]      Liang et al. (2022). Link Budget Analysis for Free-Space Optical Satellite Networks.
#              IEEE WoWMoM, pp. 471-476. https://doi.org/10.1109/WoWMoM54355.2022.00073
# [Giggenbach] Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors
#              for Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773.
#              https://doi.org/10.3390/s22186773
# [Arnon]      Arnon, S. (2005). Performance of a uSatellite Network with an Optical Preamplifier.
#              JOSA, 22(4), pp. 708-715. https://doi.org/10.1364/JOSAA.22.000708
# [LLCD]       Boroson et al. (2014). Overview and Results of the Lunar Laser Communications
#              Demonstration. SPIE Vol. 8971, 89710S. https://doi.org/10.1117/12.2045508
# [LLOT]       Sodnik et al. (2014). LLCD operations using the Lunar Lasercom OGS Terminal.
#              SPIE Vol. 8971, 89710W. https://doi.org/10.1117/12.2045510
# [Psyche]     Casanovas Ventura, M. (2021). Study: An assessment on the requirements for
#              deep space optical communications. Master's thesis, Universitat Politècnica
#              de Catalunya (ESEIAAT). Directed by J. Gutiérrez Cabello & M. Soria Guerrero.


class LiangParameters:
    # Liang et al. (2022) — LEO satellite FSO link

    lam_m = 1550e-9     # wavelength (m)                          [Table 1]

    theta_T = 15e-6     # TX beam divergence (rad)                [Table 1]
    D_R     = 0.08      # RX aperture diameter (m)                [Table 1]
    err_T   = 1e-6      # TX pointing error (rad)                 [Table 1]
    err_R   = 1e-6      # RX pointing error (rad)                 [Table 1]

    h_E  = 1.0          # ground station altitude (km)            [Table 1]
    h_S  = 550.0        # satellite altitude (km)                 [Table 1]
    h_A  = 20.0         # troposphere height (km)                 [Table 1]
    R_E  = 6378.1       # Earth radius (km)
    d_SS = 1000.0       # ISL distance (km)                       [Table 2]
    theta_E = 40.0      # elevation angle (deg)                   [Table 5]

    L_W = 0.5           # cirrus cloud concentration (cm-3)       [Table 3]
    N   = 3.128e-4      # liquid water content (g/m-3)            [Table 3]
    phi = 1.6           # particle size distribution coefficient   [Table 3]

    P_T_inter  = 15.32  # ISL TX power (dBm)                      [Table 2]
    P_T_updown = 13.98  # up/downlink TX power (dBm)              [Table 4]
    eta_T = -0.97       # TX optical efficiency (dB)              [Table 1]
    eta_R = -0.97       # RX optical efficiency (dB)              [Table 1]
    P_req = -35.5       # minimum received power (dBm)            [Table 1]

    # APD detector [Giggenbach, Table A1]
    eta = 0.8           # quantum efficiency
    M   = 10            # APD multiplication factor
    kA  = 0.2           # ionisation ratio
    Idm = 2.5e-9        # dark current (A)
    it  = 6.6e-12       # TIA thermal noise current density (A/sqrt(Hz))

    r = 10e9            # data rate (bps)                         [Table 1]

    PBGL_inter  = 0.0   # ISL background light (W) — no atmosphere between satellites
    PBGL_updown = 50e-9 # up/downlink background light (W)        [Giggenbach, Table 2]


class LLCDParameters:
    # Boroson et al. (2014) — Earth-Moon FSO downlink (space terminal → OGS)
    #
    # Compatibility notes:
    #   - Detector: SNSPD (photon counting), not APD — BER/Q output is not physically valid
    #   - Modulation: 16-PPM downlink, not OOK — BER formula does not apply
    #   - Geometry: ~384,000 km Earth-Moon link; use calculate_isl (d_SS), not calculate_updownlink

    lam_m = 1550.12e-9  # downlink wavelength (m)                 [LLCD]

    theta_T = 15e-6     # TX beam divergence (rad) — 10 cm space terminal
    D_R     = 1.016     # RX aperture diameter (m) — OGS telescope [LLOT]
    err_T   = 2e-6      # TX pointing error (rad) — MIRU-stabilised [LLCD]
    err_R   = 0.0       # RX pointing error — Arnon formula invalid at this scale; beam is ~2.9 km wide at Moon distance

    h_E  = 2.393        # OGS altitude (km)                       [LLOT]
    h_S  = 550.0        # dummy — lunar link uses d_SS; prevents divide-by-zero in channel()
    h_A  = 20.0         # troposphere height (km)
    R_E  = 6378.1       # Earth radius (km)
    d_SS = 384000.0     # Earth-Moon distance (km) — midpoint of 362,570–405,410 km range [LLCD]
    theta_E = 45.0      # elevation angle (deg) — assumed

    L_W = 0.5           # cirrus cloud concentration (cm-3)
    N   = 3.128e-4      # liquid water content (g/m-3)
    phi = 1.6           # particle size distribution coefficient

    P_T_inter  = 27.0   # space terminal TX power (dBm) — 0.5 W downlink [LLCD]
    P_T_updown = 50.8   # ground TX power (dBm) — 3×40 W uplink          [LLCD]
    eta_T = -1.7        # TX optical efficiency (dB) — TX loss + Strehl   [LLCD]
    eta_R = -3.15       # RX optical efficiency (dB) — ATC RX loss        [LLOT]
    P_req = -73.7       # minimum received power at ATC (dBm)             [LLOT]

    # SNSPD detector — placeholder values; APD noise model does not apply
    eta = 0.50          # detection efficiency                    [LLCD]
    M   = 1
    kA  = 0.0
    Idm = 0.0
    it  = 1e-12         # near-zero; kept non-zero for numerical stability

    r = 155e6           # downlink data rate (bps) — 155 Mbps    [LLOT]

    PBGL_inter  = 0.0
    PBGL_updown = 1.77e-13  # sky + lunar disc background at ATC (W) [LLOT]


class PsycheParameters:
    # NASA Psyche mission — deep-space optical downlink (spacecraft → Palomar-Hale)
    # Parameters from master's thesis link budget table [Psyche] chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://upcommons.upc.edu/server/api/core/bitstreams/2910710a-d79d-4996-8524-619dbc720c73/content
    #
    lam_m = 1550e-9     # wavelength (m)                                           [Psyche]

    theta_T = 9e-6      # TX beam divergence (rad) — derived from D_Tx = 22 cm    [Psyche]
    D_R     = 5.0       # RX aperture diameter (m) — Palomar-Hale telescope        [Psyche]
    # err_T back-calculated from thesis η_TP = 0.921 (−0.36 dB) via Arnon formula
    err_T   = 6.5e-7    # TX pointing error (rad)                                  [Psyche]
    err_R   = 0.0       # RX pointing error — thesis η_RP = 1.0; Arnon formula also invalid at planetary scale

    h_E  = 1.713        # Palomar Observatory altitude (km)
    h_S  = 550.0        # dummy — prevents divide-by-zero in channel(); use d_SS for this link
    h_A  = 20.0         # troposphere height (km)
    R_E  = 6378.1       # Earth radius (km)
    d_SS = 4.488e8      # 3 AU in km — worst-case range                            [Psyche]
    theta_E = 45.0      # elevation angle (deg) — assumed

    L_W = 0.5           # cirrus cloud concentration (cm-3) — indicative only
    N   = 3.128e-4      # liquid water content (g/m-3)      — indicative only
    phi = 1.6           # particle size distribution        — indicative only

    P_T_inter  = 36.02  # spacecraft TX power (dBm) — 4 W                         [Psyche]
    P_T_updown = 36.02  # same
    # eta_T = g_T (−0.893) + Strehl loss (−1.72) + η_Tx (−1.54) = −4.15 dB       [Psyche]
    eta_T = -4.15
    # eta_R = g_R (−0.17) + L_pol (−0.3) + NBF (−1.54) + η_Rx (−1.54) = −3.55 dB [Psyche]
    eta_R = -3.55
    P_req = -100.0      # placeholder — detector sensitivity threshold not stated in thesis

    # Detector: type not specified — APD placeholders from Liang/Giggenbach
    eta = 0.8           # quantum efficiency
    M   = 10            # APD multiplication factor
    kA  = 0.2           # ionisation ratio
    Idm = 2.5e-9        # dark current (A)
    it  = 6.6e-12       # TIA thermal noise current density (A/sqrt(Hz))

    r = 1e6             # data rate (bps) — placeholder; not stated in thesis

    PBGL_inter  = 0.0   # no background between spacecraft and ground
    PBGL_updown = 50e-9 # ground background (W) — placeholder


# Change this line to switch parameter sets
active_params = LiangParameters()


# Physical constants
q = 1.602e-19   # elementary charge (C)
h = 6.626e-34   # Planck constant (Js)
c = 2.998e8     # speed of light (m/s)

# Derived from active_params
lam_m  = active_params.lam_m
lam_um = active_params.lam_m * 1e6    # wavelength (um) — used in Mie scattering
lam_nm = active_params.lam_m * 1e9    # wavelength (nm) — used in geometrical scattering

theta_T = active_params.theta_T
D_R     = active_params.D_R
err_T   = active_params.err_T
err_R   = active_params.err_R

h_E  = active_params.h_E
h_S  = active_params.h_S
h_A  = active_params.h_A
R_E  = active_params.R_E
d_SS = active_params.d_SS

theta_E     = active_params.theta_E
theta_E_rad = math.radians(active_params.theta_E)

L_W = active_params.L_W
N   = active_params.N
phi = active_params.phi

P_T_inter  = active_params.P_T_inter
P_T_updown = active_params.P_T_updown
eta_T      = active_params.eta_T
eta_R      = active_params.eta_R
P_req      = active_params.P_req

eta = active_params.eta
M   = active_params.M
kA  = active_params.kA
Idm = active_params.Idm
it  = active_params.it

r = active_params.r
B = active_params.r / 2    # receiver bandwidth (Hz) — Nyquist: B = r/2 for OOK

PBGL_inter  = active_params.PBGL_inter
PBGL_updown = active_params.PBGL_updown

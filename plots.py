import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from parameters import *
from functions import tx_parameters, rx_parameters, channel, compute_snr, calculate_updownlink

# ---- References ----
# [Liang]      Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022).
#              Link Budget Analysis for Free-Space Optical Satellite Networks.
#              IEEE WoWMoM, pp. 471-476.
#              https://doi.org/10.1109/WoWMoM54355.2022.00073
# [Giggenbach] Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors
#              for Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773.
#              https://doi.org/10.3390/s22186773


# ==============================================================================
# PLOT 1: RECEIVER SENSITIVITY CURVE
# Sweeps received power from a low value up to the actual received power,
# and plots BER vs received power for both ISL and uplink/downlink cases.
# This shows the "sensitivity" of the receiver — i.e. how BER improves as
# signal power increases, and where the 1e-9 BER target is met.
# Returns the arrays so they can be printed as a table in main.py.
# ==============================================================================

def plot_sensitivity(P_R_inter, P_R_updown):
    # Sweep received power from -65 dBm to -20 dBm
    P_rx_dBm = np.linspace(-65, -20, 500)
    P_rx_W   = 10 ** (P_rx_dBm / 10) * 1e-3  # convert dBm to Watts

    # Compute BER at each power level for both link types
    BER_inter  = []
    BER_updown = []
    for P in P_rx_W:
        _, _, ber_i = compute_snr(P, PBGL_inter)
        _, _, ber_u = compute_snr(P, PBGL_updown)
        BER_inter.append(max(ber_i, 1e-300))
        BER_updown.append(max(ber_u, 1e-300))

    fig, ax = plt.subplots()

    ax.semilogy(P_rx_dBm, BER_inter,  label='Intersatellite (no background light)')
    ax.semilogy(P_rx_dBm, BER_updown, label='Up/Downlink (with background light)')

    # Primary BER target - Liang et al. (2022), Table 1 (OOK at 10 Gbps, from which P_req = -35.5 dBm is derived)
    ax.axhline(y=1e-12, color='red',    linestyle='--', linewidth=1, label='BER = 1e-12 target (Liang et al., 2022)')
    # Secondary BER reference - Giggenbach (2022), Section 1
    ax.axhline(y=1e-9,  color='orange', linestyle='--', linewidth=1, label='BER = 1e-9 reference (Giggenbach, 2022)')

    # Expected received power for each link type
    ax.axvline(x=P_R_inter,  color='blue',  linestyle=':', linewidth=1.5, label=f'Expected received power - ISL ({round(P_R_inter, 1)} dBm)')
    ax.axvline(x=P_R_updown, color='green', linestyle=':', linewidth=1.5, label=f'Expected received power - Up/Downlink ({round(P_R_updown, 1)} dBm)')

    ax.set_xlabel('Received Power (dBm)')
    ax.set_ylabel('Bit Error Rate (BER)')
    ax.set_title(f'Receiver Sensitivity Curve\n'
                 f'Wavelength = {lam_nm:.0f} nm, Data rate = {r/1e9:.0f} Gbps OOK  |  '
                 f'ISL: Distance = {d_SS:.0f} km, TX power = {P_T_inter} dBm  |  '
                 f'Up/Downlink: Satellite altitude = {h_S:.0f} km, TX power = {P_T_updown} dBm', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylim([1e-15, 1])

    plt.tight_layout()
    plt.savefig('sensitivity_curve.png', dpi=150)

    return P_rx_dBm, BER_inter, BER_updown


def print_sensitivity_table(P_rx_dBm, BER_inter, BER_updown, step_dBm=5):
    # Print a table of BER vs received power, sampled every step_dBm dB
    # step_dBm controls how often a row is printed (default: every 5 dBm)
    print("\nSensitivity Table (BER vs Received Power):")
    print(f"{'Power (dBm)':>12}  {'BER - ISL':>12}  {'BER - Up/Downlink':>18}")
    print("-" * 48)

    # Find indices where the power is close to a multiple of step_dBm
    P_rx_dBm = np.array(P_rx_dBm)
    for target in np.arange(math.ceil(P_rx_dBm[0] / step_dBm) * step_dBm,
                            P_rx_dBm[-1] + step_dBm, step_dBm):
        idx = np.argmin(np.abs(P_rx_dBm - target))
        print(f"{P_rx_dBm[idx]:>12.1f}  {BER_inter[idx]:>12.2e}  {BER_updown[idx]:>18.2e}")


# ==============================================================================
# PLOT 2: BER AND LINK MARGIN VS ELEVATION ANGLE
# Sweeps the elevation angle from 5 deg (near horizon) to 90 deg (directly overhead).
# As the angle decreases, the slant distance and atmospheric path both increase,
# degrading the link. This plot shows how quickly performance falls off.
# Returns the arrays so they can be printed as a table in main.py.
# ==============================================================================

def plot_vs_elevation():
    # Sweep elevation angle from 5 to 90 degrees
    angles_deg = np.linspace(5, 90, 300)

    # Pre-compute TX/RX gains — these don't change with elevation angle
    G_T, L_T = tx_parameters(theta_T, err_T)
    G_R, L_R = rx_parameters(D_R, lam_m, err_R)

    BER_list = []
    LM_list  = []
    Q_list   = []

    for angle in angles_deg:
        angle_rad = math.radians(angle)

        # Recompute channel losses for this elevation angle
        d_GS, L_PS, L_PG, L_A, _, _ = channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, angle_rad, L_W, N, phi)

        # Received power and link margin for up/downlink at this angle
        P_R = P_T_updown + eta_T + eta_R + G_T + G_R + L_T + L_R + L_A + L_PG
        LM  = P_R - P_req
        P_W = 10 ** (P_R / 10) * 1e-3  # convert dBm to Watts

        Q, _, BER = compute_snr(P_W, PBGL_updown)

        # Clamp BER to a plottable floor — erfc() returns exactly 0.0 for Q > ~37,
        # which semilogy() cannot display (log(0) = -inf). 1e-300 is near the double
        # precision limit and keeps the curve visible without distorting the result.
        BER_list.append(max(BER, 1e-300))
        LM_list.append(LM)
        Q_list.append(Q)

    # ---- Subplot 1: BER vs elevation angle ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax1.semilogy(angles_deg, BER_list, color='blue')
    # Primary BER target - Liang et al. (2022), Table 1
    ax1.axhline(y=1e-12, color='red',    linestyle='--', linewidth=1, label='BER = 1e-12 target (Liang et al., 2022)')
    # Secondary BER reference - Giggenbach (2022), Section 1
    ax1.axhline(y=1e-9,  color='orange', linestyle='--', linewidth=1, label='BER = 1e-9 reference (Giggenbach, 2022)')
    ax1.axvline(x=theta_E, color='grey', linestyle=':', linewidth=1, label=f'Current angle ({theta_E} deg)')
    ax1.set_ylabel('Bit Error Rate (BER)')
    ax1.set_title(f'Up/Downlink Performance vs Elevation Angle\n'
                  f'Wavelength = {lam_nm:.0f} nm, Data rate = {r/1e9:.0f} Gbps OOK  |  '
                  f'Satellite altitude = {h_S:.0f} km, TX power = {P_T_updown} dBm, Aperture diameter = {D_R*100:.0f} cm', fontsize=9)
    ax1.legend(fontsize=8)
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # ---- Subplot 2: Link margin vs elevation angle ----
    ax2.plot(angles_deg, LM_list, color='green')
    ax2.axvline(x=theta_E, color='grey', linestyle=':', linewidth=1, label=f'Current angle ({theta_E} deg)')
    ax2.set_xlabel('Elevation Angle (degrees)')
    ax2.set_ylabel('Link Margin (dB)')
    ax2.legend(fontsize=8)
    ax2.grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig('ber_vs_elevation.png', dpi=150)

    return angles_deg, BER_list, LM_list, Q_list


# Plots BER and link margin vs pointing error, sweeping TX and RX error separately
# so the effect of uplink (ground TX) and downlink (satellite TX) errors can be compared.
def plot_pointing_error_sensitivity():
    errors  = np.logspace(math.log10(1e-7), math.log10(2e-5), 300)
    errors_urad = errors * 1e6

    G_T_lin = 16 / theta_T**2
    G_T     = 10 * math.log10(G_T_lin)
    G_R_lin = (D_R * math.pi / lam_m)**2
    G_R     = 10 * math.log10(G_R_lin)

    d_GS, L_PS, L_PG, L_A, _, _ = channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, theta_E_rad, L_W, N, phi)

    # ---- Sweep TX pointing error (err_T), keep RX error fixed at baseline ----
    # Represents uplink: ground station (TX) pointing accuracy varies
    L_R_fixed = 10 * math.log10(math.exp(-G_R_lin * err_R**2))
    BER_TX, LM_TX = [], []
    for err in errors:
        L_T_lin = math.exp(-G_T_lin * err**2)
        if L_T_lin <= 0:
            BER_TX.append(0.5); LM_TX.append(-999); continue
        L_T = 10 * math.log10(L_T_lin)
        P_R, LM, P_W, Q, SNR, BER = calculate_updownlink(G_T, L_T, G_R, L_R_fixed, L_PG, L_A)
        BER_TX.append(max(BER, 1e-300)); LM_TX.append(LM)

    # ---- Sweep RX pointing error (err_R), keep TX error fixed at baseline ----
    # Represents downlink: ground station (RX) pointing accuracy varies
    L_T_fixed = 10 * math.log10(math.exp(-G_T_lin * err_T**2))
    BER_RX, LM_RX = [], []
    for err in errors:
        L_R_lin = math.exp(-G_R_lin * err**2)
        if L_R_lin <= 0:
            BER_RX.append(0.5); LM_RX.append(-999); continue
        L_R = 10 * math.log10(L_R_lin)
        P_R, LM, P_W, Q, SNR, BER = calculate_updownlink(G_T, L_T_fixed, G_R, L_R, L_PG, L_A)
        BER_RX.append(max(BER, 1e-300)); LM_RX.append(LM)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # ---- Subplot 1: BER vs pointing error ----
    ax1.semilogy(errors_urad, BER_TX, color='blue',  label='TX pointing error (uplink — ground station transmits)')
    ax1.semilogy(errors_urad, BER_RX, color='green', label='RX pointing error (downlink — ground station receives)')
    ax1.axhline(y=1e-12, color='red',    linestyle='--', linewidth=1, label='BER = 1e-12 target (Liang et al., 2022)')
    ax1.axhline(y=1e-9,  color='orange', linestyle='--', linewidth=1, label='BER = 1e-9 reference (Giggenbach, 2022)')
    ax1.axvline(x=err_T * 1e6, color='grey', linestyle=':', linewidth=1, label=f'Current pointing error ({err_T*1e6:.1f} µrad)')
    ax1.set_ylabel('Bit Error Rate (BER)')
    ax1.set_title(f'Pointing Error Sensitivity — TX (Uplink) vs RX (Downlink)\n'
                  f'Wavelength = {lam_nm:.0f} nm, Data rate = {r/1e9:.0f} Gbps OOK  |  '
                  f'Satellite altitude = {h_S:.0f} km, TX beam divergence = {theta_T*1e6:.0f} µrad', fontsize=9)
    ax1.legend(fontsize=8)
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # ---- Subplot 2: Link margin vs pointing error ----
    LM_TX_clipped = [max(lm, -50) for lm in LM_TX]
    LM_RX_clipped = [max(lm, -50) for lm in LM_RX]
    ax2.plot(errors_urad, LM_TX_clipped, color='blue',  label='TX pointing error (uplink)')
    ax2.plot(errors_urad, LM_RX_clipped, color='green', label='RX pointing error (downlink)')
    ax2.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Link closure threshold (0 dB)')
    ax2.axvline(x=err_T * 1e6, color='grey', linestyle=':', linewidth=1, label=f'Current pointing error ({err_T*1e6:.1f} µrad)')
    ax2.set_xlabel('Pointing Error (µrad)')
    ax2.set_ylabel('Link Margin (dB)')
    ax2.legend(fontsize=8)
    ax2.grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig('pointing_error_sensitivity.png', dpi=150)


# Plots BER vs SNR for OOK modulation, with operating points marked for ISL and up/downlink.
def plot_ber_vs_snr(SNR_inter, SNR_updown, BER_inter, BER_updown):
    # Sweep Q-factor from 0.1 to 15, derive SNR (Q^2) and BER for each value
    # For OOK: SNR = Q^2 and BER = 0.5 * erfc(Q / sqrt(2))
    Q_range   = np.linspace(0.1, 15, 500)
    SNR_range = Q_range**2
    SNR_dB    = 10 * np.log10(SNR_range)
    BER_range = np.maximum(0.5 * erfc(Q_range / np.sqrt(2)), 1e-300)

    fig, ax = plt.subplots()

    ax.semilogy(SNR_dB, BER_range, color='blue', label='Theoretical BER (OOK)')

    # Mark BER targets
    ax.axhline(y=1e-12, color='red',    linestyle='--', linewidth=1, label='BER = 1e-12 target (Liang et al., 2022)')
    ax.axhline(y=1e-9,  color='orange', linestyle='--', linewidth=1, label='BER = 1e-9 reference (Giggenbach, 2022)')

    # Mark operating points for each link type
    ax.scatter(10 * math.log10(SNR_inter),  BER_inter,  color='blue',  zorder=5, label=f'ISL operating point (SNR = {round(10*math.log10(SNR_inter), 1)} dB)')
    ax.scatter(10 * math.log10(SNR_updown), BER_updown, color='green', zorder=5, label=f'Up/Downlink operating point (SNR = {round(10*math.log10(SNR_updown), 1)} dB)')

    ax.set_xlabel('SNR (dB)')
    ax.set_ylabel('Bit Error Rate (BER)')
    ax.set_title(f'BER vs SNR — OOK Modulation\n'
                 f'Wavelength = {lam_nm:.0f} nm, Data rate = {r/1e9:.0f} Gbps  |  '
                 f'Aperture diameter = {D_R*100:.0f} cm, Satellite altitude = {h_S:.0f} km', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plt.savefig('ber_vs_snr.png', dpi=150)


def print_elevation_table(angles_deg, BER_list, LM_list, Q_list, step_deg=10):
    # Print a table of BER, link margin, and Q-factor vs elevation angle,
    # sampled every step_deg degrees (default: every 10 degrees)
    print("\nElevation Angle Sweep Table (Up/Downlink):")
    print(f"{'Angle (deg)':>12}  {'BER':>12}  {'Link Margin (dB)':>17}  {'Q-factor':>10}")
    print("-" * 58)

    angles_deg = np.array(angles_deg)
    for target in range(10, 91, step_deg):
        idx = np.argmin(np.abs(angles_deg - target))
        print(f"{angles_deg[idx]:>12.1f}  {BER_list[idx]:>12.2e}  {LM_list[idx]:>17.2f}  {Q_list[idx]:>10.3f}")

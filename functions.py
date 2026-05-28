import math
from scipy.special import erfc
from parameters import q, h, c, eta, lam_m, M, kA, Idm, it, B
from parameters import P_T_inter, P_T_updown, eta_T, eta_R, P_req, PBGL_inter, PBGL_updown

# ---- References ----
# [Liang]      Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022).
#              Link Budget Analysis for Free-Space Optical Satellite Networks.
#              IEEE WoWMoM, pp. 471-476.
#              https://doi.org/10.1109/WoWMoM54355.2022.00073
# [Giggenbach] Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors
#              for Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773.
#              https://doi.org/10.3390/s22186773
# [Arnon]      Arnon, S. (2005). Performance of a uSatellite Network with an Optical Preamplifier.
#              Journal of Optical Society of America, 22(4), pp. 708-715.
#              https://doi.org/10.1364/JOSAA.22.000708


# Calculates transmitter gain (G_T) and TX pointing loss (L_T) from beam divergence angle and pointing error.
def tx_parameters(theta_T, err_T):
    G_T_lin = 16 / theta_T**2                   # transmitter gain (linear) - [Liang] Eq. 3
    G_T     = 10 * math.log10(G_T_lin)          # transmitter gain (dB)

    L_T_lin = math.exp(-G_T_lin * err_T**2)     # TX pointing loss (linear) - [Liang] Section III-A, citing [Arnon]
    L_T     = 10 * math.log10(L_T_lin)          # TX pointing loss (dB)

    return G_T, L_T


# Calculates receiver gain (G_R) and RX pointing loss (L_R) from aperture diameter and pointing error.
def rx_parameters(D_R, lam_m, err_R):
    G_R_lin = (D_R * math.pi / lam_m)**2        # receiver gain (linear) - [Liang] Section III-A
    G_R     = 10 * math.log10(G_R_lin)          # receiver gain (dB)

    L_R_lin = math.exp(-G_R_lin * err_R**2)     # RX pointing loss (linear) - [Liang] Section III-A, citing [Arnon]
    L_R     = 10 * math.log10(L_R_lin)          # RX pointing loss (dB)

    return G_R, L_R


# Calculates slant distance, free-space path loss (ISL and up/downlink), and atmospheric loss from geometry and atmospheric parameters.
def channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, theta_E_rad, L_W, N, phi):
    # Slant distance - [Liang] Eq. 4
    R    = R_E + h_E
    H    = h_S - h_E
    d_GS = R * (math.sqrt(((R + H) / R)**2 - math.cos(theta_E_rad)**2) - math.sin(theta_E_rad))

    # Free-space path loss - [Liang] Eq. 2 (ISL) and Eq. 5 (up/downlink)
    L_PS = 10 * math.log10((lam_m / (4 * math.pi * d_SS * 1e3))**2)
    L_PG = 10 * math.log10((lam_m / (4 * math.pi * d_GS * 1e3))**2)

    # Mie scattering attenuation - [Liang] Eq. 6-7
    coef_a = -0.000545 * lam_um**2 + 0.002  * lam_um - 0.0038
    coef_b =  0.00628  * lam_um**2 - 0.0232 * lam_um + 0.00439
    coef_c = -0.028    * lam_um**2 + 0.101  * lam_um - 0.18
    coef_d = -0.228    * lam_um**3 + 0.922  * lam_um**2 - 1.26 * lam_um + 0.719
    rho    = coef_a * h_E**3 + coef_b * h_E**2 + coef_c * h_E + coef_d
    I_m    = math.exp(-rho / math.sin(theta_E_rad))

    # Geometrical scattering attenuation - [Liang] Eq. 8-10
    V       = 1.002 / (L_W * N)**0.6473
    theta_A = (3.91 / V) * (lam_nm / 550)**(-phi)
    d_A     = (h_A - h_E) / math.sin(theta_E_rad)
    I_g     = math.exp(-theta_A * d_A)

    # Total atmospheric loss - [Liang] Eq. 11
    L_A   = 10 * math.log10(I_m * I_g)
    I_m_dB = 10 * math.log10(I_m)   # Mie scattering loss in dB (for validation)
    I_g_dB = 10 * math.log10(I_g)   # geometrical scattering loss in dB (for validation)

    return d_GS, L_PS, L_PG, L_A, I_m_dB, I_g_dB


# Calculates Q-factor, SNR, and BER from received power and background light using the APD noise model.
def compute_snr(P_rx_W, PBGL):
    R_resp = eta * (q * lam_m) / (h * c)               # responsivity (A/W) - [Giggenbach] Eq. 1
    FA     = kA * M + 2 * (1 - kA) + (kA - 1) / M     # excess noise factor - [Giggenbach] Eq. 4

    sigma_t  = it * math.sqrt(B)                                                            # thermal noise - [Giggenbach] Eq. 14
    sigma_s0 = math.sqrt(B * 2 * q * M**2 * FA * (R_resp * PBGL + Idm))                    # shot noise for '0' bit - [Giggenbach] Eq. 15
    sigma_s1 = math.sqrt(B * 2 * q * (M**2 * FA * (R_resp * (2 * P_rx_W + PBGL) + Idm)))  # shot noise for '1' bit - [Giggenbach] Eq. 16

    Q   = (M * R_resp * 2 * P_rx_W) / (math.sqrt(sigma_s0**2 + sigma_t**2) + math.sqrt(sigma_s1**2 + sigma_t**2))  # Q-factor - [Giggenbach] Eq. 17
    SNR = Q**2                                                                              # SNR = Q^2 - [Giggenbach] Appendix A.3
    BER = 0.5 * erfc(Q / math.sqrt(2))                                                     # bit error rate - [Giggenbach] Eq. 13

    return Q, SNR, BER


# Calculates received power, link margin, and SNR for the intersatellite link using TX/RX gains and free-space path loss.
def calculate_isl(G_T, L_T, G_R, L_R, L_PS):
    P_R = P_T_inter + eta_T + eta_R + G_T + G_R + L_T + L_R + L_PS  # link budget equation - [Liang] Eq. 1
    LM  = P_R - P_req                                                 # link margin
    P_W = 10 ** (P_R / 10) * 1e-3                                    # convert dBm to Watts

    Q, SNR, BER = compute_snr(P_W, PBGL_inter)

    return P_R, LM, P_W, Q, SNR, BER


# Calculates received power, link margin, and SNR for the up/downlink using TX/RX gains, free-space path loss, and atmospheric loss.
def calculate_updownlink(G_T, L_T, G_R, L_R, L_PG, L_A):
    P_R = P_T_updown + eta_T + eta_R + G_T + G_R + L_T + L_R + L_A + L_PG  # link budget equation - [Liang] Eq. 1
    LM  = P_R - P_req                                                        # link margin
    P_W = 10 ** (P_R / 10) * 1e-3                                           # convert dBm to Watts

    Q, SNR, BER = compute_snr(P_W, PBGL_updown)

    return P_R, LM, P_W, Q, SNR, BER


# Prints a summary of system parameters and all computed gains and losses.
def print_system_params(lam_nm, h_S, d_SS, theta_E, d_GS, G_T, G_R, L_T, L_R, L_PS, L_PG, L_A):
    print("System parameters:")
    print(f"Wavelength: {lam_nm:.0f} nm, Satellite altitude: {h_S:.0f} km, ISL distance: {d_SS:.0f} km, Elevation angle: {theta_E:.0f} deg")
    print(f"Slant distance (ground to satellite): {round(d_GS, 2)} km")
    print("\nComputed gains and losses:")
    print(f"Transmitter gain G_T: {round(G_T, 2)} dB")
    print(f"Receiver gain G_R: {round(G_R, 2)} dB")
    print(f"TX pointing loss L_T: {round(L_T, 2)} dB")
    print(f"RX pointing loss L_R: {round(L_R, 2)} dB")
    print(f"Intersatellite FSPL L_PS: {round(L_PS, 2)} dB")
    print(f"Up/downlink FSPL L_PG: {round(L_PG, 2)} dB")
    print(f"Atmospheric loss L_A: {round(L_A, 2)} dB")


# Prints received power and link margin for both links, and states whether each link is closed.
def print_link_budget(P_R_inter, LM_inter, P_R_updown, LM_updown):
    print("\nLink Budget:")
    print("Intersatellite:")
    print(f"  Received Power: {round(P_R_inter, 2)} dBm")
    print(f"  Link Margin: {round(LM_inter, 2)} dB")
    print("\nUplink/Downlink:")
    print(f"  Received Power: {round(P_R_updown, 2)} dBm")
    print(f"  Link Margin: {round(LM_updown, 2)} dB")

    if LM_inter > 0:
        print(f"\nIntersatellite link is closed with a margin of {round(LM_inter, 2)} dB")
    elif LM_inter == 0:
        print("\nIntersatellite link is perfectly closed")
    else:
        print("\nIntersatellite link is not closed")

    if LM_updown > 0:
        print(f"\nUplink/Downlink is closed with a margin of {round(LM_updown, 2)} dB")
    elif LM_updown == 0:
        print("\nUplink/Downlink is perfectly closed")
    else:
        print("\nUplink/Downlink is not closed")


# Prints Q-factor, SNR, and BER for both link types.
def print_snr(B, r, Q_inter, SNR_inter, BER_inter, Q_updown, SNR_updown, BER_updown):
    print("\nSNR Analysis")
    print(f"Bandwidth: {round(B/1e9, 2)} GHz, Data rate: {round(r/1e9, 1)} Gbps")
    print("\nIntersatellite:")
    print(f"  Q-factor: {round(Q_inter, 3)}")
    print(f"  SNR: {round(SNR_inter, 2)} ({round(10 * math.log10(SNR_inter), 2)} dB)")
    print(f"  BER: {BER_inter:.2e}")
    print("\nUplink/Downlink:")
    print(f"  Q-factor: {round(Q_updown, 3)}")
    print(f"  SNR: {round(SNR_updown, 2)} ({round(10 * math.log10(SNR_updown), 2)} dB)")
    print(f"  BER: {BER_updown:.2e}")

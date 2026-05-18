import math
from scipy.special import erfc
from parameters import q, h, c, eta, lam_m, M, kA, Idm, it, B


def tx_parameters(theta_T, err_T):
    # Transmitter gain (Liang Section III-A)
    G_T_lin = 16 / theta_T**2
    G_T = 10 * math.log10(G_T_lin)

    # TX pointing loss (Liang Section III-A)
    L_T_lin = math.exp(-G_T_lin * err_T**2)
    L_T = 10 * math.log10(L_T_lin)

    return G_T, L_T


def rx_parameters(D_R, lam_m, err_R):
    # Receiver gain (Liang Section III-A)
    G_R_lin = (D_R * math.pi / lam_m)**2
    G_R = 10 * math.log10(G_R_lin)

    # RX pointing loss (Liang Section III-A)
    L_R_lin = math.exp(-G_R_lin * err_R**2)
    L_R = 10 * math.log10(L_R_lin)

    return G_R, L_R


def channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, theta_E_rad, L_W, N, phi):
    # Slant distance - ground station to satellite (Liang Eq. 4)
    R = R_E + h_E
    H = h_S - h_E
    d_GS = R * (math.sqrt(((R + H) / R)**2 - math.cos(theta_E_rad)**2) - math.sin(theta_E_rad))

    # Free space path loss (Liang Eq. 2 and Eq. 5)
    L_PS = 10 * math.log10((lam_m / (4 * math.pi * d_SS * 1e3))**2)
    L_PG = 10 * math.log10((lam_m / (4 * math.pi * d_GS * 1e3))**2)

    # Mie scattering coefficients - wavelength in um, altitude in km (Liang Eq. 6)
    coef_a = -0.000545 * lam_um**2 + 0.002  * lam_um - 0.0038
    coef_b =  0.00628  * lam_um**2 - 0.0232 * lam_um + 0.00439
    coef_c = -0.028    * lam_um**2 + 0.101  * lam_um - 0.18
    coef_d = -0.228    * lam_um**3 + 0.922  * lam_um**2 - 1.26 * lam_um + 0.719
    rho = coef_a * h_E**3 + coef_b * h_E**2 + coef_c * h_E + coef_d

    I_m = math.exp(-rho / math.sin(theta_E_rad))        # Mie scattering attenuation (Liang Eq. 7)
    V   = 1.002 / (L_W * N)**0.6473                     # visibility in km (Liang Eq. 8)
    theta_A = (3.91 / V) * (lam_nm / 550)**(-phi)       # geometrical scattering coefficient (Liang Eq. 9)
    d_A = (h_A - h_E) / math.sin(theta_E_rad)           # troposphere path distance (km)
    I_g = math.exp(-theta_A * d_A)                      # geometrical scattering attenuation (Liang Eq. 10)

    L_A = 10 * math.log10(I_m * I_g)                    # total atmospheric loss in dB (Liang Eq. 11)

    return d_GS, L_PS, L_PG, L_A


def compute_snr(P_rx_W, PBGL):
    R_resp = eta * (q * lam_m) / (h * c)               # responsivity (A/W)
    FA = kA * M + 2 * (1 - kA) + (kA - 1) / M         # excess noise factor (Eq. 4)

    sigma_t  = it * math.sqrt(B)                                                           # thermal noise (Eq. 14)
    sigma_s0 = math.sqrt(B * 2 * q * M**2 * FA * (R_resp * PBGL + Idm))                   # shot noise during "0" (Eq. 15)
    sigma_s1 = math.sqrt(B * 2 * q * (M**2 * FA * (R_resp * (2 * P_rx_W + PBGL) + Idm))) # shot noise during "1" (Eq. 16)

    Q   = (M * R_resp * 2 * P_rx_W) / (math.sqrt(sigma_s0**2 + sigma_t**2) + math.sqrt(sigma_s1**2 + sigma_t**2))  # Q-factor (Eq. 17)
    SNR = Q**2                                                                             # SNR = Q squared (Appendix A.3)
    BER = 0.5 * erfc(Q / math.sqrt(2))                                                    # bit error rate (Eq. 13)

    return Q, SNR, BER

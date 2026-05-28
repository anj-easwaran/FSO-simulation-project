import math
import matplotlib.pyplot as plt
from parameters import *
from functions import tx_parameters, rx_parameters, channel, compute_snr
from functions import calculate_isl, calculate_updownlink
from functions import print_system_params, print_link_budget, print_snr
from plots import plot_sensitivity, plot_vs_elevation, plot_ber_vs_snr, plot_pointing_error_sensitivity, print_sensitivity_table, print_elevation_table

# Run calculations

G_T, L_T = tx_parameters(theta_T, err_T)
G_R, L_R = rx_parameters(D_R, lam_m, err_R)
d_GS, L_PS, L_PG, L_A, I_m_dB, I_g_dB = channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, theta_E_rad, L_W, N, phi)

P_R_inter,  LM_inter,  P_inter_W,  Q_inter,  SNR_inter,  BER_inter  = calculate_isl(G_T, L_T, G_R, L_R, L_PS)
P_R_updown, LM_updown, P_updown_W, Q_updown, SNR_updown, BER_updown = calculate_updownlink(G_T, L_T, G_R, L_R, L_PG, L_A)

# ---- Print Results ----

print_system_params(lam_nm, h_S, d_SS, theta_E, d_GS, G_T, G_R, L_T, L_R, L_PS, L_PG, L_A)
print_link_budget(P_R_inter, LM_inter, P_R_updown, LM_updown)
print_snr(B, r, Q_inter, SNR_inter, BER_inter, Q_updown, SNR_updown, BER_updown)

# ---- Generate plots and print result arrays ----

plot_ber_vs_snr(SNR_inter, SNR_updown, BER_inter, BER_updown)
plot_pointing_error_sensitivity()
P_rx_dBm, BER_inter_arr, BER_updown_arr = plot_sensitivity(P_R_inter, P_R_updown)
print_sensitivity_table(P_rx_dBm, BER_inter_arr, BER_updown_arr)

angles, BER_elev, LM_elev, Q_elev = plot_vs_elevation()
print_elevation_table(angles, BER_elev, LM_elev, Q_elev)

plt.show()                                 # display all figures at once

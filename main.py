import math
from parameters import *
from functions import tx_parameters, rx_parameters, channel, compute_snr

# ---- Run calculations ----

G_T, L_T = tx_parameters(theta_T, err_T)
G_R, L_R = rx_parameters(D_R, lam_m, err_R)
d_GS, L_PS, L_PG, L_A = channel(lam_m, lam_um, lam_nm, d_SS, h_E, h_S, h_A, R_E, theta_E_rad, L_W, N, phi)

P_R_inter  = P_T_inter  + eta_T + eta_R + G_T + G_R + L_T + L_R + L_PS
P_R_updown = P_T_updown + eta_T + eta_R + G_T + G_R + L_T + L_R + L_A + L_PG

LM_inter  = P_R_inter  - P_req
LM_updown = P_R_updown - P_req

P_inter_W  = 10 ** (P_R_inter  / 10) * 1e-3
P_updown_W = 10 ** (P_R_updown / 10) * 1e-3

Q_inter,  SNR_inter,  BER_inter  = compute_snr(P_inter_W,  PBGL_inter)
Q_updown, SNR_updown, BER_updown = compute_snr(P_updown_W, PBGL_updown)

# ---- Print Results ----

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

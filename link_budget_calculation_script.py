import math
from scipy.special import erfc
#initialising
#parameters

lam_m= 1550e-9       # wavelength (m)
lam_um= lam_m * 1e6   # wavelength (um) - used in Mie coefficient equations
lam_nm= lam_m * 1e9   # wavelength (nm) - used in geometrical scattering

theta_T= 15e-6   # TX beam divergence angle (rad)
D_R= 0.08    #receiver aperture diameter (m)
err_T= 1e-6 #TX pointing error (rad)
err_R= 1e-6 # RXpointing error

h_E = 1.0 #ground station altitude (km)
h_S = 550.0 #satellite altitude 
h_A = 20.0 #troposphere height
R_E = 6378.1 #Earth radius 
d_SS = 1000.0 #ISL distance 

theta_E     = 40.0 #elevation angle (degrees)
theta_E_rad = math.radians(theta_E)   #elevation angle (radians)

L_W = 0.5 #liquid water content in air (g/m3)
N   = 3.128e-4   #cloud droplet concentration (cm-3)
phi = 1.6  #particle size coefficient

P_T_inter  = 15.32 #transmitted power ISL (dBm)
P_T_updown = 13.98 #transmitted power uplink/downlink 
eta_T = -0.97  # transmitter efficiency(dB)
eta_R = -0.97  #receiver efficiency (dB)
P_req = -35.5   #receiver sensitivity at 10 Gbps (dBm)

# TX and RX gains

G_T_lin = 16 / theta_T**2 #transmitter gain (linear)
G_R_lin = (D_R * math.pi / lam_m)**2  #receiver gain (linear)

G_T = 10 * math.log10(G_T_lin) #transmitter gain (dB)
G_R = 10 * math.log10(G_R_lin) #receiver gain (dB)

#pointing loss calculation

L_T_lin = math.exp(-G_T_lin * err_T**2) #TX pointing loss (linea scale)
L_R_lin = math.exp(-G_R_lin * err_R**2) #RX pointing loss 

L_T = 10 * math.log10(L_T_lin) #TX pointing loss (dB)
L_R = 10 * math.log10(L_R_lin) #RX pointing loss 

#Slant distance

R = R_E + h_E #radius from Earth centre to ground station (km)
H = h_S - h_E  #satellite height above ground station (km)

d_GS = R * (math.sqrt(((R + H) / R)**2 - math.cos(theta_E_rad)**2) - math.sin(theta_E_rad)) #slant distance calculation

#FSPL

L_PS_lin =(lam_m / (4 * math.pi * d_SS * 1e3))**2   #intersatellite FSPL (linear)
L_PG_lin =(lam_m / (4 * math.pi * d_GS * 1e3))**2   #up/downlink FSPL

L_PS = 10 * math.log10(L_PS_lin)   #intersatellite FSPL (dB)
L_PG = 10 * math.log10(L_PG_lin)   #up/downlink FSPL (dB)

#Total atmospheric loss

L_A = -0.48


# ---- Received Power (sum of all gains and losses) ----

P_R_inter  = P_T_inter  + eta_T + eta_R + G_T + G_R + L_T + L_R + L_PS
P_R_updown = P_T_updown + eta_T + eta_R + G_T + G_R + L_T + L_R + L_A + L_PG

# Link margin
LM_inter  = P_R_inter  - P_req
LM_updown = P_R_updown - P_req

# ---- SNR / Q-factor / BER (Giggenbach 2022) ----

# Physical constants
q = 1.602e-19   # elementary charge (C)
h = 6.626e-34   # Planck constant (Js)
c = 2.998e8     # speed of light (m/s)

# Receiver parameters - from Table A1, Giggenbach 2022
eta = 0.8       # quantum efficiency
R   = eta * (q * lam_m) / (h * c)   # responsivity (A/W), approx 1.0 at 1550nm
M   = 10        # APD multiplication factor
kA  = 0.2       # ionisation ratio
Idm = 2.5e-9    # dark current - multiplied portion (A)
it  = 6.6e-12   # TIA thermal noise density (A/sqrt(Hz))
r   = 10e9      # data rate (10 Gbps)
B   = r / 2     # receiver bandwidth (Hz)

PBGL_inter  = 0.0    # background light - zero for intersatellite (no atmosphere)
PBGL_updown = 50e-9  # background light for up/downlink - 50nW typical daytime (Table 2)

FA = kA * M + 2 * (1 - kA) + (kA - 1) / M   # excess noise factor (Eq. 4)

# Convert received power from dBm to Watts
P_inter_W  = 10 ** (P_R_inter  / 10) * 1e-3
P_updown_W = 10 ** (P_R_updown / 10) * 1e-3

def compute_snr(P_rx_W, PBGL):
    sigma_t  = it * math.sqrt(B)                                                          # thermal noise (Eq. 14)
    sigma_s0 = math.sqrt(B * 2 * q * M**2 * FA * (R * PBGL + Idm))                       # shot noise during "0" (Eq. 15)
    sigma_s1 = math.sqrt(B * 2 * q * (M**2 * FA * (R * (2 * P_rx_W + PBGL) + Idm)))     # shot noise during "1" (Eq. 16)
    Q   = (M * R * 2 * P_rx_W) / (math.sqrt(sigma_s0**2 + sigma_t**2) + math.sqrt(sigma_s1**2 + sigma_t**2))  # Q-factor (Eq. 17)
    SNR = Q**2                                                                             # SNR = Q squared (Appendix A.3)
    BER = 0.5 * erfc(Q / math.sqrt(2))                                                    # bit error rate (Eq. 13)
    return Q, SNR, BER

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
##print(f"Visibility: {round(V, 2)} km")
##print(f"Atmospheric loss L_A: {round(L_A, 2)} dB  (Mie: {round(10*math.log10(I_m), 2)} dB, Geometrical: {round(10*math.log10(I_g), 2)} dB)")

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

import math
from scipy.special import erfc

# ── Link Budget ──────────────────────────────────────────────────────────────

# Gains
P_T_inter  = 15.32  # transmitted power intersatellite [dBm]
P_T_updown = 18.9   # transmitted power uplink/downlink [dBm]
G_T = 108.52        # transmitter gain [dB]
G_R = 104.19        # receiver gain [dB]

# Losses
eta_T = -0.97   # transmitter efficiency [dB]
eta_R = -0.97   # receiver efficiency [dB]
L_T   = -0.31   # transmitter pointing loss [dB]
L_R   = -0.11   # receiver pointing loss [dB]
L_PS  = -258.18 # intersatellite free space path loss [dB]
L_PG  = -261.27 # up/downlink FSPL [dB]
L_A   = -0.48   # atmospheric attenuation (0 for intersatellite) [dB]

# Receiver sensitivity at 10 Gbps
P_req = -35.5   # [dBm]

P_R_inter  = P_T_inter  + eta_T + eta_R + G_T + G_R + L_T + L_R + L_PS           # received power intersatellite [dBm]
P_R_updown = P_T_updown + eta_T + eta_R + G_T + G_R + L_T + L_R + L_A + L_PG    # received power up/downlink [dBm]

LM_inter  = P_R_inter  - P_req  # link margin intersatellite [dB]
LM_updown = P_R_updown - P_req  # link margin up/downlink [dB]

print("── Link Budget ──")
print("Intersatellite:")
print(f"  Received Power: {P_R_inter:.2f} dBm")
print(f"  Link Margin:    {LM_inter:.2f} dB")

print("\nUplink/Downlink:")
print(f"  Received Power: {P_R_updown:.2f} dBm")
print(f"  Link Margin:    {LM_updown:.2f} dB")

### PRINTING ###
if LM_inter > 0:
    print(f"\nIntersatellite link is closed with a margin of {LM_inter:.2f} dB")
elif LM_inter == 0:
    print("\nIntersatellite link is perfectly closed")
else:
    print("\nIntersatellite link is not closed")

if LM_updown > 0:
    print(f"Uplink/Downlink is closed with a margin of {LM_updown:.2f} dB")
elif LM_updown == 0:
    print("Uplink/Downlink is perfectly closed")
else:
    print("Uplink/Downlink is not closed")

# ── SNR / Q-factor / BER  (Giggenbach 2022, Sensors 22, 6773) ───────────────
# Equation numbers below match the paper.

# Physical constants
q_e = 1.602e-19   # elementary charge [C]
h_p = 6.626e-34   # Planck constant [Ws²]
c_l = 2.998e8     # speed of light [m/s]

# APD / RFE parameters — Table A1 (Giggenbach 2022)
lam  = 1550e-9    # signal wavelength [m]
eta  = 0.8        # quantum efficiency
R    = eta * (q_e * lam) / (h_p * c_l)  # responsivity [A/W] ≈ 1.0  (Eq. 1)
M    = 10         # APD multiplication factor
kA   = 0.2        # ionisation coefficient ratio (hole/electron)
Idm  = 2.5e-9     # multiplied dark current [A]
it   = 2.1e-12    # TIA thermal noise current density [A/√Hz]
r    = 10e9       # data rate [bit/s]  (10 Gbps)
B    = r / 2      # receiver bandwidth [Hz]

# Background light power — 0 for intersatellite (no atmosphere),
# 50 nW for up/downlink (typical horizon value, Table 2 of paper)
PBGL_inter  = 0.0    # [W]
PBGL_updown = 50e-9  # [W]

# Excess noise factor FA — exact formula (Eq. 4)
FA = kA * M + 2 * (1 - kA) + (kA - 1) / M

def compute_snr(P_rx_W, PBGL):
    """
    Compute Q-factor, SNR and BER for an InGaAs-APD OOK receiver.

    Parameters
    ----------
    P_rx_W : float  — mean received optical power [W]
    PBGL   : float  — background light power seen by APD [W]

    Returns
    -------
    Q   : float  — quality factor              (Eq. 17)
    SNR : float  — electrical SNR = Q²         (Appendix A.3)
    BER : float  — bit error probability        (Eq. 13)
    """
    # Thermal noise current  (Eq. 14)
    sigma_t  = it * math.sqrt(B)

    # Shot noise during binary "0" — background + dark current only  (Eq. 15)
    sigma_s0 = math.sqrt(B * 2 * q_e * M**2 * FA * (R * PBGL + Idm))

    # Shot noise during binary "1" — signal + background + dark current  (Eq. 16)
    sigma_s1 = math.sqrt(B * 2 * q_e * (M**2 * FA * (R * (2 * P_rx_W + PBGL) + Idm)))

    # Q-factor  (Eq. 17)
    Q = (M * R * 2 * P_rx_W) / (math.sqrt(sigma_s0**2 + sigma_t**2)
                               + math.sqrt(sigma_s1**2 + sigma_t**2))

    # SNR = Q²  (Appendix A.3)
    SNR = Q**2

    # Bit error probability  (Eq. 13)
    BER = 0.5 * erfc(Q / math.sqrt(2))

    return Q, SNR, BER

# Convert received powers from dBm → W
P_inter_W  = 10 ** (P_R_inter  / 10) * 1e-3
P_updown_W = 10 ** (P_R_updown / 10) * 1e-3

Q_inter,  SNR_inter,  BER_inter  = compute_snr(P_inter_W,  PBGL_inter)
Q_updown, SNR_updown, BER_updown = compute_snr(P_updown_W, PBGL_updown)

print("\n── APD Receiver SNR Analysis (Giggenbach 2022) ──")
print(f"APD parameters : M = {M}, kA = {kA}, FA = {FA:.2f}, R = {R:.3f} A/W")
print(f"Bandwidth      : {B/1e9:.1f} GHz  |  Data rate: {r/1e9:.0f} Gbps")

print(f"\nIntersatellite (PBGL = {PBGL_inter*1e9:.0f} nW):")
print(f"  Q-factor : {Q_inter:.3f}")
print(f"  SNR      : {SNR_inter:.2f}  ({10*math.log10(SNR_inter):.2f} dB)")
print(f"  BER      : {BER_inter:.2e}")

print(f"\nUplink/Downlink (PBGL = {PBGL_updown*1e9:.0f} nW):")
print(f"  Q-factor : {Q_updown:.3f}")
print(f"  SNR      : {SNR_updown:.2f}  ({10*math.log10(SNR_updown):.2f} dB)")
print(f"  BER      : {BER_updown:.2e}")

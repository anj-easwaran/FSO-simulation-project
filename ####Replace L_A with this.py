####Replace L_A with this

#Atmospheric losses

# Mie scattering coefficients - wavelength in um, altitude in km (Liang Eq. 6)
coef_a = -0.000545 * lam_um**2 + 0.002  * lam_um - 0.0038
coef_b =  0.00628  * lam_um**2 - 0.0232 * lam_um + 0.00439
coef_c = -0.028    * lam_um**2 + 0.101  * lam_um - 0.18
coef_d = -0.228    * lam_um**3 + 0.922  * lam_um**2 - 1.26 * lam_um + 0.719

rho = coef_a * h_E**3 + coef_b * h_E**2 + coef_c * h_E + coef_d   # Mie extinction coefficient

I_m = math.exp(-rho / math.sin(theta_E_rad))   # Mie scattering attenuation (Liang Eq. 7)

V = 1.002 / (L_W * N)**0.6473   # visibility (km) (Liang Eq. 8)

theta_A = (3.91 / V) * (lam_nm / 550)**(-phi)   # geometrical scattering coefficient (km-1) (Liang Eq. 9)

d_A = (h_A - h_E) / math.sin(theta_E_rad)   # troposphere path distance (km)

I_g = math.exp(-theta_A * d_A)   # geometrical scattering attenuation (Liang Eq. 10)

L_A_lin = I_m * I_g                    # total atmospheric loss (linear) (Liang Eq. 11)
L_A     = 10 * math.log10(L_A_lin)    # total atmospheric loss (dB)
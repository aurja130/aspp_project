import diffusion.xs_setup as xs

#  (taken from JAEA NDC: https://wwwndc.jaea.go.jp/)

sigma_t_2_u238 = 11.98e-24          # cm^2
sigma_t_1_u238 = 5.884e-24          # cm^2
sigma_s_2_u238 = 9.3e-24            # cm^2
sigma_s_1_u238 = 2.789e-24          # cm^2
sigma_f_2_u238 = 16.8e-6 * 1e-24    # cm^2
sigma_f_1_u238 = 1.136e-24          # cm^2
sigma_a_2_u238 = sigma_t_2_u238 - sigma_s_2_u238 - sigma_f_2_u238       # cm^2
sigma_a_1_u238 = sigma_t_1_u238 - sigma_s_1_u238 - sigma_f_1_u238       # cm^2

sigma_t_2_u235 = 698.9e-24          # cm^2
sigma_t_1_u235 = 5.894e-24          # cm^2
sigma_s_2_u235 = 15.12e-24          # cm^2
sigma_s_1_u235 = 2.839e-24          # cm^2
sigma_f_2_u235 = 585.1e-24          # cm^2
sigma_f_1_u235 = 2.053e-24          # cm^2
sigma_a_2_u235 = sigma_t_2_u235 - sigma_s_2_u235 - sigma_f_2_u235       # cm^2
sigma_a_1_u235 = sigma_t_1_u235 - sigma_s_1_u235 - sigma_f_1_u235       # cm^2

sigma_t_2_o16 = 3.965e-24           # cm^2
sigma_t_1_o16 = 1.611e-24           # cm^2
sigma_s_2_o16 = 3.965e-24           # cm^2
sigma_s_1_o16 = 0.8978e-24          # cm^2
sigma_f_2_o16 = 0                   # cm^2
sigma_f_1_o16 = 0                   # cm^2
sigma_a_2_o16 = sigma_t_2_o16 - sigma_s_2_o16 - sigma_f_2_o16           # cm^2
sigma_a_1_o16 = sigma_t_1_o16 - sigma_s_1_o16 - sigma_f_1_o16           # cm^2

sigma_t_2_h1 = 30.60e-24            # cm^2
sigma_t_1_h1 = 0.6878e-24           # cm^2
sigma_s_2_h1 = 30.27e-24            # cm^2
sigma_s_1_h1 = 0.6878e-24           # cm^2
sigma_f_2_h1 = 0                    # cm^2
sigma_f_1_h1 = 0                    # cm^2
sigma_a_2_h1 = sigma_t_2_h1 - sigma_s_2_h1 - sigma_f_2_h1               # cm^2
sigma_a_1_h1 = sigma_t_1_h1 - sigma_s_1_h1 - sigma_f_1_h1               # cm^2

# taken from wikipedia (https://en.wikipedia.org/wiki/Neutron_cross_section)

sigma_s_2_u238 = 9e-24              # cm^2
sigma_s_1_u238 = 5e-24              # cm^2
sigma_f_2_u238 = 0.00002e-24        # cm^2
sigma_f_1_u238 = 0.3e-24            # cm^2
sigma_a_2_u238 = 2e-24              # cm^2
sigma_a_1_u238 = 0.07e-24            # cm^2

sigma_s_2_u235 = 10e-24             # cm^2
sigma_s_1_u235 = 4e-24              # cm^2
sigma_f_2_u235 = 583e-24            # cm^2
sigma_f_1_u235 = 1e-24              # cm^2
sigma_a_2_u235 = 99e-24             # cm^2
sigma_a_1_u235 = 0.09e-24           # cm^2

sigma_s_2_o16 = 4e-24               # cm^2
sigma_s_1_o16 = 3e-24               # cm^2
sigma_f_2_o16 = 0                   # cm^2
sigma_f_1_o16 = 0                   # cm^2
sigma_a_2_o16 = 0.0001e-24          # cm^2
sigma_a_1_o16 = 0.00000003e-24      # cm^2

sigma_s_2_h1 = 20e-24               # cm^2
sigma_s_1_h1 = 4e-24                # cm^2
sigma_f_2_h1 = 0                    # cm^2
sigma_f_1_h1 = 0                    # cm^2
sigma_a_2_h1 = 0.2e-24              # cm^2
sigma_a_1_h1 = 0.00004e-24          # cm^2


# material data

U235_perc = 5 / 100                 # %, enrichment, variable to be adjusted to find criticality
U238_perc = 100 - U235_perc         # %
UO2_density = 10                    # g/cm^3, approxim

# Calculate macroscopic cross-sec

A_c_1 = xs.calculate_macroscopic_xs_H2O(sigma_a_1_h1, sigma_a_1_o16)    # 1/cm
A_c_2 = xs.calculate_macroscopic_xs_H2O(sigma_a_2_h1, sigma_a_2_o16)    # 1/cm
F_c_1 = xs.calculate_macroscopic_xs_H2O(sigma_f_1_h1, sigma_f_1_o16)    # 1/cm
F_c_2 = xs.calculate_macroscopic_xs_H2O(sigma_f_2_h1, sigma_f_2_o16)    # 1/cm
S_c_12 = xs.calculate_macroscopic_xs_H2O(sigma_s_1_h1, sigma_s_1_o16)   # 1/cm
S_c_21 = xs.calculate_macroscopic_xs_H2O(sigma_s_2_h1, sigma_s_2_o16)   # 1/cm
D_c_1 = xs.calculate_diffusion_coefficient_H2O(A_c_1, S_c_12)           # cm
D_c_2 = xs.calculate_diffusion_coefficient_H2O(A_c_2, S_c_21)           # cm

A_f_1 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_a_1_u235, sigma_a_1_u238, sigma_a_1_o16)           # 1/cm
F_f_1 = 2.45 * xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_f_1_u235, sigma_f_1_u238, sigma_f_1_o16)    # 1/cm, 2.45 = average neutron yield per fission
F_f_2 = 2.45 * xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_f_2_u235, sigma_f_2_u238, sigma_f_2_o16)    # 1/cm, 2.45 = average neutron yield per fission
A_f_2 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_a_2_u235, sigma_a_2_u238, sigma_a_2_o16)           # 1/cm
S_f_21 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_s_1_u235, sigma_s_1_u238, sigma_s_1_o16)          # 1/cm
S_f_12 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_s_2_u235, sigma_s_2_u238, sigma_s_2_o16)          # 1/cm
D_f_1 = xs.calculate_diffusion_coefficient_UO2(U235_perc, U238_perc, A_f_1, S_f_12)                                                 # cm
D_f_2 = xs.calculate_diffusion_coefficient_UO2(U235_perc, U238_perc, A_f_2, S_f_21)                                                 # cm
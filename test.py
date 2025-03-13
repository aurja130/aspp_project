import numpy as np
import matplotlib.pyplot as plt
import criticality_search.criticality_solver as c_sol
import criticality_search.xs_setup as xs

# geometrical data

rod_length = 366                    # Westinghouse AP1000 has 12 ft of active fuel length

layout_1D = np.ones(round(rod_length))

# material data

U235_perc = 6                       # %, enrichment, variable to be adjusted to find criticality
U238_perc = 100 - U235_perc         # %
UO2_density = 10                    # g/cm^3, approximately

# nuclear data (taken from JAEA NDC: https://wwwndc.jaea.go.jp/)

chi_1 = 0.99                        # yield, unitless
chi_2 = 0.01                        # yield, unitless

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

# Calculate macroscopic cross-sections

A_c_1 = xs.calculate_macroscopic_xs_H2O(sigma_a_1_h1, sigma_a_1_o16)    # 1/cm
A_c_2 = xs.calculate_macroscopic_xs_H2O(sigma_a_2_h1, sigma_a_2_o16)    # 1/cm
F_c_1 = xs.calculate_macroscopic_xs_H2O(sigma_f_1_h1, sigma_f_1_o16)    # 1/cm
F_c_2 = xs.calculate_macroscopic_xs_H2O(sigma_f_2_h1, sigma_f_2_o16)    # 1/cm
S_c_12 = xs.calculate_macroscopic_xs_H2O(sigma_s_1_h1, sigma_s_1_o16)   # 1/cm
S_c_21 = xs.calculate_macroscopic_xs_H2O(sigma_s_2_h1, sigma_s_2_o16)   # 1/cm
D_c_1 = 1 / 3 / S_c_12  # cm
D_c_2 = 1 / 3 / S_c_21  # cm

A_f_1 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_a_1_u235, sigma_a_1_u238, sigma_a_1_o16)           # 1/cm
F_f_1 = 2.45 * xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_f_1_u235, sigma_f_1_u238, sigma_f_1_o16)    # 1/cm
F_f_2 = 2.45 * xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_f_2_u235, sigma_f_2_u238, sigma_f_2_o16)    # 1/cm
A_f_2 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_a_2_u235, sigma_a_2_u238, sigma_a_2_o16)           # 1/cm
S_f_21 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_s_1_u235, sigma_s_1_u238, sigma_s_1_o16)          # 1/cm
S_f_12 = xs.calculate_macroscopic_xs_UO2(UO2_density, U235_perc, U238_perc, sigma_s_2_u235, sigma_s_2_u238, sigma_s_2_o16)          # 1/cm
D_f_1 = 1 / 3 / S_f_12  # cm
D_f_2 = 1 / 3 / S_f_21  # cm

# assign macroscopic cross-sections to appropriate nodes of the layout

D_1 = xs.assign_xs(layout_1D, D_f_1, D_c_1)     # cm
D_2 = xs.assign_xs(layout_1D, D_f_2, D_c_2)     # cm
A_1 = xs.assign_xs(layout_1D, A_f_1, A_c_1)     # 1/cm
A_2 = xs.assign_xs(layout_1D, A_f_2, A_c_2)     # 1/cm
F_1 = xs.assign_xs(layout_1D, F_f_1, F_c_1)     # 1/cm
F_2 = xs.assign_xs(layout_1D, F_f_2, F_c_2)     # 1/cm
S_21 = xs.assign_xs(layout_1D, S_f_21, S_c_21)  # 1/cm
S_12 = xs.assign_xs(layout_1D, S_f_12, S_c_12)  # 1/cm

# calculate crticality

k, iterations = c_sol.find_k(D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12)   # unitless

print(k, iterations)
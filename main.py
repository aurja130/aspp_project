import diffusion.geometry_setup as geom
import diffusion.xs_setup as xs
import diffusion.criticality_solver as c_sol
import matplotlib.pyplot as plt

# SETUP

# geometrical data

rods_x = 5                  # unitless
rods_y = 5                  # unitless
rod_dia = 10                # actual value / resolution, unitless
rod_pitch = 50              # actual value / resolution, unitless
h = 1                       # resolution, cm

N_y, layout_1D, layout_2D = geom.generate_layout(rods_x, rods_y, rod_dia, rod_pitch)

# nuclear data: yield and macroscopic cross-sections

chi_1 = 0.99                # yield, unitless
chi_2 = 0.01                # unitless

# fuel

A_f_1 = 0.0092144           # 1/cm
F_f_1 = 0.0065697           # 1/cm
F_f_2 = 0.1312600           # 1/cm
A_f_2 = 0.0778104           # 1/cm
S_f_21 = 0.0013089          # 1/cm
S_f_12 = 0.0181930          # 1/cm
D_f_1 = 1.0933622           # cm
D_f_2 = 0.3266693           # cm

# coolant

A_c_1 = 0.0091005           # 1/cm
F_c_1 = 0.0065011           # 1/cm
F_c_2 = 0.1271300           # 1/cm
A_c_2 = 0.0749760           # 1/cm
S_c_21 = 0.0014334          # 1/cm
S_c_12 = 0.0164530          # 1/cm
D_c_1 = 1.1461055           # cm
D_c_2 = 0.3578957           # cm

# CALCULATION

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

phi_1, phi_2, k, iterations = c_sol.find_k(h, N_y, D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12)   # unitless, unitless

# RESULTS

fig, (left, right) = plt.subplots(figsize=(12, 6), ncols=2)

fast = left.matshow(phi_1.reshape(layout_2D.shape), cmap='plasma')
left.set_title('Fast neutron flux (1/cm^2)')
left.set_xlabel('geometrical dimension y (cm)')
left.set_ylabel('geometrical dimension x (cm)')
fig.colorbar(fast, ax=left)

thermal = right.matshow(phi_2.reshape(layout_2D.shape), cmap='ocean')
right.set_title('Thermal neutron flux (1/cm^2)')
right.set_xlabel('geometrical dimension y (cm)')
right.set_ylabel('geometrical dimension x (cm)')
fig.colorbar(thermal, ax=right)

fig.suptitle(f'Multiplication factor, k = {k}')

plt.savefig('results.png')
plt.show()
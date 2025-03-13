import numpy as np
import matplotlib.pyplot as plt
import criticality_search.setup as setup
import criticality_search.solver as solver

# geometrical data

rod_diameter = 10 # mm
rod_pitch = 15 # mm
number_of_rods_x = int(1)
number_of_rods_y = int(1)

pool_length_x = int(100) # mm
pool_length_y = int(100) # mm

# generate layout

layout_2D = setup.generate_layout(rod_diameter, rod_pitch, number_of_rods_x, number_of_rods_y, pool_length_x, pool_length_y)

plt.matshow(layout_2D)
plt.show()

N_x = int(layout_2D.shape[1])

layout_flattened = layout_2D.flatten()

# nuclear data

D_c_1 = 1.1012311
D_c_2 = 0.3332974
A_c_1 = 0.0091975
A_c_2 = 0.0767878
F_c_1 = 0.0065595
F_c_2 = 0.1295450
S_c_21 = 0.0013672
S_c_12 = 0.0179990

D_f_1 = 1.0933622
D_f_2 = 0.3266693
A_f_1 = 0.0092144
A_f_2 = 0.0778104
F_f_1 = 0.0065697
F_f_2 = 0.1312600
S_f_21 = 0.0013089
S_f_12 = 0.0181930

chi_1 = 0.95
chi_2 = 0.05

D_1 = setup.assign_xs(layout_flattened, D_f_1, D_c_1)
D_2 = setup.assign_xs(layout_flattened, D_f_2, D_c_2)
A_1 = setup.assign_xs(layout_flattened, A_f_1, A_c_1)
A_2 = setup.assign_xs(layout_flattened, A_f_2, A_c_2)
F_1 = setup.assign_xs(layout_flattened, F_f_1, F_c_1)
F_2 = setup.assign_xs(layout_flattened, F_f_2, F_c_2)
S_21 = setup.assign_xs(layout_flattened, S_f_21, S_c_21)
S_12 = setup.assign_xs(layout_flattened, S_f_12, S_c_12)

k = solver.find_k(N_x, D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12)

print(k)
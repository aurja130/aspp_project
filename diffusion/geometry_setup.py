import numpy as np
import matplotlib.pyplot as plt

def generate_layout(rods_x, rods_y, rod_dia, rod_pitch):

    one_rod = np.zeros((round(rod_pitch + 1), round(rod_pitch + 1)))
    rod_center_x = rod_pitch / 2
    rod_center_y = rod_pitch / 2
    
    for i in range(one_rod.shape[0]):
        for j in range(one_rod.shape[1]):
            if (i - rod_center_x)**2 + (j - rod_center_y)**2 <= (rod_dia / 2)**2:
                one_rod[i, j] = 1
    
    all_rods = np.tile(one_rod, (rods_x, rods_y))

    return all_rods.shape[1], all_rods.flatten(), all_rods
import numpy as np

def generate_layout(rod_diameter, rod_pitch, number_of_rods_x, number_of_rods_y, pool_length_x, pool_length_y):

    # Create a 2D array of zeros
    layout = np.zeros((pool_length_y + 1, pool_length_x + 1)) # 1 point every 1 mm

    rod_area_length_x =  (number_of_rods_x - 1) * rod_pitch + rod_diameter
    rod_area_length_y =  (number_of_rods_y - 1) * rod_pitch + rod_diameter
    initial_distance_x = (pool_length_x - rod_area_length_x) / 2 + rod_diameter / 2
    initial_distance_y = (pool_length_y - rod_area_length_y) / 2 + rod_diameter / 2


    for rod_index_x in range(0, number_of_rods_x):
        for rod_index_y in range(0, number_of_rods_y):
            rod_center_x = initial_distance_x + rod_index_x * rod_pitch
            rod_center_y = initial_distance_y + rod_index_y * rod_pitch
            
            for i in range(int(rod_center_y - rod_diameter / 2 - 2), int(rod_center_y + rod_diameter / 2 + 3)):
                    for j in range(int(rod_center_x - rod_diameter / 2 - 2), int(rod_center_x + rod_diameter / 2 + 3)):
                        if (i - rod_center_y)**2 + (j - rod_center_x)**2 <= (rod_diameter / 2)**2:
                            layout[i, j] = 1

    return layout

def assign_xs(layout, fuel_xs, moderator_xs):
    xs = np.where(layout == 1, fuel_xs, moderator_xs)
    return xs
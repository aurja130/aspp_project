import numpy as np

def calculate_average_atomic_mass(perc_U235, perc_U238):
    M_U235 = 235.0439   # g/mol
    M_U238 = 238.0508   # g/mol

    M = 100 / (perc_U235 / M_U235 + perc_U238 / M_U238)  # g/mol

    return M

def calculate_individual_density(UO2_density, perc_U235, perc_U238):
    average_M_U = calculate_average_atomic_mass(perc_U235, perc_U238)   # g/mol
    M_O = 15.999    # g/mol
    average_M_UO2 = average_M_U + M_O * 2   # g/mol
    total_weight_percent_U = average_M_U / average_M_UO2    # 1
    average_density_U = UO2_density * total_weight_percent_U    # g/cm^3
    density_O = UO2_density * (1 - total_weight_percent_U)  # g/cm^3

    density_U235 = average_density_U * perc_U235 / 100  # g/cm^3
    density_U238 = average_density_U * perc_U238 / 100  # g/cm^3

    return density_U235, density_U238, density_O

def calculate_macroscopic_xs_UO2(UO2_density, perc_U235, perc_U238, sigma_U235, sigma_U238, sigma_O):
    density_U235, density_U238, density_O = calculate_individual_density(UO2_density, perc_U235, perc_U238) # g/cm^3

    atom_density_U235 = density_U235 * 6.022e23 / 235.0439  # atoms/cm^3
    atom_density_U238 = density_U238 * 6.022e23 / 238.0508  # atoms/cm^3
    atom_density_O = density_O * 6.022e23 * 2 / 15.999  # atoms/cm^3

    total_macroscopic_xs = atom_density_U235 * sigma_U235 + atom_density_U238 * sigma_U238 + atom_density_O * sigma_O   # 1/cm

    return total_macroscopic_xs

def assign_xs(layout, fuel_xs, moderator_xs):
    xs = np.where(layout == 1, fuel_xs, moderator_xs)
    return xs

def calculate_macroscopic_xs_H2O(sigma_h, sigma_o):
    # molar_mass = 18 # g/mole
    # density = 1 # g/cm^3
    # N_avogadro = 6.022e23   # molecules/mole
    # molecular_number_density = density * N_avogadro / molar_mass    # molecules/cm^3

    # N_H = 2 * molecular_number_density  # atomic number density of hydrogen, 6.691e+22 atoms/cm^3
    # N_O = 1 * molecular_number_density  # atomic number density of oxygen, 3.346e+22 atoms/cm^3

    N_H = 6.691e+22
    N_O = 3.346e+22

    SIGMA_H2O = N_H * sigma_h + N_O * sigma_o

    return SIGMA_H2O
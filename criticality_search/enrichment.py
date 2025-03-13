def calculate_average_atomic_mass(perc_U235, perc_U238):
    """
    Calculate the average atomic mass of a uranium sample given the percentage of U-235 and U-238 isotopes.

    Parameters
    ----------
    perc_U235 : float
        The percentage of U-235 in the uranium sample.
    perc_U238 : float
        The percentage of U-238 in the uranium sample.

    Returns
    -------
    float
        The total atomic mass of the uranium sample.
    """
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

def calculate_macroscopic_cross_section(UO2_density, perc_U235, perc_U238, sigma_U235, sigma_U238, sigma_O):
    density_U235, density_U238, density_O = calculate_individual_density(UO2_density, perc_U235, perc_U238) # g/cm^3

    atom_density_U235 = density_U235 * 6.022e23 / 235.0439  # atoms/cm^3
    atom_density_U238 = density_U238 * 6.022e23 / 238.0508  # atoms/cm^3
    atom_density_O = density_O * 6.022e23 * 2 / 15.999  # atoms/cm^3

    total_macroscopic_xs = atom_density_U235 * sigma_U235 + atom_density_U238 * sigma_U238 + atom_density_O * sigma_O   # 1/cm

    return total_macroscopic_xs

SIGMA_F = calculate_macroscopic_cross_section(10, 5, 95, 585e-24, 16.8e-30, 0)

print(SIGMA_F)
# 0.6605601191478466
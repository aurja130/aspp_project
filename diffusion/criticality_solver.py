import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splg
import copy
from line_profiler import profile

def make_M_1D(h, D, A, S):
    off_diagonal = D[:-1] / 2 / h**2 + D[1:] / 2 / h**2
    main_diagonal = A + S + np.concatenate(([0], off_diagonal)) + np.concatenate((off_diagonal, [0]))

    M = sp.diags_array([- off_diagonal, main_diagonal, - off_diagonal], offsets=[-1, 0, 1], format='csr')    # 1/cm
    return M

def make_M_2D(h, N_y, D, A, S):
    off_diagonal_1 = D[:-N_y] / 2 / h**2 + D[N_y:] / 2 / h**2       # 1/cm
    off_diagonal_2 = D[:-1] / 2 / h**2 + D[1:] / 2 / h**2           # 1/cm
    off_diagonal_2[N_y-1::N_y] = 0
    
    main_diagonal = (A + S + np.pad(off_diagonal_1, (N_y, 0), mode='constant', constant_values=(0, 0))      # 1/cm
                     + np.pad(off_diagonal_1, (0, N_y), mode='constant', constant_values=(0, 0))
                     + np.pad(off_diagonal_2, (1, 0), mode='constant', constant_values=(0, 0))
                     + np.pad(off_diagonal_2, (0, 1), mode='constant', constant_values=(0, 0)))

    M = sp.diags_array([- off_diagonal_1, - off_diagonal_2, main_diagonal, - off_diagonal_2, - off_diagonal_1], offsets=[-N_y, -1, 0, 1, N_y], format='csr')    # 1/cm
    return M

@profile
def find_k(h, N_y, D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12):
    M_1 = make_M_2D(h, N_y, D_1, A_1, S_12)     # 1/cm
    M_2 = make_M_2D(h, N_y, D_2, A_2, S_21)     # 1/cm

    F_1 = sp.diags(F_1, format='csr')           # 1/cm
    F_2 = sp.diags(F_2, format='csr')           # 1/cm
    S_21 = sp.diags(S_21, format='csr')         # 1/cm
    S_12 = sp.diags(S_12, format='csr')         # 1/cm

    # initial guess
    k_old = 1.0    # unitless
    phi_1_old = np.ones(D_1.shape[0])           # 1/cm^2
    phi_2_old = np.ones(D_2.shape[0])           # 1/cm^2

    iterations = 0

    while True:

        fission_rate_old = F_1 * phi_1_old + F_2 * phi_2_old                # 1/cm^3

        b_1 = (1 / k_old) * chi_1 * fission_rate_old + S_21 * phi_2_old     # 1/cm^3
        b_2 = (1 / k_old) * chi_2 * fission_rate_old + S_12 * phi_1_old     # 1/cm^3

        phi_1_new = splg.spsolve(M_1, b_1)                                  # 1/cm^2
        phi_2_new = splg.spsolve(M_2, b_2)                                  # 1/cm^2

        fission_rate_new = F_1 * phi_1_new + F_2 * phi_2_new                # 1/cm^3

        k_new = k_old * np.dot(fission_rate_new, fission_rate_new) / np.dot(fission_rate_old, fission_rate_new)     # unitless

        iterations += 1

        if np.abs((k_new - k_old) / k_old) < 1e-6:
            break
        else:
            k_old = copy.deepcopy(k_new)
            phi_1_old = copy.deepcopy(phi_1_new)
            phi_2_old = copy.deepcopy(phi_2_new)
            continue

    return phi_1_new, phi_2_new, k_new, iterations
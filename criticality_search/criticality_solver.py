import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splg
import copy

def make_M(D, A):
    off_diagonal_2 = np.zeros(D.shape[0] - 1)
    main_diagonal = np.zeros(D.shape[0])
    off_diagonal_3 = np.zeros(D.shape[0] - 1)

    for k in range(0, D.shape[0]):
        if k == 0:  # upper left corner
            main_diagonal[k] = A[k] + (D[k] + D[k + 1]) / 2 / 1**2    # 1/cm + (cm + cm)/cm^2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2 / 1**2    # 1/cm
            continue
        elif k == D.shape[0] - 1: # lower right corner
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2 / 1**2    # 1/cm
            main_diagonal[k] = A[k] + (D[k] + D[k - 1]) / 2 / 1**2    # 1/cm
            continue
        else:   # interior
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2 / 1**2    # 1/cm
            main_diagonal[k] = A[k] + (D[k] + D[k - 1]) / 2 / 1**2 + (D[k] + D[k + 1]) / 2 / 0.1**2    # 1/cm
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2 / 1**2    # 1/cm
            continue

    M = sp.diags_array([off_diagonal_2, main_diagonal, off_diagonal_3], offsets=[-1, 0, 1], format='csr')    # 1/cm
    return M

def find_k(D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12):
    M_1 = make_M(D_1, A_1)    # 1/cm
    M_2 = make_M(D_2, A_2)    # 1/cm

    F_1 = sp.diags(F_1, format='csr')    # 1/cm
    F_2 = sp.diags(F_2, format='csr')    # 1/cm
    S_21 = sp.diags(S_21, format='csr')    # 1/cm
    S_12 = sp.diags(S_12, format='csr')    # 1/cm

    # initial guess
    k_old = 1.0    # unitless
    phi_1_old = np.ones(D_1.shape[0])    # 1/cm^2
    phi_2_old = np.ones(D_2.shape[0])    # 1/cm^2

    iterations = 0

    while True:

        b_1 = (1 / k_old) * chi_1 * (F_1 * phi_1_old + F_2 * phi_2_old) + S_21 * phi_2_old
        b_2 = (1 / k_old) * chi_2 * (F_1 * phi_1_old + F_2 * phi_2_old) + S_12 * phi_1_old

        phi_1_new = splg.spsolve(M_1, b_1)
        phi_2_new = splg.spsolve(M_2, b_2)

        k_new = k_old * np.dot(F_1 * phi_1_new + F_2 * phi_2_new, F_1 * phi_1_new + F_2 * phi_2_new) / np.dot(F_1 * phi_1_old + F_2 * phi_2_old, F_1 * phi_1_new + F_2 * phi_2_new)

        iterations += 1

        if np.abs((k_new - k_old) / k_old) < 1e-6:
            break
        else:
            k_old = copy.deepcopy(k_new)
            phi_1_old = copy.deepcopy(phi_1_new)
            phi_2_old = copy.deepcopy(phi_2_new)
            continue

    return k_new, iterations
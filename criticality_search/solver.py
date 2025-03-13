import numpy as np
import scipy as sp
import scipy.sparse.linalg as splg
import copy
import matplotlib.pyplot as plt
import pandas as pd

def make_M(N_x, D, A):
    off_diagonal_1 = np.zeros(D.shape[0] - N_x)
    off_diagonal_2 = np.zeros(D.shape[0] - 1)
    main_diagonal = np.zeros(D.shape[0])
    off_diagonal_3 = np.zeros(D.shape[0] - 1)
    off_diagonal_4 = np.zeros(D.shape[0] - N_x)

    for k in range(0, D.shape[0]):
        if k == 0:  # upper left corner
            main_diagonal[k] = A[k] + (D[k] + D[k + 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue
        elif k > 0 and k < N_x - 1: # upper edge
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - 1]) / 2 + (D[k] + D[k + 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue
        elif k == N_x - 1: # upper right corner
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue
        elif k >= N_x and k < D.shape[0] - N_x and k % N_x == 0: # left edge
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k + 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue
        elif k >= N_x and k < D.shape[0] - N_x and (k + 1) % N_x == 0: # right edge
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k - 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue
        elif k == D.shape[0] - N_x: # lower left corner
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k + 1]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            continue
        elif k > D.shape[0] - N_x and k < D.shape[0] - 1: # lower edge
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k - 1]) / 2 + (D[k] + D[k + 1]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            continue
        elif k == D.shape[0] - 1: # lower right corner
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k - 1]) / 2
            continue
        else:   # interior
            off_diagonal_1[k - N_x] = - (D[k] + D[k - N_x]) / 2
            off_diagonal_2[k - 1] = - (D[k] + D[k - 1]) / 2
            main_diagonal[k] = A[k] + (D[k] + D[k - N_x]) / 2 + (D[k] + D[k - 1]) / 2 + (D[k] + D[k + 1]) / 2 + (D[k] + D[k + N_x]) / 2
            off_diagonal_3[k] = - (D[k] + D[k + 1]) / 2
            off_diagonal_4[k] = - (D[k] + D[k + N_x]) / 2
            continue

    M = sp.sparse.diags_array([off_diagonal_1, off_diagonal_2, main_diagonal, off_diagonal_3, off_diagonal_4], offsets=[-N_x, -1, 0, 1, N_x], format='csr')
    return M

def find_k(N_x, D_1, D_2, A_1, A_2, chi_1, chi_2, F_1, F_2, S_21, S_12):
    M_1 = make_M(N_x, D_1, A_1)
    M_2 = make_M(N_x, D_2, A_2)

    F_1 = sp.sparse.diags(F_1, format='csr')
    F_2 = sp.sparse.diags(F_2, format='csr')
    S_21 = sp.sparse.diags(S_21, format='csr')
    S_12 = sp.sparse.diags(S_12, format='csr')

    # initial guess
    k_old = 1.0
    phi_1_old = np.ones(D_1.shape[0])
    phi_2_old = np.ones(D_2.shape[0])

    iterations = 0

    while True:

        b_1 = (1 / k_old) * chi_1 * (F_1 * phi_1_old + F_2 * phi_2_old) + S_21 * phi_2_old
        b_2 = (1 / k_old) * chi_2 * (F_1 * phi_1_old + F_2 * phi_2_old) + S_12 * phi_1_old

        phi_1_new = splg.spsolve(M_1, b_1)
        phi_2_new = splg.spsolve(M_2, b_2)

        k_new = k_old * np.dot(F_1 * phi_1_new + F_2 * phi_2_new, F_1 * phi_1_new + F_2 * phi_2_new) / np.dot(F_1 * phi_1_old + F_2 * phi_2_old, F_1 * phi_1_new + F_2 * phi_2_new)

        iterations += 1

        if np.abs((k_new - k_old) / k_old) < 1e-9:
            break
        else:
            k_old = copy.deepcopy(k_new)
            phi_1_old = copy.deepcopy(phi_1_new)
            phi_2_old = copy.deepcopy(phi_2_new)
            continue

    return k_new
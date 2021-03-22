#!/usr/bin/env python3.7

import fpylll
import timeit
import sys

import numpy as np
import sympy as sp
import scipy as sc

from degree_32 import degree_32_in_out as deg32_io
from degree_31 import degree_31_in_out as deg31_io
from degree_30 import degree_30_in_out as deg30_io
from degree_29 import degree_29_in_out as deg29_io
from degree_28 import degree_28_in_out as deg28_io
from degree_27 import degree_27_in_out as deg27_io
from degree_26 import degree_26_in_out as deg26_io
from degree_25 import degree_25_in_out as deg25_io
from degree_24 import degree_24_in_out as deg24_io
from degree_23 import degree_23_in_out as deg23_io
from degree_22 import degree_22_in_out as deg22_io
from degree_21 import degree_21_in_out as deg21_io
from degree_20 import degree_20_in_out as deg20_io
from degree_19 import degree_19_in_out as deg19_io
from degree_18 import degree_18_in_out as deg18_io
from degree_17 import degree_17_in_out as deg17_io
from degree_16 import degree_16_in_out as deg16_io
from degree_15 import degree_15_in_out as deg15_io
from degree_14 import degree_14_in_out as deg14_io
from degree_13 import degree_13_in_out as deg13_io
from degree_12 import degree_12_in_out as deg12_io
from degree_11 import degree_11_in_out as deg11_io
from degree_10 import degree_10_in_out as deg10_io
from degree_9 import degree_9_in_out as deg9_io
from degree_8 import degree_8_in_out as deg8_io
from degree_7 import degree_7_in_out as deg7_io
from degree_5 import degree_5_in_out as deg5_io
from degree_3 import degree_3_in_out as deg3_io

from boolean_gauss import create_input_system, iterative_column_pick, maximal_iterative_column_pick
from test_load import load_data


def test_speed():
    appr_rates = [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) for i in range(6)]
    count = (50, 50, 25, 20, 20, 20)
    counter = 0
    print("I started")
    for degree, data in [(7, deg7_io), (8, deg8_io), (9, deg9_io), (10, deg10_io), (11, deg11_io), (12, deg12_io)]:
        for io_pairs in data:
            print(counter, end='\r')
            counter += 1

            matrix, result = create_input_system(degree, 2, io_pairs)
            trans = transpose_matrix(matrix)
            
            # CVP technique
            cvp_system51 = MatrixCVP(trans, result, 5, 1)

            closest51 = cvp_system51.compute_cvp()
            cvp_best51 = cvp_system51.get_solutions(closest51)

            # CVP technique
            cvp_system21 = MatrixCVP(trans, result, 2, 1)

            closest21 = cvp_system21.compute_cvp()
            cvp_best21 = cvp_system21.get_solutions(closest21)
            
            # CVP technique
            cvp_system11 = MatrixCVP(trans, result, 1, 1)

            closest11 = cvp_system11.compute_cvp()
            cvp_best11 = cvp_system11.get_solutions(closest11)
            
            # CVP technique
            cvp_system12 = MatrixCVP(trans, result, 1, 2)

            closest12 = cvp_system12.compute_cvp()
            cvp_best12 = cvp_system12.get_solutions(closest12)
            
            # CVP technique
            cvp_system13 = MatrixCVP(trans, result, 1, 3)

            closest13 = cvp_system13.compute_cvp()
            cvp_best13 = cvp_system13.get_solutions(closest13)
            
            # add results
            appr_rates[degree - 7] = (appr_rates[degree - 7][0]  + cvp_best51[1],
                                      appr_rates[degree - 7][1]  + len(cvp_best51[0]),
                                      appr_rates[degree - 7][2]  + cvp_best51[2],
                                      appr_rates[degree - 7][3]  + cvp_best21[1],
                                      appr_rates[degree - 7][4]  + len(cvp_best21[0]),
                                      appr_rates[degree - 7][5]  + cvp_best21[2],
                                      appr_rates[degree - 7][6]  + cvp_best11[1],
                                      appr_rates[degree - 7][7]  + len(cvp_best11[0]),
                                      appr_rates[degree - 7][8]  + cvp_best11[2],
                                      appr_rates[degree - 7][9]  + cvp_best12[1],
                                      appr_rates[degree - 7][10]  + len(cvp_best12[0]),
                                      appr_rates[degree - 7][11]  + cvp_best12[2],
                                      appr_rates[degree - 7][12]  + cvp_best13[1],
                                      appr_rates[degree - 7][13] + len(cvp_best13[0]),
                                      appr_rates[degree - 7][14] + cvp_best13[2])

    for i in range(len(appr_rates)):
        appr_rates[i] = (appr_rates[i][0]  / count[i],
                         appr_rates[i][1]  / count[i],
                         appr_rates[i][2]  / count[i],
                         appr_rates[i][3]  / count[i],
                         appr_rates[i][4]  / count[i],
                         appr_rates[i][5]  / count[i],
                         appr_rates[i][6]  / count[i],
                         appr_rates[i][7]  / count[i],
                         appr_rates[i][8]  / count[i],
                         appr_rates[i][9]  / count[i],
                         appr_rates[i][10] / count[i],
                         appr_rates[i][11] / count[i],
                         appr_rates[i][12] / count[i],
                         appr_rates[i][13] / count[i],
                         appr_rates[i][14] / count[i])

    return appr_rates


def QR_decomposition_LLL(monomials, result):
    matrix = np.matrix([result] + monomials)
    trans = matrix.T

    mat_mul_trans = np.matmul(matrix, trans)

    Q, R = np.linalg.qr(mat_mul_trans)

    # TODO: edit Q to not be real-valued
    smaller_mat = 0

    # LLL + get the mapping
    smaller_mat_reduced = 0

    mapping = smaller_mat^-1 * smaller_mat_reduced


# transpose boolean/logical matrix
def transpose_matrix(matrix):
    rows = len(matrix)
    cols = len(matrix[0])

    transposed = []

    for i_col in range(cols):
        new_row = []

        for i_row in range(rows):
            new_row.append(matrix[i_row][i_col])

        transposed.append(new_row)

    return transposed


# Example matrix for monomial_value=1 and diagonal_value=1:
# Out  1 0 1 1 ... 1 0   1 0 0 0 ... 0 0 0
# V₀   1 0 0 0 ... 0 1   0 1 0 0 ... 0 0 0
# V₁   0 0 0 1 ... 1 1   0 0 1 0 ... 0 0 0
# V₂   1 1 1 1 ... 1 1   0 0 0 1 ... 0 0 0
# .    . . . . ... . .   . . . . ... . . .
# .    . . . . ... . .   0 0 0 0 ... 1 0 0
# .    . . . . ... . .   0 0 0 0 ... 0 1 0
# Vn   0 1 1 0 ... 1 1   0 0 0 0 ... 0 0 1
#
# M    2 0 0 0 ... 0 0   0 0 0 0 ... 0 0 0 
# O    0 2 0 0 ... 0 0   0 0 0 0 ... 0 0 0
# D    0 0 2 0 ... 0 0   0 0 0 0 ... 0 0 0
# U    0 0 0 2 ... 0 0   0 0 0 0 ... 0 0 0
# L    . . . . ... . .   . . . . ... . . .
# O    . . . . ... . .   0 0 0 0 ... 0 0 0
#      0 0 0 0 ... 2 0   0 0 0 0 ... 0 0 0
# 2    0 0 0 0 ... 0 2   0 0 0 0 ... 0 0 0

class MatrixLLL:
    def __init__(self, monomials, result, monomial_value=1, diagonal_value=1):
        self.no_io_pairs  = len(result)
        self.no_monomials = len(monomials)
        self.monomial_val = monomial_value
        self.diagonal_val = diagonal_value

        # matrix with 2*monomial_value on diagonal
        modulos = np.zeros((self.no_io_pairs, self.no_io_pairs), dtype=int)
        np.fill_diagonal(modulos, 2 * self.monomial_val)

        # matrix with diagonal_value on diagonal
        diagonal = np.zeros((self.no_monomials + 1, self.no_monomials + 1), dtype=int)
        np.fill_diagonal(diagonal, self.diagonal_val)

        # zero matrix
        zeros = np.zeros((self.no_io_pairs, self.no_monomials + 1), dtype=int)

        # matrix of results and monomials
        results = np.array([result] + monomials, dtype=int) * self.monomial_val

        top = np.hstack((results, diagonal))
        bottom = np.hstack((modulos, zeros))

        self.lll_matrix = np.concatenate((top, bottom))


    def lll_matrix_for_fplll(self, filename):
        before = np.get_printoptions()["threshold"]
        np.set_printoptions(threshold=np.inf)

        with open(filename, "w") as f:
            f.write(str(self.lll_matrix))

        np.set_printoptions(threshold=before)


    def compute_lll(self):
        lll_matrix_list = self.lll_matrix.tolist()
        lll_fpylll_matrix = fpylll.IntegerMatrix.from_matrix(lll_matrix_list)
        reduced = fpylll.LLL.reduction(lll_fpylll_matrix)

        return reduced


    def get_solutions(self, reduced_matrix):
        solutions = []
        for row in reduced_matrix:
            if abs(row[self.no_io_pairs]) == self.diagonal_val:
                # approximation rate
                rate = 0
                for col in range(self.no_io_pairs):
                    if (abs(row[col]) == self.monomial_val):
                        rate += 1

                abs_rate = self.no_io_pairs - rate
                rate /= self.no_io_pairs
                rate = 1 - rate

                # picked indices/monomials
                indices = []
                for col in range(self.no_io_pairs + 1, self.no_io_pairs + self.no_monomials + 1):
                    if abs(row[col]) == self.diagonal_val:
                        indices.append(col - self.no_io_pairs - 1)

                solutions.append((indices, rate, abs_rate))

        return solutions



# Example matrix for monomial_value=1 and diagonal_value=1:
# V₀   1 0 0 0 ... 0 1   1 0 0 ... 0 0 0
# V₁   0 0 0 1 ... 1 1   0 1 0 ... 0 0 0
# V₂   1 1 1 1 ... 1 1   0 0 1 ... 0 0 0
# .    . . . . ... . .   . . . ... . . .
# .    . . . . ... . .   0 0 0 ... 1 0 0
# .    . . . . ... . .   0 0 0 ... 0 1 0
# Vn   0 1 1 0 ... 1 1   0 0 0 ... 0 0 1

# Example vector for monomial_value=1 and diagonal_value=1:
# Out  1 0 1 1 ... 1 0   0 0 0 ... 0 0 0

class MatrixCVP:
    def __init__(self, monomials, result, monomial_value=1, diagonal_value=1):
        self.no_io_pairs  = len(result)
        self.no_monomials = len(monomials)
        self.monomial_val = monomial_value
        self.diagonal_val = diagonal_value

        # matrix with diagonal_value on diagonal
        diagonal = np.zeros((self.no_monomials, self.no_monomials), dtype=int)
        np.fill_diagonal(diagonal, self.diagonal_val)

        # matrix of monomials
        monos = np.array(monomials, dtype=int) * self.monomial_val

        self.cvp_matrix = np.hstack((monos, diagonal))

        # result vector
        vector = [self.monomial_val * i for i in result]

        self.cvp_vector = vector + (self.no_monomials * [0])


    def lll_matrix_for_fplll(self, filename):
        before = np.get_printoptions()["threshold"]
        np.set_printoptions(threshold=np.inf)

        with open(filename, "w") as f:
            f.write(str(self.cvp_matrix))
            f.write("\n")
            f.write(str(np.array(self.cvp_vector, dtype=int)))

        np.set_printoptions(threshold=before)


    def compute_cvp(self):
        cvp_matrix_list = self.cvp_matrix.tolist()
        cvp_fpylll_matrix = fpylll.IntegerMatrix.from_matrix(cvp_matrix_list)
        reduced_matrix = fpylll.LLL.reduction(cvp_fpylll_matrix)

        closest = fpylll.CVP.closest_vector(reduced_matrix, self.cvp_vector)

        return closest


    def get_solutions(self, closest):
        # approximation rate
        rate = 0
        for col in range(self.no_io_pairs):
            if abs(closest[col]) == self.cvp_vector[col]:
                rate += 1

        abs_rate = rate
        rate /= self.no_io_pairs

        # picked indices/monomials
        indices = []
        for col in range(self.no_io_pairs + 1, self.no_io_pairs + self.no_monomials):
            if abs(closest[col]) == self.diagonal_val:
                indices.append(col - self.no_io_pairs)

        return (indices, rate, abs_rate)



if __name__ == "__main__":
    print("TODO")

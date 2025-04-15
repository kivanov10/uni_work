# ----------------------------------------------
# PHY1038 Assignment: Linear Algebra
# URN 6534278 - CPTL Lab Group 1
# May 21st 2019
# ----------------------------------------------

import numpy as np


def determinant(mat):
    """Takes in a 3x3 matrix and returns its determinant"""
    i_det = mat[0, 0]*(mat[1, 1]*mat[2, 2] - mat[1, 2]*mat[2, 1])
    j_det = mat[0, 1]*(mat[1, 0]*mat[2, 2] - mat[1, 2]*mat[2, 0])
    k_det = mat[0, 2]*(mat[1, 0]*mat[2, 1] - mat[1, 1]*mat[2, 0])
    return i_det-j_det+k_det


def cofactor(mat):
    """Takes in a 3x3 matrix and returns its cofactor"""
    # All of the following are just the separate determinant elements needed
    det_11 = mat[1, 1]*mat[2, 2] - mat[1, 2]*mat[2, 1]
    det_12 = mat[1, 0]*mat[2, 2] - mat[1, 2]*mat[2, 0]
    det_13 = mat[1, 0]*mat[2, 1] - mat[1, 1]*mat[2, 0]
    det_21 = mat[0, 1]*mat[2, 2] - mat[0, 2]*mat[2, 1]
    det_22 = mat[0, 0]*mat[2, 2] - mat[0, 2]*mat[2, 0]
    det_23 = mat[0, 0]*mat[2, 1] - mat[0, 1]*mat[2, 0]
    det_31 = mat[0, 1]*mat[1, 2] - mat[0, 2]*mat[1, 1]
    det_32 = mat[0, 0]*mat[1, 2] - mat[0, 2]*mat[1, 0]
    det_33 = mat[0, 0]*mat[1, 1] - mat[0, 1]*mat[1, 0]
    
    return np.matrix([[1*det_11, -1*det_12,  1*det_13],
                      [-1*det_21,  1*det_22, -1*det_23],
                      [1*det_31, -1*det_32,  1*det_33]])


def inverse(mat):
    """Takes a 3x3 matrix and returns its inverse"""
    det = determinant(mat)
    cof = cofactor(mat)
    return (1/det)*np.transpose(cof)

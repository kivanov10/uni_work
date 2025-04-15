# ----------------------------------------------
# PHY1038 Assignment: Linear Algebra
# URN 6534278 - CPTL Lab Group 1
# May 21st 2019
# ----------------------------------------------

import numpy as np
import matrices as mat


print('''This program takes an input matrix in a .dat file and returns its
determinant, cofactor, the inverse and the values of x1, x2 and
x3. Make sure the matrix file is in the same folder as the .py file.\n\n''')


# Step 1:
user_matrix = input(str('Enter the name of the file you want to upload: '))
a = np.loadtxt(user_matrix)
a_det = mat.determinant(a)

print('The matrix which will be used in the calculations: \n', a)
print('\n')
print('The "determinant" function gives the result:', a_det)
print('The NumPy determinant function gives the result: ', np.linalg.det(a))
print('\n')


# Step 2:
a_cof = mat.cofactor(a)

print('This is what the cofactor is:\n', a_cof)
print('\n')


# Step 3:
inverse = mat.inverse(a)

print('The matrices inverse function:\n ', inverse)
print()
print('The NumPy inverse function gives the result:\n ', np.linalg.inv(a))
print('\n')


# Step 4:
results = np.array([[0], [3], [1]])
final = inverse*results

print('The value of x1 is equal to: ', int(final[0]))
print('The value of x2 is equal to: ', int(final[1]))
print('The value of x3 is equal to: ', int(final[2]))

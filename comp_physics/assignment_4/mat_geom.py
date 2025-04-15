#----------------------------------------------
# PHY1038 Assignment: Matrix Geometry
# URN 6534278 - CPTL Lab Group 1
# April 30th 2019
#----------------------------------------------


# Imports
import numpy as np
import matplotlib.pyplot as plt


# Functions
def vector(x, y):
    """Creates a 2D vector for use in 3D space."""
    v = np.array([[x], [y], [1]])
    return v

def rotation(vector, angle):
    """Uses a rotational matrix to rotate a vector.
    The angle argument is in degrees.
    """
    sin, cos = np.sin(np.deg2rad(angle)), np.cos(np.deg2rad(angle))
    rot_mat = np.array([[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]])
    return np.dot(rot_mat, vector)

def translation(vector, x, y):
    """Uses a translation matrix to move a vector."""
    trans_mat = np.array([[1, 0, x], [0, 1, y], [0, 0, 1]])
    return np.dot(trans_mat, vector)

def rotation_center(v1, v2, v3, angle):
    """This function translates the centre of mass of the triangle to the
    center of the coordinate system and then rotates it by a specified
    amount. After that it translates it back to the original coordinates.
    The angle argument is in degrees.
    """
    com_x = (v1[0] + v2[0] + v3[0])/3
    com_y = (v1[1] + v2[1] + v3[1])/3
    v1, v2, v3 = (translation(v1, -com_x, -com_y),
                  translation(v2, -com_x, -com_y),
                  translation(v3, -com_x, -com_y))
    v1, v2, v3 = (rotation(v1, angle),
                  rotation(v2, angle),
                  rotation(v3, angle))
    v1, v2, v3 = (translation(v1, com_x, com_y),
                  translation(v2, com_x, com_y),
                  translation(v3, com_x, com_y))
    return v1, v2, v3



# Message to user
print("""
This program uses three functions (rotation, translation and a
combination of both) to manipulate a triangle constructed by 3 vectors.
After it does its job it produces a compound plot of the original, the
rotated, the translated and the rotated at the center triangle. The latter
three have the original triangle plotted as well for comparison.
""")


# The original triangle to be used throughout this exercise
original_x = np.array([1, 2, 2, 1])
original_y = np.array([1, 1, 2, 1])

v1, v2, v3 = (vector(1, 1),
              vector(2, 1),
              vector(2, 2))


# Part 1
rot_v1, rot_v2, rot_v3 = (rotation(v1, 90),
                          rotation(v2, 90),
                          rotation(v3, 90))

# Use the rotated data points to show the new triangle
rot_x = [rot_v1[0], rot_v2[0], rot_v3[0], rot_v1[0]]
rot_y = [rot_v1[1], rot_v2[1], rot_v3[1], rot_v1[1]]

# Saving the data from rotation
np.savetxt('rotate.dat', np.c_[original_x, original_y, rot_x, rot_y],
           delimiter=' ', fmt='%0.4f', header=' x      y    rot x  rot y')

# Plotting part 1
plt.subplot(221)
plt.plot(rot_x, rot_y)
plt.plot(original_x, original_y)
plt.title('Rotational matrix comparison')


# Part 2
trans_v1, trans_v2, trans_v3 = (translation(v1, 2, 2),
                                translation(v2, 2, 2),
                                translation(v3, 2, 2))
trans_x = [trans_v1[0], trans_v2[0], trans_v3[0], trans_v1[0]]
trans_y = [trans_v1[1], trans_v2[1], trans_v3[1], trans_v1[1]]

# Saving translational data
np.savetxt('translate.dat', np.c_[original_x, original_y, trans_x, trans_y],
           delimiter=' ', fmt='%0.4f', header=' x      y    tran x  tran y')

# Plotting part 2
plt.subplot(222)
plt.plot(trans_x, trans_y)
plt.plot(original_x, original_y)
plt.title('Translation comparison')


# Part 3
rot_tran_v1, rot_tran_v2, rot_tran_v3 = rotation_center(v1, v2, v3, 90)
rot_tran_x = [rot_tran_v1[0], rot_tran_v2[0], rot_tran_v3[0], rot_tran_v1[0]]
rot_tran_y = [rot_tran_v1[1], rot_tran_v2[1], rot_tran_v3[1], rot_tran_v1[1]]

# Saving the translation_rotation data
np.savetxt('centralrotation.dat', np.c_[original_x, original_y, rot_x, rot_y],
           delimiter='   ', fmt='%0.4f',
           header=' x       y    cen_rot x  cen_rot y')

# Plotting part 3
plt.subplot(223)
plt.plot(rot_tran_x, rot_tran_y)
plt.plot(original_x, original_y)
plt.title('Translation and rotation comparison')

plt.show()

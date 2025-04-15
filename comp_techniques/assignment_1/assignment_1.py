#URN:634278
#Name: Kristian Ivanov

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# Usual values used here are  0 to 2pi and 1000 points
start = float(input('Initial finite_diff value:'))
end   = float(input('Final finite_diff value:'))
points = int(input('Number of points:'))


x_arr = np.linspace(start, end, points)
y_arr = []
h     = x_arr[1] - x_arr[0] #step calculation


def g_func(u_val):
    """
    Input:
        u_val - the t value that will be multiplied by 2 and then used in the g(2t)
                func
    Returns:
        x_val - the final t value
        y     - the y val for the g(2t) func
    """
    m = 0.6366
    x_val = u_val * 2

    m = 0.6366
    if x_val > 2*np.pi:
        while x_val > 2*np.pi:
            x_val -= 2*np.pi
    elif x_val < -2*np.pi:
        while x_val < -2*np.pi:
            x_val += 2*np.pi
    #negative region
    if x_val > np.pi and x_val <=2*np.pi:
        y = m * x_val - 1
        y -= 2 #different y-intercept for that region
        return x_val,y

    #positive region
    elif x_val >= 0 and x_val <= np.pi:
        y = -m * x_val + 1
        return x_val,y


def diag_mat(main_diag,high_diag,low_diag,ndim):
    """
    Input:
        main_diag - the diagonal of the matrix being made
        high_diag - the diagonal above the main diagonal
        low_diag  - the diagonal below the main diagonal
        ndim      - dimension of the array to be returned

    Returns:
        Diagonal matrix with the inputed values
    """
    main = np.diagflat([main_diag]*ndim)
    high = np.diagflat([high_diag]*(ndim-1),1)
    low  = np.diagflat([low_diag]*(ndim-1),-1)
    final = main+high+low
    return final


for point in x_arr:
    x_p,y_p = g_func(point)
    y_arr.append(y_p)
#     plt.scatter(x_p,y_p) # used to check if the values are reasonable
# plt.show()

finite_diff = diag_mat(-2,1,1,points) # finite difference matrix

finite_diff[0,-1] = 1 #the first and last element are the same value w(0)=w(2*pi)
finite_diff[-1,0] = 1

#the g(2t) part of the equation in matrix form
g_part = diag_mat((4*np.array(y_arr)),0,0,1)


big_mat = -finite_diff/h**2 + g_part

val,vec=np.linalg.eig(big_mat) #find eigenvalues and vectors

sort_val = np.sort(val) #sort them
low_val  = sort_val[0:5] #slice the lowest 5

#Plot making
plt.figure(figsize=(10,10))
for i in range(low_val.shape[0]):
    vec_indx = np.abs(val-low_val[i]).argmin()
    plot_y = np.array([])
    plot_y = np.append(plot_y,vec[:,vec_indx])
    plt.plot(x_arr,plot_y, label=f"Eigenvalue {i}={low_val[i]}")
plt.xlabel('t')
plt.ylabel('w(t)')
plt.legend()
plt.title('w(t) using different eigenvalues')
plt.show()

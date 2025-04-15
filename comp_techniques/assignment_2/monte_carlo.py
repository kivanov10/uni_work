#Monter Carlo approach to finding the volume of 2<=n<=15-dimensional unit spheres
#URN:6534278
import numpy as np
import math as m
import matplotlib.pyplot as plt

npts = int(input("How many sampling points:")) #1 000 000 points work well
n_arr = np.random.random([15,npts])
cnt_arr = np.zeros(13) #Count storage for later
v_func = np.zeros(13)  #Analytical result storage
vol_store = np.zeros(13) #Monte Carlo results storage
# print('n',n_arr[15])
for ndim in range(1,14): #increasing dimensions starting with 2
    print('Calculating',ndim+1, 'dimensional unit sphere...')
    v_func[ndim-1] = (m.pi**((ndim+1)/2))/m.gamma((ndim+1)/2 + 1) #analytical sol
    for j in range(npts):
        sq = 0.0   #square sum storage
        for i in range(ndim+1):
            sq += n_arr[i][j]**2 #pythagoras for multidimensions
        if sq < 1.0:
            cnt_arr[ndim-1] = cnt_arr[ndim-1]+1 #count condition
    vol_store[ndim-1] = cnt_arr[ndim-1]/npts * 2**(ndim+1) #monte carlo calc
    #multiplied by 2**(ndim+1) for every new dimension

#Plotting
plt.plot(np.arange(2,15),vol_store,label='Monte Carlo')
plt.plot(np.arange(2,15),v_func,label='Analytical')
plt.xlabel('Dimensions')
plt.ylabel('Volume')
plt.title('Comparison between analytical result and Monte Carlo')
plt.legend()
plt.show()

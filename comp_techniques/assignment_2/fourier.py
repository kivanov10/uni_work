#Fast fourier transform for a set of y values and summation for continous function
#URN:6534278
import math
import numpy as np
import matplotlib.pyplot as plt

def continous(x_points,fft_points):
    """Creates a truncated Fourier sum. Accepts linspace x_points and fast Fourier
    transform real and imaginary points to create a continuous function"""
    y_points = np.zeros(len(x_points))
    for i, x_iter in enumerate(x_points):
        y_points[i] = fft_points[0].real
        for j in range(1, m):
            a_k = (fft_points[ j] + fft_points[-j]).real
            b_k = (fft_points[-j] - fft_points[ j]).imag
            y_points[i] = y_points[i]+(a_k*np.cos(j*x_iter))+(b_k*np.sin(j*x_iter))
        y_points[i] = y_points[i]+(fft_points[-m].real)*(np.cos(m*x_iter))
    return y_points

y_val = [-0.2,-0.1, 0.3, 0.2,\
          0.4, 0.5, 0.0,-0.4,\
         -0.4,-0.2, 0.1, 0.2,\
          0.2, 0.1, 0.1,-0.1]
n     = len(y_val)
m     = n//2
x_val =  [ (i/len(y_val))*2*math.pi for i in range(len(y_val))]
x_big = np.linspace(0, 15/8*np.pi, 400)

#Fast fourier transform on the y values
fft = np.fft.fft(y_val)
#Obtainin a_k and b_k, real and imaginary parts
fft_mod = fft/n

plt.scatter(x_val,y_val,marker='o',color='red',label='Function points')
plt.plot(x_big,continous(x_big,fft_mod).real,label='Fourier summation function')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.title('Fourier Summation')
plt.show()

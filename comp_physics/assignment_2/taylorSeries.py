#----------------------------------------------
#PHY1038 Assignment: Taylor seies
#URN 6534278 - CPTL Lab Group 1
#February 19th 2019
#----------------------------------------------


#IMPORTS
import numpy as np
import math as m
import matplotlib.pyplot as plt


#FUNCTIONS
#I use a function for the derivative so it's easier to use later
def f0(n):
	f=(-np.sqrt(5)/2)**n*np.sin(-n*np.arctan(2))
	return f

def taylor(x, n ,f):
	f_taylor=(f*x**n)/m.factorial(n)
	return f_taylor


#MAIN BLOCK
#Define a range for the x axis. I chose 1000 so it can be more acurate
x=np.linspace(-2*np.pi, 2*np.pi, 1000)
#I ask the user for the order of Taylor series
n=(int(input('What order of Taylor expansion you want (interger value): ')))
#This is the function which we want to compare with (taken from the script)
y1=np.sin(x)*np.exp(-x/2)
y=np.array([]) #An empty array to append the values I get from Taylor
c=0 #I define c here so I can use it later in the for loop
for x_element in x: 
    for n_element in range (1, n+1):
        a=f0(n_element)
        b=taylor(x_element, n_element, a)
        c += b
    y=np.append(y, c)
    
    
#PLOTTING
plt.plot(x, y1) #Plot of the Taylor data
plt.plot(x, y)  #Plot of the example function
plt.show()


#FILE CREATION AND SAVING
dat=np.column_stack((x,y,y1))  # I stack the data into one array
np.savetxt('taylor.dat', dat, delimiter=' ', fmt='%8.3f')
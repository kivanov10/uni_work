#----------------------------------------------
#PHY1038 Assignment: Finite difference method
#URN 6534278 - CPTL Lab Group 1
#March 19th 2019
#----------------------------------------------

#Imports
import numpy as np #add numpy for the ln function in the integral part


#Functions
def function_x (x):      #The original function f(x)
	return 1/(x**2+3*x+2)


def der_function1 (x):	 #First derivative of f(x)
	return -(2*x+3)/(x**2+3*x+2)**2


def der_function2 (x):   #Second derivative of f(x)
    return (-2*(-3*x**2-9*x-7))/((x**2+3*x+2)**3)
            
            
def dfor(x, h):
	'''This function uses the Forward Difference method to approximate the first derivative
	of a function'''
	a=function_x(x+h)	 
	b=function_x(x)
	return (a-b)/h


def dcen1(x, h):
    '''This function uses the Centered Difference method to approximate the first derivative
    of a function'''
    a=function_x(x+h/2)
    b=function_x(x-h/2)
    return (a-b)/h


def dcen2(x, h):
    '''This function uses a derived equation from the the Centered Difference method to find
    the second derivative of f(x)'''
    a=dcen1(x+h/2,h)
    b=dcen1(x-h/2,h)
    return (a-b)/h

def trapezium(a, b, N):
    '''This function uses the trapezium rule approximation to find the integral value of f(x)
    between a and b'''
    h=(b-a)/N
    f1=function_x(a)
    f2=function_x(b)
    f_sum=0.0
    for value in range(1, N):   #sum the necessary values for the equation
        f_sum += function_x(a+value*h)
    return (h/2)*(f1+f2)+h*f_sum


def integral(a, b):
    '''This function is the derived integral of f(x)=1/(x**2+3*x+2) between a and b'''
    f_integr=np.log((b+1.0)/(b+2.0))-np.log((a+1.0)/(a+2.0))
    return f_integr


def simpson(a, b, N):
    '''This function uses the Simpson rule for the function f(x)'''
    h=(b-a)/N
    fa=function_x(a)  #To be used in the final equation
    fb=function_x(b)  #Same as for fa
    f_sum_odd, f_sum_even= 0.0, 0.0 #define the sums for the two loops
    for value_odd in range (1, N, 2):   
            f_sum_odd += function_x(a+value_odd*h)
    for value_even in range (2, N-1, 2):   #the -1 in N-1 is not necessary in this case
            f_sum_even += function_x(a+value_even*h)
    return (h/3)*(fa+fb+4*f_sum_odd+2*f_sum_even)
    

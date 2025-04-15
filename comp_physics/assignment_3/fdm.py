#----------------------------------------------
#PHY1038 Assignment: Finite difference method
#URN 6534278 - CPTL Lab Group 1
#March 19th 2019
#----------------------------------------------


#Importing the custom module
import fdm_mod as md
import scipy.integrate as spi


#Welcome message to explain the purpose of the program
print ('This program is designed to show how close different methods of approximation come close to actual values.')
print ('This first section shows the Forward difference and the Centered difference method in comparison')
print ('to the first and second derivative of f(x)=1/(x**2+3*x+2).')


#First part - derivatives with forward and centered method
#User inputs
x=float(input('Enter a value for x: '))
h=float(input('Enter a small value for h: '))


#Use of functions from module.py
actual_diff1=md.der_function1(x)  #The precise value of the first derivative
actual_diff2=md.der_function2(x)
dfor1=md.dfor(x, h)                      #The first derivative using Forward Difference method
dcen1=md.dcen1(x, h)
dcen2=md.dcen2(x, h)                      #the first derivative using Centered Difference method0


#Results print
print ('Actual value of the derivative: ', actual_diff1)
print ('Forward difference method value: ', dfor1) 
print ('The difference between this method and the actual value is: ',abs(dfor1-actual_diff1))
print ()                #Add some blank lines to improve readiblity
print ('Centered difference method value: ', dcen1) 
print ('The difference between this method and the actual value is: ',abs(dcen1-actual_diff1))
print ()
print ('Actual value of the second derivative: ', actual_diff2)
print ('Centered difference method second differential value: ', dcen2)
print ('The difference between this method and the actual value is: ',abs(dcen2-actual_diff2))
print ()


#Second part - Trapezium and Simpsons integral approximations
#Message for explaining the second part
print ('In this section the program shows how close do the Trapezium rule and Simpson rule approximation get to the actual value of the')
print ('integral of the same function f(x) between points a and b')


#User inputs
x=float(input('Enter a new value for x: '))
a=float(input('Enter the first limit of the integral a: '))
b=float(input('Enthe the second limit of the integral b: '))
N=int(input('Enter an even number of increments you would like the integral to be split into: '))


#Even N condition loop
while N%2==1:
    print ('You have to enter an even number for the increments for the program to work in this section.')
    N=int(input('Enter a valid value for increments: '))
    

#Use of functions from module.py
integral_trapezium=md.trapezium(a, b, N)
integral_simpson=md.simpson(a, b, N)
integral_actual=md.integral(a, b)
integral_scipy=spi.quad(md.function_x, a, b)


#Results print
print ('The actual value of the integral is: ', integral_actual)
print ('The trapezium rule approximation gives a value of: ', integral_trapezium)
print ('The difference with the actual value is: ', abs(integral_trapezium-integral_actual))
print ()
print ('The Simpson rule approximation gives a value of: ', integral_simpson)
print ('The difference with the actual value is: ', abs(integral_simpson-integral_actual))
print ()
print ('Scipy gives a value of: ', integral_scipy[:1])   #Show only the value of the integral approximation
print ('The difference with the actual value is: ', abs(integral_scipy[:1]-integral_actual))

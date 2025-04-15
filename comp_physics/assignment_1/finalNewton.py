#---------------------------------------------
#PHY1038 Assignment: Newton-Raphson root finder
#URN 6534278 - CPTL Lab Group 1
#February 19th 2019
#----------------------------------------------

'''Importing the SciPy'''
from scipy import optimize

'''Newton function'''
#Use a function for the polynomial, its derivative and Newton method for better
#readability and usability (especially for the SciPy part)
def f1(x0, a, b, c, d):
    f=a*x0**3+b*x0**2+c*x0+d
    return f
def f2(x0, a, b, c, d):
    fprime=a*3*x0**2+b*2*x0+c
    return fprime
def newton(x0, fo, fd):
    i=0    
    while abs(fo)>limiter:
        x=x0-(fo/fd)
        x0=x
        fo=f1(x,a,b,c,d)
        fd=f2(x,a,b,c,d)
        i=i + 1    #Counter value to show number of iterations
        if fd==0:
            print('No root could be found.')
            break
        if i==max_iter:
            print('Limit of iterations reached.')
            break
    print ('The root of the polynomial is: ',x0)
    print ('It took',i,'iterations to get the result')

'''Welcome message'''
print ('This program computes roots for 3rd degree polynomials using Newton-Raphson method\n3rd degree polynomial: a.x^3+b.x^2+c.x+d=0')

''' Hard set values'''
limiter=10**-6     #Limit given by assignment. Can be easily added to the 'newton'
                   #function to add more user interactivity. Increasing it to
                   #10**-8 would bring it closer to SciPy Newton value
max_iter=100       #A limit to how many times the Newton method can be used in the
                   #program. Can be added to the 'newton' function, as well

'''User input'''
#Obtain user values for the polynomial
a=(float(input('Enter in the value of a: ')))
b=(float(input('Enter in the value of b: ')))
c=(float(input('Enter in the value of c: ')))
d=(float(input('Enter in the value of d: ')))
x0=(float(input('Enter the first x value: ')))

'''Computing the Newton root'''
f_original=f1(x0,a,b,c,d)
f_derivative=f2(x0,a,b,c,d)
f_newton=newton(x0, f_original, f_derivative )
f_scipy=optimize.newton(f1 , x0, f2, args=(a,b,c,d), tol=10**-6)
if f_original==0:  #If the user guesses the root
    print ('The inputed value of x is a root of the polynomial')
print ('This is what the SciPy Newton function computes: ', f_scipy)
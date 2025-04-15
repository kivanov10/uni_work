#Metropolis method for creating an XOR neural network gate
#URN:6534278

import numpy as np
import math as m
import matplotlib.pyplot as plt


k=1.0 #Steepness
a=0.5 #Learning rate

def sigma(x_val):
    """Sigmoid function"""
    k=1.0
    return 1.0/(1+m.exp(-k*x_val))

def neural(weight):
    """Neural network function. Calculates node values by using the supplied
    weights. Returns the a "happiness" integer, the output values rounded and
    the error."""
    err = 0.0
    happy = 0
    outs = [0 for i in range(8)]
    for i in range(8):
        i1 = tt[i][0]
        i2 = tt[i][1]
        i3 = tt[i][2]

        h4 = sigma(weight[0]+weight[3]*i1+weight[5]*i2+weight[7]*i3)
        h5 = sigma(weight[1]+weight[4]*i1+weight[6]*i2+weight[8]*i3)
        o6 = sigma(weight[2]+weight[9]*h4+weight[10]*h5)

        err = err + 0.5*(o6-tt[i][3])**2

        outs[i] = int(round(o6))
        if outs[i] == tt[i][3]:
                happy+=1

    return happy, outs, err


def metro(beta=100.0):
    """Metropolis algorithm. Accepts different starting beta values.
    Automatically restarts if it reaches max iterations with a new set of
    weights."""
    done = False
    iter_max = 100000
    while done == False: # This is here to restart the process

        w_a = 2*np.random.random(11)-1
        w_b = w_a.copy()
        happy_hold, outs_hold, err_hold = neural(w_a)
        for i in range(iter_max): #max iters
            beta = beta + 1
            for elem in range(len(w_a)):
                #change a weight
                w_b[elem] = w_b[elem]+a*(2*np.random.random(1)-1)
            #new model
            happy_new, outs_new, err_new = neural(w_b)
            if err_new < 1e-3 or happy_new == 8: #satisfactory conditions
                print('Iterations done:',i)
                done = True
                return w_b
            err_delta = err_new-err_hold

            if err_delta < 0:
                err_hold=err_new
                w_a = w_b.copy()

            elif err_delta > 0:
                prob      = np.random.random(1)
                if prob < m.exp(-beta*err_delta):
                    err_hold=err_new
                    w_a = w_b.copy()
                else: #Revert back to previous config
                    w_b = w_a.copy()
            if i==iter_max-1:
                print('Algorithm failed, restarting with new weights.')

def checker(w):
    truth_check = 0
    for i in range(8):
        i1 = tt[i][0]
        i2 = tt[i][1]
        i3 = tt[i][2]
        h4 = sigma(w[0]+w[3]*i1+w[5]*i2+w[7]*i3)
        h5 = sigma(w[1]+w[4]*i1+w[6]*i2+w[8]*i3)
        o6 = sigma(w[2]+w[9]*h4+w[10]*h5)
        err = 0.5*(o6-tt[i][3])**2
        if round(o6) == tt[i][3]:
            truth_check += 1
        print('Algorithm values:',round(o6), 'True values:',tt[i][3], 'Error:',err)
        if truth_check == 8:
            print('All algorithm values correspond to the truth table ones')
    print('Weights:')
    print('w_14 =',w[3])
    print('w_24 =',w[5])
    print('w_34 =',w[7])
    print('w_15 =',w[4])
    print('w_25 =',w[6])
    print('w_35 =',w[8])
    print('w_45 =',w[9])
    print('w_56 =',w[10])
    print('Node biases:')
    print('w_04 =',w[0])
    print('w_05 =',w[1])
    print('w_06 =',w[2])


#               I1 I2 I3 O
tt = np.array([[0, 0, 0, 0],
               [0, 0, 1, 1],
               [0, 1, 0, 1],
               [0, 1, 1, 0],
               [1, 0, 0, 1],
               [1, 0, 1, 0],
               [1, 1, 0, 0],
               [1, 1, 1, 1]])

w_final=metro()
checker(w_final)

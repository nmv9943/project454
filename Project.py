import numpy
#import openpyxl
#import cmath
#import math
import os
import numdifftools as nd #remove eventually
from numpy.linalg import inv

#shortcut for block comment = ctrl + 4

#==============================================================================
# Input File
#==============================================================================



scenarioName = "1"

#basefolder = "/Users/ninavincent/Desktop/Git/project454/"
basefolder = os.getcwd()
database = basefolder + "/data/"
lineCSV = database + "LineData" + scenarioName + ".csv"
busCSV = database + "BusData" + scenarioName + ".csv"



#==============================================================================
# Constants
#==============================================================================


Eps = 0.001
MVAbase = 100
#something else = 20 but can't read my notes


#N = 2
#Matrix = [[0.] * N for i in range(N)]


           
#==============================================================================
# Functions about read/write files
#==============================================================================








#==============================================================================
#  Functions about creating Y matrix
#==============================================================================






#==============================================================================
#  Functions about forming power system equations
#==============================================================================

#==============================================================================
# def Pk():
#     for i in 1:N :
#         Pk[k,i] = V[k]*V[i]*(G[k,i]*cos(theta[k,i])+B[k,i]*sin(theta[k,i]))
#     return Pk
# def Qk():
#==============================================================================
 





#==============================================================================
#  Functions about Newton-Raphson
#==============================================================================
# Function that computes Jacobian
def jacobian(x,y):
    return nd.Jacobian(x,y)
    
    
    
def NewtonRaphson(f,g,x0,y0,Eps):
    delta = [[1],[1]]
    while delta > [[Eps],[Eps]]:
        fx = f(x0,y0)
        gx = g(x0,y0)
        fxgx = [[fx],[gx]]
        print(fxgx)
        J = jacobian(x0,y0)
        Jinv = inv(J)
        delta = - numpy.dot(Jinv,fxgx)
        x0 = x0 + delta
    print('Root is at: ',x0)
    print('f(x) at root is: ',f(x0))
    print('g(x) at root is: ',g(x0))
    return x0


# NR test
x0test = 4
y0test = 9
def ftest(x,y):
    return x+y-15
def gtest(x,y):
    return x*y-50
print(NewtonRaphson(ftest,gtest,x0test,y0test,Eps))
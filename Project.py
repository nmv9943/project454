import numpy
#import openpyxl
#import cmath
#import math
import os
from numpy import linalg as LA
#import numdifftools as nd
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


N = 2
Matrix = [[0.] * N for i in range(N)]
print(Matrix)

print(numpy.dot(Matrix, Matrix )  )


           
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
    return [[1,1],[y,x]]
    
    
    
def NewtonRaphson(f,g,x0,y0,Eps):
    delta = [[1],[1]]
    while abs(LA.norm(delta)) > Eps:
        fx = f(x0,y0)
        gx = g(x0,y0)
        print("fx",fx)
        print("gx",gx)
        fxgx = [fx,gx]
        print("fxgx")
        print(fxgx)
        J = jacobian(x0,y0)
        print("J")
        print(J)
        Jinv = inv(J)
        print("Jinv")
        print(Jinv)
        delta = - numpy.dot(Jinv,fxgx)
        print("delta")
        print(delta)
        x0 = x0 + delta[0]
        y0 =  y0 + delta[1]
        print("x0")
        print(x0)
    print('Root is at: ',x0)
    print('f(x,y) at root is: ',f(x0,y0))
    print('g(x,y) at root is: ',g(x0,y0))
    return (x0,y0)


# NR test
x0test = 4
y0test = 9
def ftest(x,y):
    return x+y-15
def gtest(x,y):
    return x*y-50
print(NewtonRaphson(ftest,gtest,x0test,y0test,Eps))
#this test works
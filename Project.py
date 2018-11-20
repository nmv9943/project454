import numpy
from math import cos
from math import sin
import math
import cmath
#import openpyxl
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

#these are the mismatch equations
#not sure they're right

#delta P
def dPk(V,theta,P):
    dPki = [[0.] * N for i in range(N)]
    for i in range(N):
        for k in range(N):
            dPki[k,i] = V[k]*V[i]*(G[k,i]*cos(theta[k]-theta[i])+B[k,i]*sin(theta[k]-theta[i]))-P[k]
    dP =  [sum(i) for i in dPki]
    return dP
#delta Q
def dQk(V,theta,Q):
    dQki = [[0.] * N for i in range(N)]
    for i in range(N):
        for k in range(N):
            dQki[k,i] = V[k]*V[i]*(G[k,i]*sin(theta[k]-theta[i])+B[k,i]*cos(theta[k]-theta[i]))-Q[k]
    dQ =  [sum(i) for i in dQki]
    return dQ




#==============================================================================
#  Functions about Newton-Raphson
#==============================================================================
# Function that computes Jacobian
def jacobian(theta,V):
    return [[5*V*cos(theta),5*sin(theta)],
            [5*V*sin(theta),10*V-5*cos(theta)]]


# Computes the solution of a 2-D system of equations using the Newton Raphson method
# inputs: functions P and Q, initial theta0, initial V0, tolerance (epsilon)
# outputs: solutions for theta and V
def NewtonRaphson(P,Q,theta0,V0,Eps):
    delta = [[1],[1]]   #initializing the delta matrix
    while abs(LA.norm(delta)) > Eps:
        Px = P(theta0,V0)
        Qx = Q(theta0,V0)
        PxQx = [Px,Qx]
        print(PxQx)
        J = jacobian(theta0,V0)
        print("J")
        print(J)
        Jinv = inv(J)
        print("Jinv")
        print(Jinv)
        delta = - numpy.dot(Jinv,PxQx)
        print("delta")
        print(delta)
        theta0 = theta0 + delta[0]
        V0 =  V0 + delta[1]
    print('Root is at: ',theta0)
    print('f(x,y) at root is: ',P(theta0,V0))
    print('g(x,y) at root is: ',Q(theta0,V0))
    return (theta0,V0)


# NR test
theta0test = 0
V0test = 1
def Ptest(theta,V):
    return 5*V*sin(theta)-1.689
def Qtest(theta,V):
    return 5*V**2-5*V*cos(theta)-1.843
print(NewtonRaphson(Ptest,Qtest,theta0test,V0test,Eps))
#this test works
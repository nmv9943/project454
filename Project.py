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



#==============================================================================
# Functions about read/write files
#==============================================================================

def lineRead():
    infile = open(lineCSV, 'r') 
    
    headers = [] # Extracts the column headers from the .csv file
    node_from,node_to = [],[] 
    R_values, X_values, B_values,Fmax_values = [],[],[],[]
    
    # Extracts the column headers from the LineData1.csv file before iterating 
    # through and organizing the data into the appropriate list
    grab_headers = False
    for line in infile:
        entries = line.split(',')
        if grab_headers:
            node_from.append(int(entries[0]))
            node_to.append(int(entries[1]))
            R_values.append(float(entries[2]))
            X_values.append(float(entries[3]))
            B_values.append(float(entries[4]))
            Fmax_values.append(entries[5].rstrip())
        else:
            headers = entries
            grab_headers = True
    infile.close()
    
    return(headers, node_from, node_to, R_values, X_values, B_values, Fmax_values)

def busRead():
    infile = open('BusData1.csv', 'r')
    
    headers,Bus_num, P, Q, Type, Gen, V = [],[],[],[],[],[],[]
    
    # Extracts the column headers from the BusData1.csv file before iterating 
    # through and organizing the data into the appropriate list
    grab_headers = False
    for line in infile:
        entries = line.split(',')
        if grab_headers:
            Bus_num.append(int(entries[0]))
            P.append(entries[1])
            Q.append(entries[2])
            Type.append(entries[3])
            Gen.append(entries[4])
            V.append(entries[5].rstrip())
        else:
            headers = entries
            grab_headers = True
    infile.close()    
    
    #values = [headers, Bus_num, P_MW, Q_MVAR, Type, P_gen, V_set]
    
    P_MW = {}
    Q_MVAR = {}
    bus_type = {}
    P_gen = {}
    V_set = {}
    for bus in Bus_num:
        P_MW[bus] = [P.pop(0)]
        Q_MVAR[bus] = [Q.pop(0)] 
        bus_type[bus] = [Type.pop(0)] 
        P_gen[bus] = [Gen.pop(0)] 
        V_set[bus] = [V.pop(0)]
    
    return(P_MW, Q_MVAR, bus_type, P_gen, V_set)

#==============================================================================
#  Functions about creating Y matrix. Use in conjunction with lineRead()
#==============================================================================

def admittance_matrix(headers, node_from, node_to, R_values, X_values, B_values, Fmax_values):
    
    # Gets the total number of buses
    bus_amount = max(node_to)
    
    # Forms a matrix with a number of rows equal to the amount of buses
    Y_matrix = [[] for i in range (0,bus_amount)]   
    
    # Presets each row with a number of zeros equal to the amount of buses 
    for row in Y_matrix:
        for i in range(0,bus_amount):
            row.append(0)
    
    # Iterates through the list of line parameters to form the admittance matrix Y
    for start_node in node_from:
        end_node = node_to.pop(0)
        line_R = R_values.pop(0)
        line_X = X_values.pop(0)
        shunt = B_values.pop(0)/2
        
        # Adds line admittances to the appropriate element in the list
        Y_matrix[start_node-1][end_node-1] = -1/complex(line_R,line_X)
        Y_matrix[end_node-1][start_node-1] = -1/complex(line_R,line_X)
        
        # Forms the diagonal entries (self-admittances) by summing all of 
        # admittances that terminate on the present node 
        Y_matrix[start_node-1][start_node-1] += 1/complex(line_R,line_X)
        Y_matrix[end_node-1][end_node-1] += 1/complex(line_R,line_X)
        
        # Adds shunt admittances to the diagonal entries 
        Y_matrix[start_node-1][start_node-1] += complex(0,shunt)
        Y_matrix[end_node-1][end_node-1] += complex(0,shunt)

    return(Y_matrix)

#==============================================================================
#  Functions about forming power system equations
#==============================================================================
#testing
#Ybus = [[-5j,5j],[5j,-5j]]
#thetaknown = [0]


#create P known
def Pknown(inp):
    Pknown = []
    for row in range(length(inp)):
        if ("G" in imp[4,row]) or ("D" in imp[4,row]):
            Pknown = Pknown.append(row)
    return Pknown
#create Q known



#these are the mismatch equations
#do for implicit equations only (PQ for PQ buses, P for PV buses).
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
def NewtonRaphson(P,Q,Pknown,Qknown,theta0,V0,Eps):
    delta = [[1],[1]]   #initializing the delta matrix
    while abs(LA.norm(delta)) > Eps:
        Px = P(theta0,V0,Ybus,Pknown)
        Qx = Q(theta0,V0,Ybus,Qknown)
        PxQx = [Px,Qx]
        #print(PxQx)
        J = jacobian(theta0,V0)
        #print("J")
        #print(J)
        Jinv = inv(J)
        delta = - numpy.dot(Jinv,PxQx)
        #print("delta")
        #print(delta)
        theta0 = theta0 + delta[0]
        V0 =  V0 + delta[1]
    #print('Root is at: ',theta0)
    #print('f(x,y) at root is: ',P(theta0,V0))
    #print('g(x,y) at root is: ',Q(theta0,V0))
    return (theta0,V0)


# NR test
Ybus = [[-5,5],[5,-5]]
theta0test = 0
V0test = 1
Pknown = 1.689
Qknown =  1.843

#mismatch equations, these work for 2D B only matricies
def Ptest(theta,V,Ybus,Pknown):
    n=2
    P2 =  - Pknown
    for i in range(N-1):
        P2 = P2 + Ybus[i][i-1]*V*sin(theta)
    #    print(P2)
    #P2 = Ybus[n-1][n-2]*V*sin(theta)-Pknown  #not sure if correct side of matrix, but should be same
    return P2
def Qtest(theta,V,Ybus,Qknown):
    n = 2
    Q2 =  -Ybus[n-1][n-1]*V**2  - Qknown
    for i in range(N-1):
        Q2 = Q2 -Ybus[i-1][i]*V*cos(theta)
    #Q2 = -Ybus[n-1][n-1]*V**2-Ybus[n-2][n-1]*V*cos(theta)-Qknown
    return Q2
print(NewtonRaphson(Ptest,Qtest,Pknown,Qknown,theta0test,V0test,Eps))
#this test works
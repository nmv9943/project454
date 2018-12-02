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
print(database)
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
    return([headers, node_from, node_to, R_values, X_values, B_values, Fmax_values,])



def busRead():
    infile = open(busCSV, 'r')
    
    headers,Bus_num, P_MW, Q_MVAR, Type, P_gen, V_set = [],[],[],[],[],[],[]
    
    # Extracts the column headers from the BusData1.csv file before iterating 
    # through and organizing the data into the appropriate list
    grab_headers = False
    for line in infile:
        entries = line.split(',')
        if grab_headers:
            Bus_num.append(int(entries[0]))
            P_MW.append(entries[1])
            Q_MVAR.append(entries[2])
            Type.append(entries[3])
            P_gen.append(entries[4])
            V_set.append(entries[5].rstrip())
        else:
            headers = entries
            grab_headers = True
    infile.close()    
    
    #values = [headers, Bus_num, P_MW, Q_MVAR, Type, P_gen, V_set]
    
    busDict = {}
    for bus in Bus_num:
        busDict[bus] = [P_MW.pop(0), Q_MVAR.pop(0), Type.pop(0), P_gen.pop(0), V_set.pop(0)]
    
    return busDict


#==============================================================================
#  Functions about creating Y matrix
#==============================================================================






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

# Jacobian section
def J11_same(Vk,Pk,Qk,Gkk,Bkk):
    return -Vk**2*Bkk-Qk
def J11_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*Vj*(Gkj*sin(thetakj)-Bkj*cos(thetakj))

def J12_same(Vk,Pk,Qk,Gkk,Bkk):
    return Pk/Vk-Vk*Gkk
def J12_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*(Gkj*cos(thetakj)+Bkj*sin(thetakj))

def J21_same(Vk,Pk,Qk,Gkk,Bkk):
    return Pk-Vk**2*Gkk
def J21_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return -Vk*Vj*(Gkj*cos(thetakj)+Bkj*sin(thetakj))

def J22_same(Vk,Pk,Qk,Gkk,Bkk):
    return Qk/Vk-Vk*Bkk
def J22_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*(Gkj*sin(thetakj)-Bkj*cos(thetakj))

#function that computes each jacobian part (aka J11, J21, J12, J22)
def Jpartsolve(G,B,P,Q,V,theta,J_same,J_diff,rangej,rangek):
    Jpart = [[0.] * (len(rangej)) for i in rangek]
    for j in rangej:
        j1 = j+1
        Vj = V[j]
        for k in rangek:
            k1 = k+1
            Vk = V[k1]
            Bkj = B[k1][j1]
            Gkj = G[k1][j1]
            if j==k:
                Jpart[k][k]=J_same(Vk,P[k1],Q[k1],Gkj,Bkj)
            else:
                Vj = V[j1]
                thetakj = theta[k1] - theta[j1]
                Jpart[k][j]=J_diff(Vk,Vj,Gkj,Bkj,thetakj)
    return numpy.array(Jpart)

#This function computes the Jacobian of a power system
def jacobian(G,B,P,Q,V,theta):
    rangeP = range(N - 1)
    rangePQ = range(N - m)
    J11 = Jpartsolve(G,B,P,Q,V,theta,J11_same,J11_diff,rangeP,rangeP)
    J21 = Jpartsolve(G,B,P,Q,V,theta,J21_same,J21_diff,rangeP,rangePQ)
    J12 = Jpartsolve(G,B,P,Q,V,theta,J12_same,J12_diff,rangePQ,rangeP)
    J22 = Jpartsolve(G,B,P,Q,V,theta,J22_same,J22_diff,rangePQ,rangePQ)
    J = numpy.concatenate((numpy.concatenate((J11,J12),axis=1),
                           numpy.concatenate((J21,J22),axis=1)))
    return J

#test data
N = 5 #number buses
m = 2 #number of PQ buses
Pt=[0,-0.96,-0.35,-0.16,0.24]
Qt=[0,-0.62,-0.14,-0.08,-0.35]
Vt=[1.05,1,1,1,1.02]
thetat=[0,0,0,0,0]
rangek = range(N-1)
rangej = range(N-1)
Ybus = [[2.6923-13.4115j,-1.9231+9.6154j,0,0,-0.7692+3.8462j],
        [-1.9231+9.6154j,3.6538-18.1942j,-0.9615+4.8077j,0,-0.7692+3.8462j],
        [0, -0.9615+4.8077j, 2.2115-11.0027j, -0.7692+3.8462j,-0.4808+2.4038j],
        [0,0,-0.7692+3.8462j, 1.1538-5.6742j, -0.3846+1.9231j],
        [-0.7692 + 3.8462j, -0.7692 +3.8462j, -0.4808 + 2.4038j, -0.3846 + 1.9231j,2.4038 -11.8942j]]
Gt = [list(numpy.array(x).real) for x in Ybus]
Bt = [list(numpy.array(x).imag) for x in Ybus]

def Pcomputed(G,B, V, theta):
    Pcomp = [0.] * (N)
    print(Pcomp)
    for k in range(N):
        tmp = 0
        for i in range(N):
            thetaki = theta[k]-theta[i]
            tmp = tmp + V[k]*V[i]*(G[k][i]*cos(thetaki)+B[k][i]*sin(thetaki))
        Pcomp[k] = tmp
    return numpy.array(Pcomp)

def Qcomputed(G,B,V,theta):
    Qcomp = [0.] * (N)
    print(Qcomp)
    for k in range(N):
        tmp = 0
        for i in range(N):
            thetaki = theta[k]-theta[i]
            tmp = tmp + V[k]*V[i]*(G[k][i]*sin(thetaki)-B[k][i]*cos(thetaki))
        Qcomp[k] = tmp
    return numpy.array(Qcomp)

def rad2deg(rad):
    return rad/(2*math.pi)*360

# Computes the solution of a system of equations using the Newton Raphson method
# inputs: initial P and Q, initial theta0, initial V0, tolerance (epsilon)
# outputs: solutions for theta and V
def NewtonRaphson(Ybus,P,Q,theta0,V0,Eps):
    G = [list(numpy.array(x).real) for x in Ybus]
    B = [list(numpy.array(x).imag) for x in Ybus]
    delta =[[1.] * (N-m-1)] #initializing the delta matrix
    while abs(LA.norm(delta)) > Eps:
        Pcomp = Pcomputed(G, B, V0, theta0)
        Qcomp = Qcomputed(G, B, V0, theta0)
        print("Pcomp")
        print(Pcomp)
        print("Qcomp")
        print(Qcomp)
        dP = numpy.subtract(Pcomp,numpy.array(P))
        dP = dP[numpy.array(range(1,N))] #remove the first one
        print("dP")
        print(dP)
        dQ = numpy.subtract(Qcomp,numpy.array(Q))
        dQ = dQ[numpy.array(range(1,N-m+1))] #only save m+1 to N eventually
        print("dQ")
        print(dQ)
        J = jacobian(G,B,Pcomp,Qcomp,V0,theta0)
        print("J")
        print(J)
        Jinv = inv(J)
        dPQ = numpy.concatenate((dP,dQ))
        print("dPQ")
        print(dPQ)
        delta = - numpy.dot(Jinv,dPQ)
        print("delta")
        print(delta)
        dtheta = delta[numpy.array(range(N-1))]
        dtheta = numpy.insert(dtheta, 0, 0.)
        print("dtheta")
        print(dtheta)
        theta0 = theta0 + dtheta
        print("theta0")
        print(theta0)
        print(2*N-m)
        dV = delta[numpy.array(range(N-1,2*N-m-1))] #this bad code
        print("dV")
        print(dV)
        dV = numpy.insert(dV, 0, 0.)
        print(range(N-m,N))
        dV = numpy.insert(dV, range(N-m+1,N), 0.) #eventually change to other way
        print("dV")
        print(dV)
        V0 =  V0 + dV
        print('V0')
        print(V0)
        print("Norm")
        print(LA.norm(delta))
    Pcomp = Pcomputed(G, B, V0, theta0)
    Qcomp = Qcomputed(G, B, V0, theta0)
    print("Pcomp")
    print(Pcomp)
    print("Qcomp")
    print(Qcomp)
    theta0 = rad2deg(theta0) #converting to degrees instead of radians
    return (theta0,V0)

print(NewtonRaphson(Ybus,Pt,Qt,thetat,Vt,Eps))


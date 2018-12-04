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


Eps = 0.00001
MVAbase = 100
#something else = 20 but can't read my notes


#N = 2
#Matrix = [[0.] * N for i in range(N)]
N = 5 #number of buses
#m = #number of PV buses
Nl = 4 #number of lines

#==============================================================================
# Misc functions
#==============================================================================
def rad2deg(rad):
    return rad/(2*math.pi)*360

def deg2rad(deg):
    return deg/360*2*math.pi

#==============================================================================
# Functions about read/write files
#==============================================================================

def lineRead():
    infile = open(lineCSV, 'r') 
    
    headers = [] # Extracts the column headers from the .csv file
    node_from,node_to = [],[] 
    R_values, X_values, B_values,Fmax_values, Z_values = [],[],[],[],[]
    
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
            Fmax_values.append(float(entries[5].rstrip()))
            Z_values.append(complex(float(entries[2]),float(entries[3])))
        else:
            headers = entries
            grab_headers = True
    infile.close()
    
    return(node_from, node_to, R_values, X_values, B_values, Fmax_values, Z_values)

def busRead():
    infile = open(busCSV, 'r')
    
    headers = []
    Bus_num, P, Q, Type, Gen, V = [],[],[],[],[],[]
    
    # Extracts the column headers from the .csv file before iterating through
    # and organizing the data into the appropriate list
    grab_headers = False
    for line in infile:
        entries = line.split(',')
        if grab_headers:
            Bus_num.append(int(entries[0]))
            
            # Set swing bus real power to 0.
            if entries[1] == '':
                P.append(0)
            else:
                P.append(float(entries[1]))
            
            # Set swing bus reactive power to 0.
            if entries[2] == '':
                Q.append(0)
            else:    
                Q.append(float(entries[2]))
                
            Type.append(entries[3])
            
            # Set all PQ bus real power generation to 0. 
            if entries[3] == 'G':
                Gen.append(None)
            elif entries[4] == '':
                Gen.append(0)
            else:
                Gen.append(float(entries[4]))
            
            # Preset all empty V set points to 1    
            if entries[5].rstrip() == '':
                V.append(1)
            else:
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

#Alex can you write these
def lineWrite():
    return

def busWrite():
    return


#==============================================================================
#  Functions about creating Y matrix. Use in conjunction with lineRead()
#==============================================================================

def admittance_matrix(node_from, node_to, R_values, X_values, B_values, Fmax_values, Z_values):
    
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
        line_Z = complex(line_R, line_X)
        
        # Adds line admittances to the appropriate element in the list
        Y_matrix[start_node-1][end_node-1] = -1/line_Z
        Y_matrix[end_node-1][start_node-1] = -1/line_Z
        
        # Forms the diagonal entries (self-admittances) by summing all of 
        # admittances that terminate on the present node 
        Y_matrix[start_node-1][start_node-1] += 1/line_Z
        Y_matrix[end_node-1][end_node-1] += 1/line_Z
        
        # Adds shunt admittances to the diagonal entries 
        Y_matrix[start_node-1][start_node-1] += complex(0,shunt)
        Y_matrix[end_node-1][end_node-1] += complex(0,shunt)

    return(Y_matrix)

#==============================================================================
#  Functions about forming power system equations
#==============================================================================
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

#==============================================================================
#  Functions about the Jacobian
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


#==============================================================================
#  Functions about Newton-Raphson
#==============================================================================

#test data
N = 5 #number buses
m = 2 #number of PQ buses
Pt=[0,-0.96,-0.35,-0.16,0.24]
Qt=[0,-0.62,-0.14,-0.08,-0.35]
Vt=[1.05,1,1,1,1.02]
thetat=[0,0,0,0,0]
Ybus = [[2.6923-13.4115j,-1.9231+9.6154j,0,0,-0.7692+3.8462j],
        [-1.9231+9.6154j,3.6538-18.1942j,-0.9615+4.8077j,0,-0.7692+3.8462j],
        [0, -0.9615+4.8077j, 2.2115-11.0027j, -0.7692+3.8462j,-0.4808+2.4038j],
        [0,0,-0.7692+3.8462j, 1.1538-5.6742j, -0.3846+1.9231j],
        [-0.7692 + 3.8462j, -0.7692 +3.8462j, -0.4808 + 2.4038j, -0.3846 + 1.9231j,2.4038 -11.8942j]]
Gt = [list(numpy.array(x).real) for x in Ybus]
Bt = [list(numpy.array(x).imag) for x in Ybus]

# Computes the solution of a system of equations using the Newton Raphson method
# inputs: initial P and Q, initial theta0, initial V0, tolerance (epsilon)
# outputs: solutions for theta and V

def NewtonRaphson(Ybus,P,Q,theta0,V0,Eps):
    G = [list(numpy.array(x).real) for x in Ybus]
    B = [list(numpy.array(x).imag) for x in Ybus]
    delta =[[1.] * (N-m-1)] #initializing the delta matrix
    Pcomp = Pcomputed(G, B, V0, theta0)
    Qcomp = Qcomputed(G, B, V0, theta0)
    dP = numpy.subtract(Pcomp, numpy.array(P))
    dP = dP[numpy.array(range(1, N))]  # remove the first one
    dQ = numpy.subtract(Qcomp, numpy.array(Q))
    dQ = dQ[numpy.array(range(1, N - m + 1))]  # only save m+1 to N eventually
    dPQ = numpy.concatenate((dP, dQ))
    while abs(LA.norm(dPQ)) > Eps:
        J = jacobian(G,B,Pcomp,Qcomp,V0,theta0) #should this be dP & dQ or Pcomp and Qcomp?
        Jinv = inv(J)
        delta = - numpy.dot(Jinv,dPQ)
        dtheta = delta[numpy.array(range(N-1))]
        dtheta = numpy.insert(dtheta, 0, 0.)
        theta0 = theta0 + dtheta
        dV = delta[numpy.array(range(N-1,2*N-m-1))] #this bad code
        dV = numpy.insert(dV, 0, 0.)
        dV = numpy.insert(dV, range(N-m+1,N), 0.) #eventually change to other way
        V0 =  V0 + dV
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
        dPQ = numpy.concatenate((dP, dQ))
    theta0 = rad2deg(theta0) #converting to degrees instead of radians
    return (theta0,V0,Pcomp,Qcomp)



NRsol = NewtonRaphson(Ybus,Pt,Qt,thetat,Vt,Eps)
theta_sol = NRsol[0]
V_sol = NRsol[1]
P_sol_pu = NRsol[2]
Q_sol_pu = NRsol[3]

P_sol_MVA = P_sol_pu * MVAbase
Q_sol_MVA = Q_sol_pu * MVAbase
print(P_sol_MVA)
print(Q_sol_MVA)



buses = busRead()
N = len(buses[1])
print(N)
lines=lineRead()
Nl = len(lines[0])
node_from = lines[0]
node_to = lines[1]
line_R = lines[2]
line_X = lines[3]
print(Nl)
print(line_R)
print(line_X)

#==============================================================================
#  Functions about solving the power flow on lines
#==============================================================================


def lineIflow(V,node_from,node_to,R,X):
    Iline = [[0.] * N for x in range(N)]
    for line in range(Nl):
        nfrom = node_from[line]
        nto = node_to[line]
        Z = complex(R[line],X[line])
        Iline[nfrom-1][nto-1]= (V[nto-1]-V[nfrom-1])/Z
        Iline[nto-1][nfrom-1]= -Iline[nfrom-1][nto-1]
    return Iline

Vtest = [1.05,1.045,1.01,1.0,.95,.9,1.05,1.1,1.0,1.05,1.0,1.05]

def lineSflow(Iline,V):
    Sline = [[0.] * N for x in range(N)]
    for k in range(N):
        for i in range(N):
            if k!=i:
                Sline[k][i] = V[k]*(Iline[k][i]).conjugate()
    return Sline


def LineFlowTable(V, node_from, node_to, R, X):
    lineI = lineIflow(V, node_from, node_to, R, X)
    lineS = lineSflow(lineI,V)
    LineFlows = [[0.] * 10 for x in range(Nl)]
    for line in range(Nl):
        nfrom = node_from[line]
        nto = node_to[line]
        LineFlows[line][0]= nfrom
        LineFlows[line][1]= nto
        S = lineS[nfrom-1][nto-1]
        print(S)
        print(abs(S))
        print(S.real)
        print(S.imag)
        LineFlows[line][2] = S #S value
        LineFlows[line][3] = abs(S) #S magnitude
        LineFlows[line][4] =  S.real #P value
        LineFlows[line][5] =  S.imag #Q value
    return LineFlows

node_from = lines[0]
node_to = lines[1]
R = lines[2]
X = lines[3]
LineFlowTable = LineFlowTable(Vtest,node_from, node_to, R, X)
print(LineFlowTable)
#(node_from, node_to, R_values, X_values, B_values, Fmax_values, Z_values)

#def BusPInj():
#    for k in range(N):
#        for i in range(N):
#            Sinj


def solve_all():
    #read in data files
    buses = busRead()
    lines = lineRead()
    N = len(buses[1]) #number of buses
    Nl = len(lines[1]) #number of lines
    #create new objects for lines
    node_from = lines[0]
    node_to = lines[1]
    line_R = lines[2]
    line_X = lines[3]
    line_B = lines[4]
    line_Fmax = lines[5]
    line_Z = lines[6]
    #create new objects for buses
    P = buses[0]
    Q = buses[1]
    bus_type = buses[2]
    P_gen = buses[3]
    V_set = buses[4]


    Ybus = admittance_matrix(node_from,node_to,line_R,line_X,line_B,line_Fmax,line_Z)
    NR = NewtonRaphson(Ybus, P, Q, theta0, V0, Eps)
    V = NR[1]
    LineFlowTable = LineFlowTable(V,node_from,node_to,line_R,line_X)

    #write functions
    busWrite()
    lineWrite()

solve_all()
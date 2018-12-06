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
import csv

#shortcut for block comment = ctrl + 4

#==============================================================================
# Input File
#==============================================================================

scenarioName = "1" #this is base case
#scenarioName = "2" #this is contingency case #2
#scenarioName = "3" #this is contingency case #3

#basefolder = "/Users/ninavincent/Desktop/Git/project454/"
basefolder = os.getcwd()
database = basefolder + "/data/"

#files to read in
lineCSV = database + "LineData" + scenarioName + ".csv"
busCSV = database + "BusData" + scenarioName + ".csv"
#files to write out
busWritefile = database +"BusDataOutput" + scenarioName + ".csv"
lineWritefile = database +"LineDataOutput" + scenarioName + ".csv"


#==============================================================================
# Constants
#==============================================================================


MVAbase = 100.
Eps = .1/MVAbase
#something else = 20 but can't read my notes
Vmax = 1.05
Vmin = 0.95
nround = 5 #number of decimals to round to


#N = 2
#Matrix = [[0.] * N for i in range(N)]
#N = 5 #number of buses
#m = #number of PV buses
#Nl = 4 #number of lines

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
                P.append(float(entries[1])) #change to pu
            
            # Set swing bus reactive power to 0.
            if entries[2] == '':
                Q.append(0)
            else:    
                Q.append(float(entries[2])) #change to pu
                
            Type.append(entries[3])
            
            # Set all PQ bus real power generation to 0. 
            if entries[3] == 'G':
                Gen.append(0)
            elif entries[4] == '':
                Gen.append(0)
            else:
                Gen.append(float(entries[4]))
            
            # Preset all empty V set points to 1    
            if entries[5].rstrip() == '':
                V.append(float(1))
            else:
                V.append(float(entries[5].rstrip()))
        else:
            headers = entries
            grab_headers = True
    infile.close()

    P_MW = P
    Q_MVAR = Q
    bus_type = Type
    P_gen = Gen
    V_set = V
    
    return(Bus_num,P_MW, Q_MVAR, bus_type, P_gen, V_set)

#==============================================================================
#  Functions about initializing data or things to do at the beginning
#==============================================================================

#converts power to pu units
def toPU(P_MW_gen,P_MW,Q_MVAR):
    P_gen_pu = [x / MVAbase for x in P_MW_gen]
    P_load_pu = [x / MVAbase for x in P_MW]
    Q_load_pu = [x / MVAbase for x in Q_MVAR]
    return (P_gen_pu, P_load_pu,Q_load_pu)


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
def Pcomputed(G,B, V, theta,N):
    Pcomp = [0.] * (N)
    for k in range(N):
        tmp = 0
        for i in range(N):
            thetaki = theta[k]-theta[i]
            tmp = tmp + V[k]*V[i]*(G[k][i]*cos(thetaki)+B[k][i]*sin(thetaki))
        Pcomp[k] = tmp
    return numpy.array(Pcomp)

def Qcomputed(G,B,V,theta,N):
    Qcomp = [0.] * (N)
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
    return Pk/Vk+Vk*Gkk
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
    Jpart = [[0.] * (len(rangej)) for i in range(len(rangek))]

    for j in rangej:
        Vj = V[j]
        matrix_j = rangej.index(j)

        for k in rangek:
            matrix_k = rangek.index(k)
            Vk = V[k]
            Bkj = B[k][j]
            Gkj = G[k][j]
            if j==k:
                Jpart[matrix_k][matrix_j]=J_same(Vk,P[k],Q[k],Gkj,Bkj)
            else:
                Vj = V[j]
                thetakj = theta[k] - theta[j]
                Jpart[matrix_k][matrix_j] = J_diff(Vk, Vj, Gkj, Bkj,thetakj)
    return numpy.array(Jpart)

#This function computes the Jacobian of a power system
def jacobian(G,B,P,Q,V,theta,PQPVbuses,PQbuses):
    J11 = Jpartsolve(G,B,P,Q,V,theta,J11_same,J11_diff,PQPVbuses,PQPVbuses)
    J21 = Jpartsolve(G,B,P,Q,V,theta,J21_same,J21_diff,PQPVbuses,PQbuses)
    J12 = Jpartsolve(G,B,P,Q,V,theta,J12_same,J12_diff,PQbuses,PQPVbuses)
    J22 = Jpartsolve(G,B,P,Q,V,theta,J22_same,J22_diff,PQbuses,PQbuses)
    J = numpy.concatenate((numpy.concatenate((J11,J12),axis=1),
                           numpy.concatenate((J21,J22),axis=1)))
    return J


#==============================================================================
#  Functions about Newton-Raphson
#==============================================================================

#supporting function
def MaxMismatch(dP,dQ,PQbuses,PQPVbuses):
    #print largest P mismatch and the bus
    absdP = [abs(item) for item in dP]
    maxPmismatch = max(absdP)
    i =absdP.index(maxPmismatch)
    print("Largest P mismatch: " + str(maxPmismatch) +" pu at bus: "+str(PQPVbuses[i]))
    #print largest Q mimatch and the bus
    absdQ = [abs(item) for item in dQ]
    maxQmismatch = max(absdQ)
    i =absdQ.index(maxQmismatch)
    print("Largest Q mismatch: " + str(maxQmismatch) +" pu at bus: "+str(PQbuses[i]))
    return

# Computes the solution of a system of equations using the Newton Raphson method
# inputs: initial P and Q, initial theta0, initial V0, tolerance (epsilon)
# outputs: solutions for theta and V
def NewtonRaphson(bus_types,Ybus,Pknown,Qknown,theta0,V0,Eps,N):
    G = [list(numpy.array(x).real) for x in Ybus]
    B = [list(numpy.array(x).imag) for x in Ybus]
    #which are PQ buses and PV buses?
    PQbuses = [i for i in range(len(bus_types)) if bus_types[i] == "D"] #PQ bus if only load
    PVbuses = [i for i in range(len(bus_types)) if (bus_types[i] == "G" or bus_types[i] == "GD" or bus_types[i] == "DG")] #PV bus if all but swing
    del PVbuses[0] #removing swing bus
    PQPVbuses = [i for i in range(len(bus_types)) if (bus_types[i] == "D" or bus_types[i] == "G" or bus_types[i] == "GD" or bus_types[i] == "DG")]
    del PQPVbuses[0] #removing swing bus

    m = len(PVbuses)+1 #m value (there are m-1 PVbuses)

    dPQ =numpy.array([1. for i in range(N-m-1)]) #initializing the dPQ matrix
    theta0 = numpy.array(theta0)
    V0 = numpy.array(V0)
    iteration = 1
    print('======== Newton Raphson Mismatch Information ========')
    while abs(max(dPQ.min(),dPQ.max(), key=abs)) > Eps: #find maximum mismatch to test against eps and iterate if greater than
        print("Iteration: " + str(iteration))
        Pcomp = Pcomputed(G, B, V0, theta0,N)
        Qcomp = Qcomputed(G, B, V0, theta0,N)
        #mismatches
        dP = numpy.subtract(Pcomp,numpy.array(Pknown))
        dP = dP[PQPVbuses]

        dQ = numpy.subtract(Qcomp,numpy.array(Qknown))
        dQ = dQ[PQbuses]
        MaxMismatch(dP, dQ, PQbuses, PQPVbuses) #printing max mismatch value (absolute value)

        dPQ = numpy.concatenate((dP, dQ))#mismatch vector

        #solve jacobian
        J = jacobian(G,B,Pcomp,Qcomp,V0,theta0,PQPVbuses,PQbuses) #should this be dP & dQ or Pcomp and Qcomp?
        Jinv = inv(J) #inverse jacobian
        delta = - numpy.dot(Jinv,dPQ) #correction vector
        dtheta = delta[numpy.array(range(N - 1))]
        theta0[PQPVbuses] = theta0[PQPVbuses] + dtheta

        dV = delta[numpy.array(range(N-1,2*N-m-1))]
        V0[PQbuses] =  V0[PQbuses] + dV
        iteration = iteration + 1
    [P_new, Q_new] = solveExplicit(G, B, theta0, V0, Pknown, Qknown, N)
    return (theta0,V0,P_new,Q_new)


def solveExplicit(G,B,theta,V,Pknown,Qknown,N):
    Pcomp = Pcomputed(G, B, V, theta, N)
    Qcomp = Qcomputed(G, B, V, theta, N)
    return(Pcomp,Qcomp)

def V_violation(V,N):
    Vviolation = ["No" for i in range(N)]
    for i in range(N):
        if (V[i]>Vmax or V[i]<Vmin):
            Vviolation[i] = "Yes"
    return Vviolation

#==============================================================================
#  Functions about solving the power flow on lines
#==============================================================================


def calclineflow(V_pu,theta, node_from,node_to,Z,B,Fmax,Nl):
    Pline = [0. for i in range(Nl)]
    Qline = [0. for i in range(Nl)]
    Sline = [0. for i in range(Nl)]
    Smagline = [0. for i in range(Nl)]
    Pline_MW = [0. for i in range(Nl)]
    Qline_MVAR = [0. for i in range(Nl)]
    Smag_MVA = [0. for i in range(Nl)]
    violation = ["No" for i in range(Nl)]
    for line in range(Nl):
        nfrom = node_from[line]
        nto = node_to[line]
        Vdrop = V_pu[nfrom-1]-V_pu[nto-1]
        I_line = Vdrop/Z[line]
        I_shunt = V_pu[nfrom-1]*(complex(0,B[line]/2))
        #I_total = I_line + I_shunt
        I_total = I_line
        Sline[line] = Vdrop*I_total.conjugate()
        Pline[line] = Sline[line].real
        Qline[line] = Sline[line].imag
        Smagline[line] = abs(Sline[line])
        # rescale powers to not pu
        Smag_MVA[line] = numpy.round(Smagline[line] * MVAbase,nround)
        Pline_MW[line] = numpy.round(Pline[line]* MVAbase,nround)
        Qline_MVAR[line] = numpy.round(Qline[line] * MVAbase,nround)
        if Smagline[line] > Fmax[line]:
            violation[line] = "Yes"
    return (Pline_MW,Qline_MVAR,Smag_MVA,violation)



#==============================================================================
#  Functions about writing output files
#==============================================================================
def lineWrite(listOfArrays):
    toPrint = zip(*listOfArrays)
    myFile = open(lineWritefile,'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerow(["node_from", "node_to", "P_MW", "Q_MVAR", "S_MVA", "Fmax_violation?"])
        writer.writerows(toPrint)
    return


def busWrite(listOfArrays):
    #[Bus_num,bus_type,theta,V,P_MW_new,Q_MVA_new]
    toPrint = zip(*listOfArrays)
    myFile = open(busWritefile,'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerow(["BusNumber","Type","theta (deg)","P_MW","Q_MVAR","V (pu)","V_Violation"])
        writer.writerows(toPrint)
    return

def LineFlowTable(V, theta, node_from, node_to, Z,Fmax):
    #print(Nl)
    #print(Fmax)
    LineFlows = calclineflow(V,theta, node_from,node_to,Z,Fmax)
    Pline = LineFlows[0]
    Qline = LineFlows[1]
    Sline = LineFlows[2]
    LineFlows = [[0.] * 10 for x in range(Nl)]
    for line in range(Nl):
        #print(line)
        nfrom = node_from[line]
        #print(nfrom)
        nto = node_to[line]
        LineFlows[line][0]= nfrom
        LineFlows[line][1]= nto
        S = lineS[nfrom-1][nto-1]
        LineFlows[line][2] = S #S value
        LineFlows[line][3] = abs(S) #S magnitude
        LineFlows[line][4] =  S.real #P value
        LineFlows[line][5] =  S.imag #Q value
    return LineFlows


#==============================================================================
#  Functions about writing results to the output window
#==============================================================================

def busPrint(Bus_num,bus_type,theta,P_MW_new,Q_MVAR_new,V,VViolation):
    # Printing the bus information
    print('======== Bus Information ========')    
    for i in range(max(Bus_num)):
        print('Bus ' + str(i+1))
        print('Voltage: ' + str(V[i]) + ' p.u.')
        print('Angle: ' + str(theta[i]) + ' degrees')
        if (bus_type[i] == 'G' or bus_type[i] == 'DG'):
            print('Real Power Produced: ' + str(P_MW_new[i]) + ' MW')
            print('Reactive Power Produced: ' + str(Q_MVAR_new[i]) + ' MVAR')
            print('Voltage violation? ' + VViolation[i])
        print('')
    return
        
def linePrint(node_from,node_to,Pline,Qline,Smagline,MVAviolation):
    # Printing the line information
    print('======== Line Information ========')
    for i in range(max(node_to)):
        print('')
        print('Line ' + str(node_from[i]) + ' to ' + str(node_to[i]))
        print('Active Power: ' + str(Pline[i]) + ' MW')
        print('Reactive Power: ' + str(Qline[i]) + ' MVAR')
        print('Apparent Power: ' + str(Smagline[i]) + ' MVA')
        print('Apparent Power flow violation? ' + MVAviolation[i])
    return


def solve_all():
    #read in data files
    buses = busRead()
    lines = lineRead()

    N = len(buses[1]) #number of buses
    Nl = len(lines[1]) #number of lines

    #create new objects for lines
    [node_from,node_to,line_R,line_X,line_B,line_Fmax,line_Z]=lines
    #create new objects for buses
    [Bus_num,P_MW,Q_MVAR,bus_type,P_MW_gen,V_set] = buses

    #conversion to pu units
    [P_gen_pu, P_load_pu,Q_load_pu] = toPU(P_MW_gen,P_MW,Q_MVAR)

    #initializations
    Q_gen_pu = [0. for i in range(N)] #initialize Q generation matrix
    theta = [0. for i in range(N)]  # initialize theta

    #create Ybus matrix
    Ybus = admittance_matrix(node_from.copy(), node_to.copy(), line_R.copy(), line_X.copy(), line_B.copy(), line_Fmax.copy(), line_Z.copy())

    #calculate P & Q injections
    P_inj = numpy.subtract(P_gen_pu,P_load_pu)
    Q_inj = numpy.subtract(Q_gen_pu,Q_load_pu)

    #perform Newton Raphson method to solve for implicit equations
    [theta, V, P_pu_new, Q_pu_new] = NewtonRaphson(bus_type, Ybus, P_inj, Q_inj, theta, V_set, Eps,N)

    #check for violations
    VViolation = V_violation(V,N)

    #calculate power flows on transmission lines
    [Pline, Qline, Smagline, MVAviolation] = calclineflow(V,theta, node_from,node_to,line_Z,line_B,line_Fmax,Nl)

    #LineFlowTable = LineFlowTable(V,theta,node_from,node_to,line_Z)
    theta = rad2deg(theta) #converting to degrees instead of radians

    #multiply by MVAbase so not pu anymore
    P_MW_new = numpy.round(numpy.multiply(P_pu_new,MVAbase),nround)
    Q_MVAR_new = numpy.round(numpy.multiply(Q_pu_new,MVAbase),nround)
    
    # Writing results to the output window
    busPrint(Bus_num,bus_type,theta,P_MW_new,Q_MVAR_new,V,VViolation)
    linePrint(node_from,node_to,Pline,Qline,Smagline,MVAviolation)
    
    # Writing results to two output files in the data folder
    busWrite([Bus_num,bus_type,theta,P_MW_new,Q_MVAR_new,V,VViolation]) 
    lineWrite([node_from,node_to,Pline,Qline,Smagline,MVAviolation])
    return

solve_all()

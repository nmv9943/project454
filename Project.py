# Nina Vincent, Alexander Le, and Zufan Mehari
# EE 454: Power System Analysis
# December 7, 2018
# Power Flow Project

import numpy
from math import cos
from math import sin
import math
import cmath
import os
from numpy.linalg import inv
import csv

#==============================================================================
#  Function that puts all functions together and solves
#==============================================================================

# This function puts all functions together and solves for all V,theta,P,Q of generators, P,Q,S of lines, etc.
# Writes solution to csv files / prints
def solve_all():
    #read in data files
    buses = busRead()
    lines = lineRead()

    N = len(buses[1]) #number of buses
    Nl = len(lines[1]) #number of lines

    #create new objects for lines
    [node_from,node_to,line_R,line_X,line_B,line_Fmax,line_Z]=lines
    #create new objects for buses
    [Bus_num,P_MW_load,Q_MVAR_load,bus_type,P_MW_gen,V_set] = buses

    #conversion to pu units
    [P_gen_pu, P_load_pu,Q_load_pu] = toPU(P_MW_gen,P_MW_load,Q_MVAR_load)

    #initializations
    Q_gen_pu = [0. for i in range(N)] #initialize Q generation matrix
    theta = [0. for i in range(N)]  # initialize theta

    #create Ybus matrix (note: I have to copy each because otherwise it deletes certain values within this function)
    Ybus = admittance_matrix(node_from.copy(), node_to.copy(), line_B.copy(), line_Z.copy())

    # calculate (initial) P & Q injections
    P_inj = numpy.subtract(P_gen_pu,P_load_pu)
    Q_inj = numpy.subtract(Q_gen_pu,Q_load_pu)

    # perform Newton Raphson method to solve the power flow equations
    [theta, V, P_pu_new, Q_pu_new] = NewtonRaphson(bus_type, Ybus, P_inj, Q_inj, theta, V_set, Eps,N)

    # check for voltage violations
    VViolation = V_violation(V,N)

    # calculate power flows on transmission lines and violations
    [Smagline_from2to,Pline_from2to,Qline_from2to,Smagline_to2from,Pline_to2from,Qline_to2from,MVAviolation] = calc_LineFlow(V,theta, node_from,node_to,line_Fmax,Nl,line_Z,line_B)

    #converting to degrees instead of radians
    theta = rad2deg(theta)

    # multiply by MVAbase so not pu anymore
    P_MW = numpy.multiply(P_pu_new,MVAbase)
    Q_MVAR = numpy.multiply(Q_pu_new,MVAbase)

    #solving for generator P & Q:
    P_MW_gen_solved = numpy.add(P_MW,P_MW_load)
    Q_MVAR_gen_solved = numpy.add(Q_MVAR,Q_MVAR_load)

    #rounds these values
    [Smagline_from2to,Pline_from2to,Qline_from2to,Smagline_to2from,Pline_to2from,Qline_to2from,V,theta,P_MW_gen_solved,Q_MVAR_gen_solved]=roundAll(Smagline_from2to,Pline_from2to,Qline_from2to,Smagline_to2from,Pline_to2from,Qline_to2from,V,theta,P_MW_gen_solved,Q_MVAR_gen_solved)

    # Writing results to the output window
    busPrint(Bus_num,bus_type,theta,P_MW_gen_solved,Q_MVAR_gen_solved,V,VViolation)
    linePrint(node_from,node_to,Smagline_from2to,MVAviolation)
    
    # Writing results to two output files in the data folder
    busWrite([Bus_num,bus_type,theta,P_MW_gen_solved,Q_MVAR_gen_solved,V,VViolation])
    lineWrite([node_from,node_to,Pline_from2to,Qline_from2to,Smagline_from2to,Pline_to2from,Qline_to2from,Smagline_to2from,MVAviolation])
    return

# actually running everything
solve_all()

#==============================================================================
# Input File
#==============================================================================

scenarioName = "1" #this is base case
#scenarioName = "2" #this is contingency case #2
#scenarioName = "3" #this is contingency case #3

#setting base folder
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

MVAbase = 100. #this is the MVA base
Eps = .1/MVAbase #this is epsilon (aka, tolerance of convergence for Newton Raphson (NR))

#eventually will check that V is within this range (pu) at each bus
Vmax = 1.05
Vmin = 0.95

#number of decimals to round to
nround = 2
Vround = 3


#==============================================================================
# Misc functions
#==============================================================================

#converts radians to degrees
def rad2deg(rad):
    return rad/(2*math.pi)*360

#converts degrees to radians
def deg2rad(deg):
    return deg/360*2*math.pi

#==============================================================================
# Functions about read/write files
#==============================================================================

# This function reads in the line data from the corresponding .csv file and writes out arrays with:
# node_from (the start bus of the line), node_to (the end bus of the line),
# R_values (the pu values of series resistance), X_values (the pu values of series reactance),
# B_values (the pu values of shunt charging), F_max_values (the MVA limit of the line), and
# Z_values (the pu values of series impedance)
def lineRead():
    infile = open(lineCSV, 'r') 

    node_from,node_to = [],[] 
    R_values, X_values, B_values,Fmax_values, Z_values = [],[],[],[],[]
    
    # Extracts the column headers from the LineData file before iterating
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


# This function reads in the bus data from the corresponding .csv file and writes out arrays with:
# Bus_num (the bus number), P_load (The active power of the load at that bus),
# Q_load (The reactive power of the load at that bus), Type (the bus type, D=load, G=generation, DG=generation and load)
# V_set (The voltage at that bus (if it has a set voltage))
def busRead():
    infile = open(busCSV, 'r')
    Bus_num, P_load, Q_load, Type, P_gen, V_set = [],[],[],[],[],[] #initialization
    # Extracts the column headers from the BusData file before iterating through
    # and organizing the data into the appropriate list
    grab_headers = False
    for line in infile:
        entries = line.split(',')
        if grab_headers:
            Bus_num.append(int(entries[0]))
            # Initialize swing bus real power at 0.
            if entries[1] == '':
                P_load.append(0)
            else:
                P_load.append(float(entries[1])) #change to pu
            # Initialize swing bus reactive power at 0.
            if entries[2] == '':
                Q_load.append(0)
            else:    
                Q_load.append(float(entries[2])) #change to pu
            Type.append(entries[3])
            # All buses without generation -> set real power generation to 0.
            if entries[3] == 'G':
                P_gen.append(0)
            elif entries[4] == '':
                P_gen.append(0)
            else:
                P_gen.append(float(entries[4]))
            # Initialize all empty V set points to 1
            if entries[5].rstrip() == '':
                V_set.append(float(1))
            else:
                V_set.append(float(entries[5].rstrip()))
        else:
            grab_headers = True
    infile.close()
    
    return(Bus_num,P_load, Q_load, Type, P_gen, V_set)

#==============================================================================
#  Functions about initializing data or things to do at the beginning
#==============================================================================

# Converts power to pu units based off the MVA base
def toPU(P_MW_gen,P_MW,Q_MVAR):
    P_gen_pu = [x / MVAbase for x in P_MW_gen]
    P_load_pu = [x / MVAbase for x in P_MW]
    Q_load_pu = [x / MVAbase for x in Q_MVAR]
    return (P_gen_pu, P_load_pu,Q_load_pu)


#==============================================================================
#  Functions about creating Y matrix
#==============================================================================

# This function computes the admittance matrix (see details within)
# Inputs: node_from, node_to, B and Z values
# Outputs: Y_matrix
def admittance_matrix(node_from, node_to, B_values, Z_values):
    
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
        line_Z = Z_values.pop(0)
        shunt = B_values.pop(0)/2
        
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
#  Functions about solving Newton-Raphson
#==============================================================================

# --------Functions about forming power system equations ----------------------

# This function solves for the computed active power
# Input: G,B,V,theta,N
# Output: Pcomp (computed P injection)
def Pcomputed(G,B, V, theta,N):
    Pcomp = [0.] * (N) #intializing
    # loop through each bus index k, then sum all P power flow equations against all buses
    for k in range(N):
        tmp = 0
        for i in range(N):
            thetaki = theta[k]-theta[i]
            tmp = tmp + V[k]*V[i]*(G[k][i]*cos(thetaki)+B[k][i]*sin(thetaki))
        Pcomp[k] = tmp
    return numpy.array(Pcomp)

# This function solves for the computed reactive power
# Input: G,B,V,theta,N
# Output: Qcomp (computed Q injection)
def Qcomputed(G,B,V,theta,N):
    Qcomp = [0.] * (N) #intializing
    #loop through each bus index k, then sum all Q power flow equations against all buses
    for k in range(N):
        tmp = 0
        for i in range(N):
            thetaki = theta[k]-theta[i]
            tmp = tmp + V[k]*V[i]*(G[k][i]*sin(thetaki)-B[k][i]*cos(thetaki))
        Qcomp[k] = tmp
    return numpy.array(Qcomp)

# This function solves the computed equations (puts above 2 together)
# Used within calc_Mismatches function and eventually finding explicit equations
# Inputs: G,B,theta,V, N
# Outputs: Pcomp (Pcomputed), Qcomp (Qcomputed)
def solveComputed(G,B,theta,V,N):
    Pcomp = Pcomputed(G, B, V, theta, N)
    Qcomp = Qcomputed(G, B, V, theta, N)
    return(Pcomp,Qcomp)

# This function calculates the mismatches and computed values
# Inputs: G,B,V0,theta0,N,Pknown,Qknown,PQbuses,PQPVbuses
# Outputs: Pcomputed, Qcomputed, and mismatch vector dPQ
def calc_Mismatches(G,B,V0,theta0,N,Pknown,Qknown,PQbuses,PQPVbuses):
    # -----solve computed values ---------
    #Pcomp = Pcomputed(G, B, V0, theta0, N)
    #Qcomp = Qcomputed(G, B, V0, theta0, N)
    [Pcomp,Qcomp]=solveComputed(G,B,theta0,V0,N)

    # ------ calculate mismatches -------
    #calculate mismatches for P and Q vectors
    dP = numpy.subtract(Pcomp, numpy.array(Pknown))
    dP = dP[PQPVbuses] #only keep PQPV bus values for P (those are known)

    dQ = numpy.subtract(Qcomp, numpy.array(Qknown))
    dQ = dQ[PQbuses] #only keep PQ bus values for Q (those are known)

    # printing max mismatch value (absolute value)
    MaxMismatch(dP, dQ, PQbuses, PQPVbuses)

    # combine mismatches to get total mismatch vector
    dPQ = numpy.concatenate((dP, dQ))
    return (Pcomp,Qcomp,dPQ)


# ----------- Functions about creating the Jacobian ---------------------

#This returns an element in the J11 section of the Jacobian matrix if j=k
def J11_same(Vk,Pk,Qk,Gkk,Bkk):
    return -Vk**2*Bkk-Qk
#This returns an element in the J11 section of the Jacobian matrix if j!=k
def J11_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*Vj*(Gkj*sin(thetakj)-Bkj*cos(thetakj))

#This returns an element in the J12 section of the Jacobian matrix if j=k
def J12_same(Vk,Pk,Qk,Gkk,Bkk):
    return Pk/Vk+Vk*Gkk
#This returns an element in the J12 section of the Jacobian matrix if j!=k
def J12_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*(Gkj*cos(thetakj)+Bkj*sin(thetakj))

#This returns an element in the J21 section of the Jacobian matrix if j=k
def J21_same(Vk,Pk,Qk,Gkk,Bkk):
    return Pk-Vk**2*Gkk
#This returns an element in the J21 section of the Jacobian matrix if j!=k
def J21_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return -Vk*Vj*(Gkj*cos(thetakj)+Bkj*sin(thetakj))

#This returns an element in the J22 section of the Jacobian matrix if j=k
def J22_same(Vk,Pk,Qk,Gkk,Bkk):
    return Qk/Vk-Vk*Bkk
#This returns an element in the J22 section of the Jacobian matrix if j!=k
def J22_diff(Vk,Vj,Gkj,Bkj,thetakj):
    return Vk*(Gkj*sin(thetakj)-Bkj*cos(thetakj))

# Function that computes each Jacobian part (aka J11, J21, J12, J22)
# Inputs: G,B,P,Q,V,theta, 2 equations for solving each element in the respective J part (one if j=k, one if j!=k)
# Also, the rangej (which buses to include in j in the jacobian part) and rangek (which buses to include in k when solving jacobian part)
# Output: one section of the jacobian matrix
def Jpartsolve(G,B,P,Q,V,theta,J_same,J_diff,rangej,rangek):
    Jpart = [[0.] * (len(rangej)) for i in range(len(rangek))] #initialization

    for j in rangej: #j is bus j
        matrix_j = rangej.index(j) #this is the position of j in the Jpart matrix you are at
        for k in rangek: #this is bus k
            matrix_k = rangek.index(k) #this is the position of k in the Jpart matrix you are at

            # If j=k, solve the J_same function
            # If j!=k, solve the J_diff function
            if j==k:
                Jpart[matrix_k][matrix_j]=J_same(V[k],P[k],Q[k],G[k][j],B[k][j])
            else:
                thetakj = theta[k] - theta[j] #difference between theta k and j
                Jpart[matrix_k][matrix_j] = J_diff(V[k], V[j], G[k][j], B[k][j],thetakj)
    return numpy.array(Jpart)

# This function computes the complete Jacobian of the power system
# Inputs: G,B,P,Q,V,theta, PQPVbuses (the index of all PQ and PV buses), PQbuses (the index of all PQ buses)
# Output: the solved Jacobian
def jacobian(G,B,P,Q,V,theta,PQPVbuses,PQbuses):
    J11 = Jpartsolve(G,B,P,Q,V,theta,J11_same,J11_diff,PQPVbuses,PQPVbuses) #solves J11 section of jacobian
    J21 = Jpartsolve(G,B,P,Q,V,theta,J21_same,J21_diff,PQPVbuses,PQbuses) #solves J21 section of jacobian
    J12 = Jpartsolve(G,B,P,Q,V,theta,J12_same,J12_diff,PQbuses,PQPVbuses) #solves J12 section of jacobian
    J22 = Jpartsolve(G,B,P,Q,V,theta,J22_same,J22_diff,PQbuses,PQbuses) #solves J22 section of jacobian
    J = numpy.concatenate((numpy.concatenate((J11,J12),axis=1),
                           numpy.concatenate((J21,J22),axis=1))) #combines them into one J matrix
    return J



# ----------------- Other misc supporting Newton Raphson functions ----------------------------------


# This function solves and prints the maximum absolute mismatch of P & Q (and which bus it's at)
# for each NR iteration
# Input: mismatch values dP and dQ, PQbuses (which buses are PQ), and PQPVbuses (which buses are PQ or PV)
# Output: nothing but prints the values.
def MaxMismatch(dP,dQ,PQbuses,PQPVbuses):
    #print largest P mismatch and the bus
    absdP = [abs(item) for item in dP] #take absolute value of each dP
    maxPmismatch = max(absdP) #finds the max absolute dP
    i =absdP.index(maxPmismatch) #backsolves for which bus this occurs at
    print("Largest P mismatch: " + str(maxPmismatch) +" pu at bus: "+str(PQPVbuses[i])) #prints the values

    #print largest Q mimatch and the bus
    absdQ = [abs(item) for item in dQ] #take absolute value of each dP
    maxQmismatch = max(absdQ) #finds the max absolute dP
    i =absdQ.index(maxQmismatch) #backsolves for which bus this occurs at
    print("Largest Q mismatch: " + str(maxQmismatch) +" pu at bus: "+str(PQbuses[i])) #prints the values
    print('')

    return

# Function that finds which are PQ and PV buses
# given a bus_types array,
# returns which index are PVbuses (PV buses), which index are PQbuses (PQ buses), and which index are PQPV buses (PQ or PV buses)
def bus_Types(bus_types):
    # PQ bus if only load
    PQbuses = [i for i in range(len(bus_types)) if bus_types[i] == "D"]

    #PV bus if only gen or gen and load
    PVbuses = [i for i in range(len(bus_types)) if (bus_types[i] == "G" or bus_types[i] == "GD" or bus_types[i] == "DG")]
    del PVbuses[0] #removing swing bus

    #PQ or PV bus only if not swing bus
    PQPVbuses = [i for i in range(len(bus_types)) if (bus_types[i] == "D" or bus_types[i] == "G" or bus_types[i] == "GD" or bus_types[i] == "DG")]
    del PQPVbuses[0] #removing swing bus

    return(PVbuses,PQbuses,PQPVbuses)

#This function determines if there is a violation in the voltage of each bus
# inputs: V (solved bus voltage vector), N (number of buses)
# outputs: Vviolation (a vector of yes or no's if the voltage is out of the given Vmin, Vmax range)
def V_violation(V,N):
    Vviolation = ["No" for i in range(N)] #initializing vector to "No" - not violating
    for i in range(N): #for each bus
        # if outside range, update violation to "Yes"
        if (V[i]>Vmax or V[i]<Vmin):
            Vviolation[i] = "Yes"
    return Vviolation


# ----------------- functions directly solving Newton-Raphson ----------------------------------


# Computes the solution of a system of equations using the Newton Raphson method
# This puts together all other functions pertaining to NR
# inputs: initial P and Q, initial theta0, initial V0, tolerance (epsilon)
# outputs: solutions for theta and V
def NewtonRaphson(bus_types,Ybus,Pknown,Qknown,theta0,V0,Eps,N):
    #extracts G and B values from Ybus (G is real element, B is imaginary element)
    G = [list(numpy.array(x).real) for x in Ybus]
    B = [list(numpy.array(x).imag) for x in Ybus]
    # which are PQ buses and PV buses?
    [PVbuses,PQbuses, PQPVbuses] = bus_Types(bus_types)
    m = len(PVbuses)+1 #m value (there are m-1 PVbuses)
    #initializations
    theta0 = numpy.array(theta0)
    V0 = numpy.array(V0)
    iteration = 0 #start at iteration 0 (this is for printing mismatch information)

    print('======== Newton Raphson Mismatch Information ========')
    print("Iteration: " + str(iteration))
    #calculate starting mismatch
    [Pcomp,Qcomp,dPQ] = calc_Mismatches(G, B, V0, theta0, N, Pknown, Qknown, PQbuses, PQPVbuses)

    # find maximum mismatch to check for convergence and iterate if greater than eps
    while abs(max(dPQ.min(),dPQ.max(), key=abs)) > Eps:
        iteration = iteration + 1
        print("Iteration: " + str(iteration))
        #solve NR for this iteration
        [theta0,V0] = NR_correction(G, B, Pcomp, Qcomp,dPQ, V0, theta0,PQbuses, PQPVbuses,N,m)
        [Pcomp,Qcomp,dPQ] = calc_Mismatches(G, B, V0, theta0, N, Pknown, Qknown, PQbuses, PQPVbuses)

    #solve explicit equations
    ## Note - because it's easier when coding, I'm just resolving implicit and explicit equations, then throwing away the previously known values
    [P_new,Q_new] = solveComputed(G, B, theta0, V0, N)
    P_new[PQPVbuses] = Pknown[PQPVbuses]
    Q_new[PQbuses] = Qknown[PQbuses]
    return (theta0,V0,P_new,Q_new)



# This function is an iteration of solving implicit equations for updated V and theta values using the Newton Raphson method
# Inputs: G, B, Pcomp, Qcomp,dPQ, V0, theta0,PQbuses, PQPVbuses,N,m
# Output: updated V0 and theta0 values
def NR_correction(G, B, Pcomp, Qcomp,dPQ, V0, theta0,PQbuses, PQPVbuses,N,m):
    #  solve for correction vector
    J = jacobian(G, B, Pcomp, Qcomp, V0, theta0, PQPVbuses, PQbuses)  # solve Jacobian
    Jinv = inv(J)  # inverse jacobian
    delta = - numpy.dot(Jinv, dPQ)  # correction vector

    #separate correction into delta V and theta
    dtheta = delta[numpy.array(range(N - 1))] #dtheta is first section of the vector (PQPV buses)
    dV = delta[numpy.array(range(N - 1, 2 * N - m - 1))] #dV is second section of the vector (only for PQbuses)

    #update the theta and V values of specified buses
    theta0[PQPVbuses] = theta0[PQPVbuses] + dtheta
    V0[PQbuses] = V0[PQbuses] + dV
    return (theta0,V0)


#==============================================================================
#  Functions about solving the power flow on lines
#==============================================================================

# This function calculate the power flow over lines and if there is a violation
# Input: V_pu (per unit voltage array), theta (V angle array), node_from array, node_to array, Fmax (maximum MVA limit array), Nl (number lines), Z_line (line impedance array), B_line (shunt suseptance)
# Output:
def calc_LineFlow(V_pu,theta, node_from,node_to,Fmax,Nl,Z_line,B_line):
    #initializing new parameters
    Pline_to2from = [0. for i in range(Nl)]
    Qline_to2from = [0. for i in range(Nl)]
    Smagline_to2from = [0. for i in range(Nl)]
    Pline_from2to = [0. for i in range(Nl)]
    Qline_from2to = [0. for i in range(Nl)]
    Smagline_from2to = [0. for i in range(Nl)]
    MVAviolation = ["No" for i in range(Nl)]
    V_complex = [0. for i in range(len(V_pu))]
    #creating a complex V
    for i in range(len(V_pu)):
        V_complex[i] = cmath.rect(V_pu[i], theta[i])
    #iterating through each line to calculate power flows
    for line in range(Nl):
        nfrom = node_from[line]
        nto = node_to[line]
        #first check the flow from -> to
        [Smag,P,Q] = calc_flows_oneline(nfrom,nto,V_complex,Z_line[line],B_line[line])
        Smagline_from2to[line] = Smag
        Pline_from2to[line] = P
        Qline_from2to[line] = Q
        #next check the flow to -> from
        [Smag,P,Q] = calc_flows_oneline(nto,nfrom,V_complex,Z_line[line],B_line[line])
        Smagline_to2from[line] = Smag
        Pline_to2from[line] = P
        Qline_to2from[line] = Q
        # This next part checks if either (from to to) or (to to from) violates the MVA limit of the line and says yes if apparent power S is greater than Fmax
        if (Smagline_from2to[line] > Fmax[line] or Smagline_to2from[line] > Fmax[line]):
            MVAviolation[line] = "Yes"
    return (Smagline_from2to,Pline_from2to,Qline_from2to,Smagline_to2from,Pline_to2from,Qline_to2from,MVAviolation)

# This function calculates the line flow for one line (is iterated through in calc_LineFlow)
# Inputs: nfrom, nto vector, complex voltage vector, Z, B of the specific line
# outputs magnitude of S over the line, P over the line, Q over the line
def calc_flows_oneline(nfrom,nto,V_complex,Z,B):
    Vdrop = V_complex[nfrom - 1] - V_complex[nto - 1] #voltage drop between the 2 bus voltages
    I_line = Vdrop / Z  # the series current
    if B == 0:  # avoiding 0/0 with if statement
        I_shunt = 0
    else:
        I_shunt = V_complex[nfrom - 1] * (complex(0, B) / 2)  # the shunt current
    I_total = I_line + I_shunt #adding series and shunt currents
    S = V_complex[nfrom - 1] * I_total.conjugate() * MVAbase  # multiply V by conjugate and rescale powers to not pu
    P = S.real  # active power
    Q = S.imag  # reactive power
    Smag = abs(S)  # apparent power
    return (Smag,P,Q)

#==============================================================================
#  Functions about wrapping up data
#==============================================================================

# rounds all values that need rounding
def roundAll(Smagline_from2to,Pline_from2to,Qline_from2to,Smagline_to2from,Pline_to2from,Qline_to2from,V,theta,P_MW_gen_solved,Q_MVAR_gen_solved):
    return (numpy.around(Smagline_from2to,nround), numpy.around(Pline_from2to,nround), numpy.around(Qline_from2to,nround), numpy.around(Smagline_to2from,nround), numpy.around(Pline_to2from,nround), numpy.around(Qline_to2from,nround), numpy.around(V,Vround), numpy.around(theta,nround), numpy.around(P_MW_gen_solved,nround), numpy.around(Q_MVAR_gen_solved,nround))


#==============================================================================
#  Functions about writing results to the output window
#==============================================================================

# This function prints the bus information to terminal
# Prints Bus number, bus type, theta, P produced by generator (MW), Q produced by generator (MVAr), V (pu), and VViolation (if voltage is violated)
def busPrint(Bus_num,bus_type,theta,P_MW_new,Q_MVAR_new,V,VViolation):
    # Printing the bus information
    print('======== Bus Information ========')    
    for i in range(max(Bus_num)):
        print('Bus ' + str(i+1))
        print('Voltage: ' + str(V[i]) + ' p.u.')
        print('Voltage violation? ' + VViolation[i])
        print('Angle: ' + str(theta[i]) + ' degrees')
        if (bus_type[i] == 'G' or bus_type[i] == 'DG'):
            print('Real Power Produced by Generator: ' + str(P_MW_new[i]) + ' MW')
            print('Reactive Power Produced by Generator: ' + str(Q_MVAR_new[i]) + ' MVAR')
        print('')
    return

# This function prints the line information to terminal
# Prints Line From to To, apparent power in that flow direction and if there is a violation in max MVA
def linePrint(node_from,node_to,Smagline,MVAviolation):
    # Printing the line information
    print('======== Line Information ========')
    for i in range(max(node_to)):
        print('')
        print('Line ' + str(node_from[i]) + ' to ' + str(node_to[i]))
        print('Apparent Power (in that direction): ' + str(Smagline[i]) + ' MVA')
        print('Apparent Power flow violation (in either direction)? ' + MVAviolation[i])
    return

#==============================================================================
#  Functions about writing output files
#==============================================================================

# This function writes the computed line data to a .csv file (lineWritefile set at beginning) and writes out columns with:
# node_from (the start bus of the line), node_to (the end bus of the line),
# P (MW) (Active power on the line in MW), Q (MVAr) (Reactive power on the line in MVAR)
# S (MVA) (Apparent power on the line in MVA), and Fmax_violation? (If the computed S_MVA value exceeds the given MVA line limit)
def lineWrite(listOfArrays):
    toPrint = zip(*listOfArrays)
    myFile = open(lineWritefile,'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerow(["", "", "From to To", "From to To", "From to To", "To to From","To to From", "To to From",""])  # write headers
        writer.writerow(["node_from", "node_to", "P (MW)", "Q (MVAr)", "S (MVA)","P (MW)", "Q (MVAr)", "S (MVA)", "Fmax_violation?"]) #write headers
        writer.writerows(toPrint)
    return

# This function writes the computed bus data to a .csv file (busWritefile set at beginning) and writes out columns with:
# BusNumber (the bus number), Type (bus type), P_gen (MW) is the generated active power, Q_gen is generated reactive power,
# V (pu) is per unit voltage at that bus, V_Violation? is if the voltage violates the reasonable limits set at the top
def busWrite(listOfArrays):
    toPrint = zip(*listOfArrays)
    myFile = open(busWritefile,'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerow(["BusNumber","Type","theta (deg)","P_gen (MW)","Q_gen (MVAr)","V (pu)","V_Violation?"]) #write headers
        writer.writerows(toPrint)
    return



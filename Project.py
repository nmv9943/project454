import numpy
import openpyxl
import cmath
import math
import os

########################
# Input File
########################


scenarioName = "1"

#basefolder = "/Users/ninavincent/Desktop/Git/project454/"
basefolder = os.getcwd()
database = basefolder + "/data/"
lineCSV = database + "LineData" + scenarioName + ".csv"
busCSV = database + "BusData" + scenarioName + ".csv"


########################
# Constants
########################

Eps = 0.001
MVAbase = 100
#something else = 20 but can't read my notes


#N = 2
#Matrix = [[0.] * N for i in range(N)]


###################################
# Functions about read/write files
##################################






###################################
# Functions about creating Y matrix
##################################




###################################
# Functions about forming power system equations
##################################




###################################
# Functions about Newton-Raphson
##################################



# Function that computes Jacobian


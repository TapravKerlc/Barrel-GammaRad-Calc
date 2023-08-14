import csv
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
import pylab
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import math
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib import pyplot
import matplotlib.pyplot as plt
from pylab import *
import sys
import time
from scipy.integrate import nquad
import random
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


################## Variables ##############################################

R = 0.1 #radij
H = 1.7 #height

E1 = 2 #začetna Energija v MeV
currentE = E1
currentKot = 0

tries = 5000#število delcev

dh = 0.01 #velikost koraka

odloženaE = 0 #variabla: kok energije so pustil delci notr
counterPobegli = 0

################## Generacija delcev v geometriji #########################
m = np.random.random(500)
n = np.random.random(500)
k = np.random.random(500)
l = np.random.random(500)

phi = n*2*np.pi
theta = m*2*np.pi
r = R*np.sqrt(k)
h = H*l

x = r*np.cos(phi)
y = r*np.sin(phi)
z = h

#################### Prvi graf ############################################

#fig = plt.figure(figsize=(8,8))

#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(x,y,z) # plot the point (2,3,4) on the figure


#plt.show()

#################### Drugo ################################################

# csv file name
filename = "voda.csv"
 
# initializing the titles and rows list
fields = []
rows = [] #0,5 MeV je 21ti
 
# reading csv file
with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
     
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        rows.append(row)
 
    # get total number of rows
    #print("Total no. of rows: %d"%(csvreader.line_num))

numLines=csvreader.line_num

# printing the field names
#print('Field names are:' + ', '.join(field for field in fields))
 
#  printing first 5 rows
#print('\nFirst 5 rows are:\n')
fixdrows0 = []
fixdrows1 = []
fixdrows2 = []
fixdrows3 = []
currentErow = 0 #na kateri energiji se trenutno nahajamo
for row in rows[:]:
    # parsing each column of a row
    #print(row)
    stringtorow = str(row)
    stringtorow0 = float(stringtorow[2:11]) #float verzija energije
    fixdrows0.append(stringtorow0)
    stringtorow1 = float(stringtorow[13:21]) #float verzija CS
    fixdrows1.append(stringtorow1)
    stringtorow2 = float(stringtorow[24:32]) #float verzija Fotoefekt
    fixdrows2.append(stringtorow2)
    stringtorow3 = float(stringtorow[35:43]) #float verzija tvorba parov
    fixdrows3.append(stringtorow3)
    #print(stringtorow3)
    
for row in range(80):
        if currentE>=fixdrows0[row] and currentE<fixdrows0[row+1]:
            currentErow = row + 1

currentCSpath = fixdrows1[currentErow]*10
#print(currentCSpath)
currentFEpath = fixdrows2[currentErow]*10
currentTPpath = fixdrows3[currentErow]*10

def compton(entryE, začetnKot):
    alfa = entryE/1.022 #1.022= 2mc^2
    epsilon0 = 1/(1+2*alfa)
    alfa1 = (1-epsilon0**2)/2
    alfa2 = -math.log(epsilon0)    
    bully = True
    
    ####Sledi ful zakomplicirana funkcija za določanje epsilona, skupej s preverjanjem, da je prava rešitev (bully = True pomeni napačna rešitev)
    
    while bully:
        xi1 = np.random.random(1)
        xi2 = np.random.random(1)
        if xi1 < alfa1/(alfa1+alfa2):
            epsilonČ = math.sqrt(epsilon0**2+2*(alfa1+alfa2)*xi1)
        else:
            epsilonČ = epsilon0*math.exp((alfa1+alfa2)*xi1-alfa1)
            
        if xi2 < (1-(epsilonČ*(1-epsilonČ**2))/((1+epsilonČ**2)*(alfa*epsilonČ)**2)): #preverjamo če rešitev ustreza
            bully = False
    
    novKot = math.acos(1-((1-epsilonČ)/(alfa*epsilonČ))) + začetnKot
    exitE = epsilonČ*entryE
    global odloženaE
    odloženaE = odloženaE+(1-epsilonČ)*entryE
    #print(exitE, novKot)
    return exitE, novKot

def potdoroba (phi, theta, r, h):
    yProjekcijaPoti = (R-r)*np.sin(phi) #višina od spawna delca do izhoda
    #print(h, yProjekcijaPoti)

    if H-(yProjekcijaPoti+h) >= 0 and h+yProjekcijaPoti > 0: #če gre skozi plašč soda
        dRoba = np.sqrt(yProjekcijaPoti**2+(R-r)**2)*np.sin(0.5*theta)*(R-r)/R #ubistvu pitagora
        #print("plašč")
        
    elif (h+yProjekcijaPoti <= 0): #če gre skozi dno soda        
        dRoba = np.sqrt(h**2+(h*np.sin(phi))**2)*np.sin(0.5*theta)*(R-r)/R #koren(višina * sinus (višina))
        #print("                             dno")
        
    else:   #če gre skozi pokrov soda    
        dRoba = np.sqrt((H-h)**2+ ((H-h)*np.cos(90-phi))**2)*np.sin(0.5*theta)*(R-r)/R #koren(kateta trikotnika do vrha soda (y) + kateta x na vrhu soda) = hipotenuza oz pot delca
        #print("         vrh")
        
    #print(dRoba, comptonPopravek)
    #dRoba = dRoba - comptonPopravek
    return dRoba
        


for x in range(tries):
    
    ###########################particular particle stats#################################
    currentE = E1
    currentKot = 0    
    comptonPopravek = 0
    
    m = np.random.random(1)
    n = np.random.random(1)
    k = np.random.random(1)
    l = np.random.random(1)
    
    phi = n*2*np.pi
    theta = m*2*np.pi
    r = R*np.sqrt(k)
    h = H*l    
    
    while currentE > 0:
        dRoba = potdoroba(phi, theta, r, h)
        for row in range(80):
            if currentE>=fixdrows0[row] and currentE<fixdrows0[row+1]:
                currentErow = row + 1
        
        currentCSpath = fixdrows1[currentErow]*10
        currentFEpath = fixdrows2[currentErow]*10
        #currentTPpath = fixdrows3[currentErow]
        
        dPE = -math.log(np.random.random(1))/currentFEpath
        dCS = -math.log(np.random.random(1))/currentCSpath
        
        #vstavi pot do roba tle dRob = potdoroba (xlega, ylega, zlega, usmerjenostKot1, usmerjenostKot2)
        #dRoba = 0 #POBRIŠI ME
        #print(dRoba, dCS, dPE)
        if (dRoba <= dCS) and (dRoba <= dPE): #escape
           counterPobegli += 1
           currentE = 0
        elif (dCS <= dPE): #compton
           currentE, theta = compton(currentE, theta)
           comptonPopravek = dCS
        else: #photoeffect
           odloženaE = odloženaE + currentE
           currentE = 0
        
        #print("Odložena E: ", odloženaE)


    #for col in row:
     #   print("%10s"%col,end=" "),
#print(fixdrows0)        
        
   # print('\n')
    #print(rows[1])
print(counterPobegli, "Jih pobegne")
print(odloženaE, " MeV")
print(odloženaE/tries, "energije na delec")
print("Razmerje je", odloženaE/(tries*E1), "%")
###DOZA

doza = odloženaE/(np.pi * R**2*H)
    
print(doza, " Gy")

    
    #vektor kam gre, vektor do sredine, dot product. Za pokrov ali ne: če je višina valja-višina izhoda < 0, poj gre skoz pokrov. Obrneš situacijo in je še laži.
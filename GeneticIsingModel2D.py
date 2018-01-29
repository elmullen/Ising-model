# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime

# this algorithm investigates genetic fitness over sucsesive generations- this simplified model of reality requires that populations size(x) is held constant 

# two parts to this investigation: 
# part 1 = getic algorithm fitness vrs generation - implement for fixed temprature- temprature effects how dominat genes are 
# part 2= how this implemenation improves rate of convergence to a state of equilibrium for the isisng model 

#defining lattice structure 
#----------------------------------------------------------------------------------------------------------------------------------------
#creates 2D lattice 
def make_lattice(rows,columns): 
    lattice = np.zeros((rows, columns))
    for i in range(rows):
        for j in range(columns):
                lattice[i][j] = np.random.choice((-1,1)) # spin glases- spin states 180 deg seperation denoted by -1 and 1 
    return lattice





# conversions between lattice and string 
#--------------------------------------------------------------------------------------------------------------------------------
#re-shape returened lattice- to single line- turn into array
def convert_string(lattice):
    lattice = np.asarray(lattice) # turns into an array 
    string = np.reshape(lattice, ( 1, -1))[0] # reshapes array into rows 
    for i in range(0,len(string)):
        if(string[i] == -1.):
            string[i] = 0. # conversion between spin state and binary
    return string

#re-shape returned array- to lattice 
def convert_lattice(string):
    string = np.asarray(string) 
    sqrt = np.sqrt(string.size)# square root gives dimensons of lattice
    sqrt= int(sqrt) # define as integer value 
    lattice = np.reshape(string, (sqrt, sqrt)) # reshape to give lattice of defined integer value 
    for i in range(0,lattice.shape[0]):
        for j in range(0,lattice.shape[1]):
            if(lattice[i][j] == 0):
                lattice[i][j] = -1 # conversion between spin state and binary
    return lattice  

# note on notation: let x be the population number, x/2 is the number of offspring resulting from the genetic algorithm, x/2 + x is the number of parents and offspring



# mutation and recombination- genetic stuff 
#------------------------------------------------------------------------------------------------------------------------------------

# code for recombination 
def recombination(strings):
    #print strings
    Z = len(strings) # counts the number of strings  
    N = strings[0].shape[0]
    n= Z/2 # divdes the number of strings in half 
    for i in range(0,n,2):# selects strings in groups of 2  
        sindex = np.random.randint(0,N) # select a random integer(sindex) on each of the two selected strings 
        lenght = np.random.randint(sindex,N) # length with starting point at sindex is same magnitude for both strings
        #print "starting index and length : ", sindex, lenght
# define two strings 
        newString1 = strings[i]
        newString2 = strings[i+1]
       # print newString1, newString2 - debug code 
        ex_data = 0 # holder variable- holds the value of one string as it is changed 
# change some value of string 1 with a value of new string 2- ex data holds the string 
        for j in range(sindex,lenght):
            ex_data = newString1[j] 
            newString1[j] = newString2[j] # updates selected lenght from string 1 with selected lenght from string 2
            newString2[j] = ex_data # updates selected lenght from string 2 with selected lenght from string 1
        #print "end"
        #print newString1, newString2
        strings.append(mutation(newString1)) # mutate and update string 1 
        strings.append(mutation(newString2)) # mutate and update string 2 
        #print i         
    #print strings
    return strings 

# 1/L**2 chance of mutation at each string site for offspring

def mutation(string):
    prob = 1 - 1./(string.shape[0]**2)
    #print prob
    for i in range(0,string.shape[0]): 
        if prob < np.random.random():
            if string[i] == 1:
                string[i] = 0
            else:
                string[i] = 1
    return string
# change to energy 
# interpretation of this part of algorithm- some genes are more favorable than others- spin states which contribute to energy reduction represent the effect of dominant alles- thus the issing model can be seen as a conversion from genotype to phenotype. 
def energyfavorability(lattice,t):
# finding energy value/ before flip acording to issing 
    for i in range(0,lattice.shape[0]):
       for j in range(0,lattice.shape[1]): 
           right = lattice[(i+1) % lattice.shape[0]][j]
           left = lattice[(i-1) % lattice.shape[0]][j]
           bottom = lattice[i][(j+1) % lattice.shape[1]]
           top = lattice[i][(j-1) % lattice.shape[1]]
           energy = -1 * lattice[i][j] * (left + right + top + bottom)
           flip_energy = -energy
           energy_difference = flip_energy - energy
           if np.exp(-(energy_difference)/(t)) > np.random.random():
               lattice[i][j] *= -1 # energy 
    return lattice

# population now consists of: x/2 + x (parents and offspring) 

#energy 
#----------------------------------------------------------------------------------------------------------------
# recombination and mutation for the current generation 
def genetic_elements(lattice,t):
    strings = []
    n = len(lattice)
    #print "small n", n
    for i in range(0,n):
        strings.append(convert_string(lattice[i])) # convert to string from lattice 
    #return
    #print "strings", len(strings)
    rec = recombination(strings) # genetic recombination and offspring mutation
    #print "rec length", len(rec)
    lattice = []
    for i in range(0,int(n+(n*0.5))):
        #print i, n
        lattice.append(energyfavorability(convert_lattice(rec[i]),t)) #turn into lattice again and apply energy favourability function
        #print lattice[i]
    return lattice


# survival of the fitest 
#------------------------------------------------------------------------------------------------------------------------------------
# determining the fitness of each lattice

def fitness_per_lattice(L,S):
    fitness = 2*L**2 - S # formula for calculating genetic fitness 
    return fitness 
    
# determining which which lattices are kept and go progress to next generation 
def fitness_keep(fitnesses,lattices):     
    fitsum = []    
    value = 0 
    results = []
    for a in range(0, len(fitnesses)):
        for z in range(0,a):
            value += fitnesses[z] # add each value of fitness to the previous one 
        fitsum.append(value) # record these values 
    #print fitsum
   # print int(fitsum[len(fitsum)-1])
    #return
    indexlist = []
    while len(results) < len(fitnesses)*(2./3.): # selects the same number of survivors as origonal population size(x) from X+X/2
        #R = fitsum[int(len(fitsum))-1]
        r = np.random.randint(0,int(fitsum[len(fitsum)-1])) # select a random number 
        #print r
       # print R
        for x in range(0,len(fitsum)):
            if(r > fitsum[x] and x not in indexlist): # if this value lies wihin the fitsum keep the last value of fitness last added to the fitness array that makes this statement true
                results.append(lattices[x])
                indexlist.append(x)
                break
    #print i
    #print results
    return results

#------------------------------------------------------------------------------------------------------------------------------------------------

def calc_energies(lattice): # function used to calculate lattice energies, both the average of each site(which is used for the varience calculation later) 
    totalenergy = 0
    increment = 0
    magnetisation = 0
    mag_list = []
    energylist = []
   # print lattice.shape[0],lattice.shape[1]
    for i in range(0,lattice.shape[0]):
        for j in range(0,lattice.shape[1]):
            right = lattice[(i+1) % lattice.shape[0]][j]
            left = lattice[(i-1) % lattice.shape[0]][j]
            bottom = lattice[i][(j+1) % lattice.shape[1]]
            top = lattice[i][(j-1) % lattice.shape[1]]
            energy = -1 * lattice[i][j] * (left + right + top + bottom)
            magnetisation += lattice[i][j]
            mag_list.append(lattice[i][j])
            totalenergy += energy
            energylist.append(energy)
            increment += 1
    #print totalenergy
    return totalenergy/increment, energylist, magnetisation/increment, mag_list
    
def implementation():
    #initialises all of the variables
    starttime = datetime.datetime.now()
    lattice_size = 25#lattice dimension
    lattices = []
    for i in range(0,4): # range selected depends on origonal population size- the greater the population size the more realisitic the algorithm 
        lattices.append(make_lattice(lattice_size,lattice_size))
    #holding arrays for data
   # fitness_keep = [] 
    tlist = []
    energylist = []
    energyvar = []
    magnetisation = []
    magnetic_variance = []
    sfitnesses = []
    fitnesskeep = []
    t = 0.01 
    tmax = 5.0 # if t= 0.01 tmax= 0.02 for const temprature to get fitness vrs energy 
    #print lattices
    while t < tmax:
        z = 0 
# temprature loop 
        while z <= 5: # number of generations code runs for 
# calculating the fitness of each lattice for enlarged population: x + x/2 (parent and offspring) 
            fitnesskeep = [] # number of survivors must equal to origonal population number(x)- x is held constant 
            sfitnesses = [] # x +x/2  
            lattices = genetic_elements(lattices,t)
            for i in range(0,len(lattices)):
                #print fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i]))
                sfitnesses.append(fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i])[0]))  
            
            lattices = fitness_keep(sfitnesses,lattices)
            for i in range(0,len(lattices)):
                #print fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i]))
                fitnesskeep.append(fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i])[0]))
            #fitness = sum(fitnesskeep)
            #print fitness
            #print lattices
            z += 1 # loop adds one to generation number every time it is run 
     #      print z
        print "temperature", t
        dataholder = calc_energies(lattices[0])
        energylist.append(dataholder[0])
        energyvar.append(np.var(dataholder[1]))
        magnetisation.append(dataholder[2])
        magnetic_variance.append(np.var(dataholder[3]))
        tlist.append(t)
        t+=0.01
        
    for i in range(0,len(magnetisation)):
        if(magnetisation[i] <= 0):
            print "Critical temperature : ", tlist[i]
            break
    specific_Heat_capacity = []    
    #print len(tlist), len(energyvar)
    for i in range(len(energyvar)):
        specific_Heat_capacity.append(energyvar[i] / (tlist[i] ** 2))
        
    magneticSuscep = []

    for i in range(len(magnetic_variance)):
        magneticSuscep.append(magnetic_variance[i] / tlist[i])
        
    endtime = datetime.datetime.now()
    print endtime - starttime
    #print tlist, energylist
    
    plt.figure(1)
    plt.xlabel("Temperature(J/Kb)")
    plt.ylabel("Average lattice Energy(J)")
    plt.title("Average energy vs Temperature(" + str(lattice_size) + "x" + str(lattice_size) + " lattice)")
    plt.plot(tlist,energylist , 'b--')
    
    plt.figure(2)
    plt.xlabel("Temperature(J/Kb)")
    plt.ylabel("Specific Heat Capacity(Cv)(J/K*2)")
    plt.title("Specific Heat Capacity vs Temperature (" + str(lattice_size) + "x" + str(lattice_size) + " lattice)")
    plt.plot(tlist,specific_Heat_capacity , 'b--')
    
    plt.figure(3)
    plt.xlabel("Temperature(J/Kb)")
    plt.ylabel("Average Magnestism(|<m>| per site [])")
    plt.title("Average Lattice Magnetism vs Temperature (" + str(lattice_size) + "x" + str(lattice_size) + " lattice)")
    plt.plot(tlist,magnetisation , 'b--')
    plt.show()
           
    plt.figure(4)
    plt.xlabel("Temperature(J/Kb)")
    plt.ylabel("Magnetic Susceptability (mu/Kb)")
    plt.title("Magnetic Susceptability vs Temperature (" + str(lattice_size) + "x" + str(lattice_size) + " lattice)")
    plt.plot(tlist,magneticSuscep , 'b--')
    plt.show()
            
def fitnesstest():
    z = 1
    lattice_size = 100#lattice dimension
    lattices = []
    for i in range(0,8): # range selected depends on origonal population size- the greater the population size the more realisitic the algorithm 
        lattices.append(make_lattice(lattice_size,lattice_size))
    fitnesslist = []
    generations = []
# temprature loop 
    while z <= 10: # number of generations code runs for 
# calculating the fitness of each lattice for enlarged population: x + x/2 (parent and offspring) 
        fitnesskeep = [] # number of survivors must equal to origonal population number(x)- x is held constant 
        sfitnesses = [] # x +x/2  
        lattices = genetic_elements(lattices,1)
        for i in range(0,len(lattices)):
            #print fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i]))
            sfitnesses.append(fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i])[0]))  
        
        lattices = fitness_keep(sfitnesses,lattices)
        for i in range(0,len(lattices)):
            #print fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i]))
            fitnesskeep.append(fitness_per_lattice(lattices[i].shape[0],calc_energies(lattices[i])[0]))
        fitness = sum(fitnesskeep)
        #print fitness
        #print lattices
        generations.append(z)
        fitnesslist.append(fitness)
        z += 1 # loop adds one to generation number every time it is run 
     #      print 
    plt.figure(1)
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.title("Fitness vs Generation (" + str(lattice_size) + "x" + str(lattice_size) + " lattice)")
    plt.plot(generations,fitnesslist , 'b--')
    plt.show()
implementation() 
#fitnesstest()
#print np.random.random()

import numpy as np
import random
import matplotlib.pyplot as plt
from mpi4py import MPI #importing multiple processor module.
import datetime

#Prerequisites for execution.
# 1) Working MPI implementation installed.
# 2) MPICC compiler wrapper
# 3) MPI4PY library

#To install prerequisites run the following code
#sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev
#sudo apt-get install python-pip python-dev build-essential
#sudo pip install mpi4py

#To execute the code run the following command 
# mpiexec -n 4 python 1DIsing.py
# In the -n [number] case, number stands for number of processor cores to run the code over.
# If a different number than 4 is available, modify this argument.

#variables needed for multiprocessor code.
#MPI= message passing interface
comm = MPI.COMM_WORLD # comunicator 
size = float(comm.Get_size()) # returns the number of processes in the communicator 
rank = float(comm.Get_rank()) # processes given rank in comumnicator 

# generate lattice 
def make_lattice(rows): 
    lattice = np.zeros((rows))
    for i in range(rows):
        lattice[i] = 1 # lattice is filled with +1 spins initially. 
    return lattice

def H_energy(lattice,lattice_size_x, x):
#Calculates energy of a lattice site, usage of modulo operator instead of cascading if statements
#saves computation time.
    right = lattice[(x+1) % lattice_size_x]
    left = lattice[(x-1) % lattice_size_x]
#here the hamiltonian is computed as: 
    return -1 * lattice[x] * (left + right)
 
def spin_energy(lattice, lattice_size_x, i, T): 
    energy = H_energy(lattice, lattice_size_x, i) # get energy of current site
    flip_energy = -energy #flip it
    energy_difference = flip_energy - energy #get energy difference
    if energy_difference < 0: #if energy is favourable flip spin
        lattice[i] *= -1 #flips spin
    elif np.exp(-(energy_difference)/(T)) > random.random() : #if probability passes flip spin
        lattice[i] *= -1 #flips spin
    energy = flip_energy 
    return lattice, energy #returns updated lattice and the energy of the site

# converging lists- multiprocessing 
def MPI_converge(clist):
    holderlist = []
    for a in range(0,len(clist)):
            holderlist += clist[a]
    return holderlist
    
def Metropolis_MonteCarlo():
    #initialises all of the variables
    starttime = datetime.datetime.now()
    #setting up lattice dimensions
    lattice_size_x =5
    #lattice generation
    lattice = make_lattice(lattice_size_x)
    #loop variables    
    samples = 25 #multiplier for the number of times the data is sampled from an equiliblirated lattice 
    Temp_inital= 0.01 # initial temperature
    T_step = 0.01 # temperature step
    Temp_max = 10  # base max temperature
    
    #holding arrays for data
    Tlist = [] 
    lattice_energies = []
    lattice_energy_variance = []
    lattice_magnetism = []
    lattice_magnetism_variance = []

       
    #Setting up Temperature loop variables over MPI
    T = Temp_inital + (rank * (Temp_max / size)) 
    Temp_max = (rank+1)*(Temp_max /size)
   # print "Temp = ", T
    #print "Temp_max = ", Temp_max #DEBUG STATEMENTS
    while T< Temp_max:
        i = 0# iteration variable
        #equilibrium loop is initialised- implementing monte carlo random sampling to equilibriate

# number of points selected for equiliberation- greater than number of lattice points avalible- ensure every lattice point multiple times- randomisation may not include all points otherwise      
        while i< (lattice_size_x) * 50: 
# random sampling according to monte carlo method to equiliberate the lattice 
            x = random.randint(0,lattice_size_x-1) 
            lattice = spin_energy(lattice, lattice_size_x, x, T)[0]
            i += 1
# equlibrium achieved 
        i = 0
        energy = 0
        magnetism = 0
        energy_array = []
        magnetism_array = []
        #samples values by picking many sites and averaging over the lattice
        while i < (lattice_size_x) * samples:  
            x = random.randint(0,lattice_size_x-1)
#collect data in arrays 
            energy_holder = H_energy(lattice,lattice_size_x, x) # energy value for a lattice site 

            magnetism_array.append(lattice[x])# magnetic value for a lattice site 
            energy_array.append(energy_holder) # energy value for a lattice site 
            energy += energy_holder #  add energy to total energy for averaging 
            magnetism += lattice[x]# add magnetism to total magnetism for averaging 
            i += 1
        energy = energy / ((lattice_size_x) * samples) # average of energy over the crystal
        magnetism = magnetism / (((lattice_size_x) * samples)) # average of magnetisation over the crystal
        lattice_energies.append(energy)
        lattice_energy_variance.append(np.var(energy_array))# np.var returns the variance of an array of elements- in this case energy varience 
        lattice_magnetism.append(magnetism)
        lattice_magnetism_variance.append(np.var(magnetism_array)) # magnetic variance 
        Tlist.append(T)
        T+= T_step
        
    
    specific_Heat_capacity = []    
    
    for i in range(len(lattice_energy_variance)):
        specific_Heat_capacity.append(lattice_energy_variance[i] / (Tlist[i] ** 2))# statistical thermodynamics defintion of specific heat capacity

    magnetic_susceptibility = []

    for i in range(len(lattice_magnetism_variance)):
        magnetic_susceptibility.append(lattice_magnetism_variance[i] / Tlist[i]) # statistical thermodynamics defintion of magnetic susceptability 
   
    #MPI variable/list gathering.
    MPI_Energy = []
    MPI_Tlist = []
    MPI_Mag = []
    MPI_magnetic_susceptibility = []
    MPI_specific_Heat_capacity = []
    
    MPI_Energy = comm.gather(lattice_energies, root = 0) #communicator holds MPI functions, in this case the gather function 
    MPI_Tlist = comm.gather(Tlist,root = 0)
    MPI_specific_Heat_capacity = comm.gather( specific_Heat_capacity, root = 0)
    MPI_Mag = comm.gather(lattice_magnetism, root = 0)
    MPI_magnetic_susceptibility = comm.gather(magnetic_susceptibility, root = 0)
    if rank == 0: #this only runs on the first core.
        #print "head processor reached"
#converge the data
        MPI_Tlist = MPI_converge(MPI_Tlist)
        MPI_Energy = MPI_converge(MPI_Energy)
        MPI_specific_Heat_capacity = MPI_converge(MPI_specific_Heat_capacity)
        MPI_Mag = MPI_converge(MPI_Mag)
        MPI_magnetic_susceptibility = MPI_converge(MPI_magnetic_susceptibility)
        for a in range(0,len(MPI_Mag)):
            if(MPI_Mag[a] <= 0): #gets the critical temperature from magnetism curve
                print "Critical Temperature = ", MPI_Tlist[a]
                break

        endtime = datetime.datetime.now()
        print "efficiency evaluation(running time(Hrs:min:sec)) = ", endtime - starttime #testing for efficiency

        # results to plot 

        plt.figure(1)
        plt.xlabel("Temperature(J/KbT)")
        plt.ylabel("Average Magnestism(|<m>| per site [])")
        plt.title("Average Lattice Magnetism vs Temperature("+ str(lattice_size_x)  +" lattice)")
        plt.plot(MPI_Tlist, MPI_Mag, 'b--')

        plt.figure(2)
        plt.xlabel("Temperature(J/KbT)")
        plt.ylabel("Average lattice Energy(J)")
        plt.title("Average energy vs Temperature(" + str(lattice_size_x) +  " lattice)")
        plt.plot(MPI_Tlist, MPI_Energy, 'b--')
        
        plt.figure(3)
        plt.plot(MPI_Tlist, MPI_magnetic_susceptibility, 'b--')         
        plt.xlabel("Temperature(J/KbT)")
        plt.title("Magnetic Susceptability vs Temperature(J)("+ str(lattice_size_x) + " lattice)")
        plt.ylabel("Magnetic Susceptability (mu/K)")

        plt.figure(4)
        plt.plot(MPI_Energy, MPI_magnetic_susceptibility, 'b--')
        plt.xlabel("Average lattice Energy(J)")
        plt.title("Average lattice Energy vs Magnetic Susceptability)("+ str(lattice_size_x) + " lattice)")
        plt.ylabel("Magnetic Susceptability (mu/K)")
 	
        plt.figure(5)
        plt.plot(MPI_Tlist,MPI_specific_Heat_capacity, 'b--')
        plt.xlabel("Temperature(J/KbT)")
        plt.title("Specific Heat vs Temperature("+ str(lattice_size_x) +  " lattice)")
        plt.ylabel("Specific Heat Capacity(Cv)(J/K*2)")
        plt.show()
            
Metropolis_MonteCarlo() 
    

    





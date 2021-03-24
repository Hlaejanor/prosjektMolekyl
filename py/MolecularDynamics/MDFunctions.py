# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:33:06 2021

@author: jenst
"""
import numpy as np
import matplotlib.pyplot as plt
import math as math
from MolecularDynamics import constants
from MolecularDynamics import MDGenerators

def MoveToCentreOfBox(particles):
    c = [constants.box /2,constants.box /2,constants.box /2]
    
    for part in particles:
        part.x[0] += c

def LennardJonesPotential(dist, eps, sig, shifted = 0):

        return  4*eps*(
                (sig / dist)**12 - (sig/dist)**6 
                )- shifted
    
    
def CreateRandomVelocityFromTemperature(particleCount, temperature):
    if temperature is None:
         raise Exception("Temperature was none")
    t0 = constants.epsEV / constants.kBoltz
    
    c = np.random.normal(0, np.sqrt(temperature/t0), size=(particleCount,3))
         
    return  c
     


def CreateRandomVelocityFromAvgSpeed(n, speedSigma, N):
     if speedSigma is None:
         raise Exception("Average speed in sigma was none")
     v0 =  np.random.normal(0, math.sqrt(speedSigma), size=(n,3))
     return v0


def ComputeTemperature(t, particles):
    print(f"Will compute temperature from looking at velocity of atoms at time {t}, mass {massMetric}")
    totalKineticEnergySquared  = 0
    for p in particles:
        vel = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        totalKineticEnergySquared += vel**2
    
    N = len(particles) 
    T = massMetric / (3 * bltz * N)  * totalKineticEnergySquared
    print(f"Temperature at time {t} is {T} kelvin")



def ComputeTemperature(t, particles):
    print(f"Will compute temperature from looking at velocity of atoms at time {t}, mass {massMetric}")
    totalKineticEnergySquared  = 0
    for t in range(0, N):
        for p in particles:
            vel = p.v[t][0]**2 + p.v[t][1]**2 + p.v[t][2]**2
            totalKineticEnergySquared += vel
    
    N = len(particles) 
    T = massMetric / (3 * bltz * N)  * totalKineticEnergySquared
    print(f"Temperature at time {t} is {T} kelvin")

def ArgonTimeToSeconds(t):
    print(f"Converting from reduced time untis to seconds for Argon gas")
    tau = np.sqrt(constants.massMetric / constants.epsEV)
    
    return t * tau, "seconds"

def SigmaToAngstrom(series):
    print(f"Converting timeseries in sigma units to Ångstrom")
    
    return series * constants.sigmaAnsgtrom, "Å"
    


def ArgonTPrimeToKelvin(tPrimeSeries):
    print(f"Converting from reduced units to kelving for Argon gas")
    t0 = constants.epsEV / constants.kBoltz
    print(f"Should be 19.7 : {t0}")
    if np.abs(t0 - 119.7 ) > 0.1: 
        raise Exception(f"Argon T0 should be 119.7, was {t0}")
    
    return tPrimeSeries * t0
    



def ComputeMeanSquaredDisplacement(t,  particles):
    print(f"Will compute mean squared displacement, returns time-series of squared vectors")
    N = len(t)
    msd = np.array([[0., 0., 0.] for i in range(N)])# Velocity
    
    for t in range(1, N):
        #print(f"Time {t}")
        for p in particles:
            msd[t] += p.d[t][0]**2 + p.d[t][1]**2 + p.d[t][2]**2              #Square the displacement from initial position, and thats it
    
    
    pCount = len(particles)
    msd = msd / (pCount) 
    diffusion_constant = (msd / (6*t))
    return  msd, diffusion_constant

def ComputeTemperatureReduced(N, particles):
    print(f"Will compute reduced temperature")
    tPrime  = np.zeros(N)
    
    for t in range(0, N):
        #print(f"Time {t}")
        for p in particles:
            tPrime[t] += p.v[t][0]**2 + p.v[t][1]**2 + p.v[t][2]**2 
    
    
    pCount = len(particles)
    
    tPrime = tPrime / (3 *pCount)
    
    return tPrime


def ComputeVelocityAutoCorrelation(t, particles):
    print(f"Will compute velocity autocorrelation")
   
    vac = np.zeros(len(t)) # Will hold a scalar for each timestep
    diffu_coeff = np.zeros(len(t))
    N = len(particles)          #The number of particles to sum over
    size = len(t)               #To get the index 
    
    for p in particles:
        scale = np.dot(p.v[0] , p.v[0])
        for i in range(0, size):
            vac[i] += np.dot(p.v[i], p.v[0]) / scale
            
    vac = vac / N
    
    dt = t[1] - t[0] # Assuming linearly spaced time-series
    
    diffu_coeff = (vac * dt) / 3
    
    return vac, diffu_coeff



## Make a random position with mean at given distance (in sigma units)
def getPositionInBox(count):
    
    slots = int(count**(1/3)) +2
    dist = constants.box / slots
    
    
    pos_x = np.linspace(0, constants.box, slots)
    pos_y = np.linspace(0, constants.box, slots)
    pos_z = np.linspace(0, constants.box, slots)
    
    
    x = np.array([[0., 0., 0.] for i in range(count)]) # Position
    index = 0
    for i in range(0, slots):
        for j in range(0, slots):
            for k in range(0, slots):
                
                if index >= count :
                    break ;
                   
                x[index] = [pos_x[i], pos_y[j], pos_z[k]]
                index += 1
    
    return x

## Make a random direction vector (normalized)
def getRandomDirectionWithTemperature(T):

    direction = np.array(np.random.normal(0, math.sqrt(T), size=(N,3)), np.random.normal(0, math.sqrt(T), size=(N,3)), np.random.normal(0, math.sqrt(T), size=(N,3)))
    
    return direction

def getRandomNormalizedDirection():
    direction= np.array([np.random.uniform()- 0.5, np.random.uniform() - 0.5, np.random.uniform() - 0.5])
    return direction / np.linalg.norm(direction)
    
## Make a random velocity (scalar), combine with getRandomDirection() to get velocity vector
def getRandomSpeed():
    
    return np.random.uniform()*0.8


## Get length of vector
def getLength(vec):
    return np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

## Get distance between two particles       
def getDistance(p1, p2):
    
    distvec = np.subtract(p1, p2)
    return np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)


def distance(a, b):
    dx = abs(a[0] - b[0])
    x = min(dx, abs(constants.box - dx))
     
    dy = abs(a[1] - b[1])
    y = min(dy, abs(constants.box - dy))
     
    dz = abs(a[2] - b[2])
    z = min(dz, abs(constants.box - dz))
 
    return np.sqrt(x**2 + y**2 + z**2)

def rdf(bin_edges, r, V):
    """
    bin_edges = edges of bins. Typically np.linspace(0, rc, num_bins+1)
    for some cut-off rc.
    r = Nx3-array of positions of atoms at a given timestep.
    V = volume of system.
    """
    N = r.shape[0]
    bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_sizes = bin_edges[1:] - bin_edges[:-1]
    n = np.zeros_like(bin_sizes)
    for i in range(N):
        if i % 1000 == 0:
            print(f"Binning {i} of {N} distances, {100*(i/N):11.3f}% commplete")
        dr = np.linalg.norm(r - r[i], axis=1) # Distances from atom i.
        n += np.histogram(dr, bins=bin_edges)[0] # Count atoms within each
        # distance interval.
    n[0] = 0
    
    # Equation (7) on the preceding page:
    rdf = V / N**2 * n / (4 * np.pi * bin_centres**2 * bin_sizes)
    return rdf


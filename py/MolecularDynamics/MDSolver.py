# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:34:12 2021

@author: jenst
"""
        
import numpy as np 
import multiprocessing as mp
import math as math
from MolecularDynamics import MDFunctions, constants
## Create particless


## Velocity Verlet
def velocityVerlet(p, n,  dt, particles, eps, sig):
    print(f"VVslow : frame {n}")
    scaleBy = 24
    acc_vec = np.array([0.0, 0.0, 0.0])
    pot = 0.
    
    p.x[n+1] = p.x[n] + dt * p.v[n] + (dt**2 * p.a[n])/2.0
    p.x[n+1] = p.x[n+1] % constants.box
    for op in particles:
        if op is p:
            continue
        
        move =  dt * op.v[n] + (dt**2 * op.a[n])/2.0
        op.x[n+1] = op.x[n] + move                                              # Modify the position
        op.d[n+1] += move                                                       #Modify the displacement from x0
        op.x[n+1] = op.x[n+1] % constants.box
        
        distvec = np.subtract(p.x[n+1] , op.x[n+1])
        
        
        
        dist = np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)          #Compute the distance between the two 
        # pointTo = np.linalg.norm(distvec)
        forcescalar =scaleBy * (2*(dist**-12.0) - dist**-6.0)                              # Compute the Lennard Jones potential
        acc_vec += forcescalar * (distvec / dist**2)                           # Add to the acceleration vector
        pot += (4*(dist**-12.0 - dist **-6.0)) 
    
    p.p[n+1] = pot
 
    p.a[n+1] = acc_vec
    p.v[n+1] = p.v[n] + dt * (p.a[n+1] + p.a[n])/2
    #p.x[n+1] = p.x[n] + dt * p.v[n] + (dt**2 * p.a[n+1])/2.0
    
    
    return p



def velocityVerletFastParalell(dt, n, particles, beginPindex, pCount):
   
    length= len(particles)
    scaleBy = 24
   
    
    #Move all the particles for this iteration
    for op in particles:
        move = dt * op.v[n] + (dt**2 * op.a[n])/2.0
        op.x[n+1] = op.x[n] + move
        op.d[n+1] += move
        op.x[n+1] = op.x[n+1] % constants.box


    i = 0
    j = beginPindex
    
    #COMPUTE THE FORCES BETWEEN THE PARTICLESzn
    #prepare an array to hold all pull and push vectors acting between particles.
    #acc_vec_ij = np.zeros((length,length,3))
    #print(f"Apply forces : ignore where range less than {cutOffRange}")
    while i < length:
            j = i+1
            p1 = particles[i]
            
            while j < pCount + beginPindex:
             
                p2 = particles[j]
                if i == j :
                    # The particle cannot interact with itself
                    continue
                
                distvec = np.subtract(p1.x[n+1] , p2.x[n+1])                               # Find the vector between them from i to j
                distvec = distvec - np.round(distvec / constants.box) * constants.box  # Compute the distance to the atom, but the shortest direction passing through periodic boundary conditions
                dist =  np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)
                if dist < cutOffRange and dist > 0 :
                    
                    potential = (4*(dist**-12.0 - dist**-6.0))
                    forceVector_ij =  ((2*(dist**-12.0) - dist**-6.0 - cutOffPotential) * (distvec / dist**2))  #Let the potential weaken with the square of the distance
                    
                    p1.a[n+1] +=  forceVector_ij * scaleBy                         # The pull for the j particle is always the opposite direction
                    p2.a[n+1] -=  forceVector_ij * scaleBy  
                    p1.p[n+1] += potential
                    p2.p[n+1] += potential
                j+=1
            i+=1
    
    i = 0  
    
    
    
    while i < length:
        p = particles[i]
        p.v[n+1] =  p.v[n] + dt * ((p.a[n+1] + p.a[n])/2)
        i+=1
    

def velocityVerletFast2(dt, n, particles, cutOffRange, cutOffPotential):
   
    length= len(particles)
    scaleBy = 24
   
    
    #Move all the particles for this iteration
    for op in particles:
        move = dt * op.v[n] + (dt**2 * op.a[n])/2.0
        op.x[n+1] = op.x[n] + move
        op.d[n+1] += move
        op.x[n+1] = op.x[n+1] % constants.box


    
    i = 0
    j = 0
    
    #COMPUTE THE FORCES BETWEEN THE PARTICLESzn
    #prepare an array to hold all pull and push vectors acting between particles.
    #acc_vec_ij = np.zeros((length,length,3))
    #print(f"Apply forces : ignore where range less than {cutOffRange}")
    while i < length:
            j = i+1
            p1 = particles[i]
            
            while j < length:
             
                p2 = particles[j]
                if i == j :
                    # The particle cannot interact with itself
                    continue
                
                distvec = np.subtract(p1.x[n+1] , p2.x[n+1])                               # Find the vector between them from i to j
                distvec = distvec - np.round(distvec / constants.box) * constants.box  # Compute the distance to the atom, but the shortest direction passing through periodic boundary conditions
                dist =  np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)
                if dist < cutOffRange and dist > 0 :
                    
                    potential = (4*(dist**-12.0 - dist**-6.0)) - cutOffPotential
                    forceVector_ij =  ((2*(dist**-12.0) - dist**-6.0 ) * (distvec / dist**2))  #Let the potential weaken with the square of the distance
                    
                    p1.a[n+1] +=  forceVector_ij * scaleBy                         # The pull for the j particle is always the opposite direction
                    p2.a[n+1] -=  forceVector_ij * scaleBy  
                    p1.p[n+1] += potential
                    p2.p[n+1] += potential
                j+=1
            i+=1
    
    i = 0  
    
     
    
    while i < length:
         p = particles[i]
         p.v[n+1] =  p.v[n] + dt * ((p.a[n+1] + p.a[n])/2)
        
         i+=1
    
    


##  Forward Euler
def euler(p, n,  dt, particles, eps, sig):
    scaleBy = 24
    acc_vec = np.array([0.0, 0.0, 0.0])
    pot = 0.
    #for op in particles:
    plen = len(particles)
    for op in particles:
        if op is p:
            continue
        
        distvec = np.subtract(p.x[n] , op.x[n])
        dist = np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)
        # pointTo = np.linalg.norm(distvec)
        forcescalar = 2*(dist**-12) - dist**-6 
        acc_vec = forcescalar * (distvec / dist**2)
        pot = MDFunctions.LennardJonesPotential(dist, eps, sig)
        
        #Initial particle
        p.p[n+1] += pot
        p.a[n+1] += scaleBy * acc_vec
        p.v[n+1] += p.v[n] + dt * p.a[n+1]
        p.x[n+1] += p.x[n] + dt * p.v[n]
        
        
   
    p.swap(n+1)


def eulerCromer(p, n,  dt, particles, eps, sig):
    scaleBy = 24
    acc_vec = np.array([0.0, 0.0, 0.0])
    pot = 0.
    for op in particles:
        if op is p:
            continue
        
        distvec = np.subtract(p.x[n] , op.x[n])
        dist = np.sqrt(distvec[0]**2 + distvec[1]**2 + distvec[2]**2)
        forcescalar = 2*(dist**-12) - dist**-6 
        acc_vec += forcescalar * (distvec / dist**2.0)
        pot += MDFunctions.LennardJonesPotential(dist, eps, sig)
    
    p.p[n+1] = pot
    p.a[n+1] = scaleBy * acc_vec
    p.v[n+1] = p.v[n] + dt * p.a[n+1]
    p.x[n+1] = p.x[n] + dt * p.v[n+1]  
    return p



def solveParalell(f, particles, T, N):
    cpus = mp.cpu_count()
    pool = mp.Pool(cpus)
    dt = T/N
    print(f"Beginning fast simulation {N} with delta {dt}")
    if dt > 0.1 :
        raise Exception("The delta t is too high, simulation wil fail")
   
          
    if not constants.cutOffRange :
        raise Exception("Cutoffrange must be specified in constants")
            
    if constants.eps is None or constants.sig is None :
        raise Exception("To specify a cutoff range, please specify Epsilon and Sigma in constants")
    
    cutOffPotential = 2*(constants.cutOffRange**-12.0) - constants.cutOffRange**-6.0
    print(f"Cutoff range {constants.cutOffRange} shifts potential to {cutOffPotential}")
    
    
    pCount = int(len(particles) / cpus)+1
    cpuParams =[]
    
    for i in range(0, cpus):
        cpuParams.append([i*pCount, (i+1)*pCount])
    cpuParams[cpus-1][1] = pCount - i*pCount - len(particles) 
    
    results = [pool.apply(f, args=(n, particles, row[0], row[1])) for row in cpuParams]
    
        
    for n in range(N-1):
        #if  n == 3:
            #print("Trap")
        print(f"Frame {n}")
        f(dt, n, particles, constants.cutOffRange, cutOffPotential, beginParticleIndex, endParticleIndex)
    print(f"End Fast simulation {N} with delta {dt}")
    pool.close()


def solveFaster(f, particles, T, N):
    dt = T/N
    print(f"Beginning fast simulation {N} with delta {dt}")
    if dt > 0.1 :
        raise Exception("The delta t is too high, simulation wil fail")
   
          
    if not constants.cutOffRange :
        raise Exception("Cutoffrange must be specified in constants")
            
    if constants.eps is None or constants.sig is None :
        raise Exception("To specify a cutoff range, please specify Epsilon and Sigma in constants")
    
    cutOffPotential = 2*(constants.cutOffRange**-12.0) - constants.cutOffRange**-6.0
    print(f"Cutoff range {constants.cutOffRange} shifts potential to {cutOffPotential}")
    
        
    for n in range(N-1):
        #if  n == 3:
            #print("Trap")
        print(f"Frame {n}")
        f(dt, n, particles, constants.cutOffRange, cutOffPotential)
    print(f"End Fast simulation {N} with delta {dt}")
    

def solve(f, particles,T, N, eps = 1, sig = 1):
    dt = T/N
    print(f"Beginning slow simulation {N} with delta {dt}")
    if eps is None or sig is None :
        raise Exception("To specify a cutoff range, please specify Epsilon and Sigma to after CutoffRange to compute the shifted potential")
    
    if len(particles)<= 0:
        raise Exception(f"Cannot simulate {len(particles)} particle system")
    
    for n in range(N-1):
        i = 0
        for part in particles:
            f(part, n, dt, particles, eps, sig)
            i += 1
    
    print(f"End slow simulation {N} with delta {dt}")
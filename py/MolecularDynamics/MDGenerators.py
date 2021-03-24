# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:03:43 2021

@author: jenst
"""

from MolecularDynamics import MDParticle
from MolecularDynamics import MDPlot
from MolecularDynamics import MDGenerators
from MolecularDynamics import MDSolver
from MolecularDynamics import MDFileWriter
from MolecularDynamics import MDFunctions

        
import numpy as np 
## Create particless
def createParticles(count, atomType, avgDistanceSigma, N , mass, temperature = None):
    particles  = []
    x0 = [0,0,0]
    if temperature is not None:
        v0 = MDFunctions.CreateRandomVelocityFromTemperature(count, temperature)
    else:
        v0 = np.array([[0., 0., 0.] for i in range(N)])  # Velocity
        
  

    xPositions = MDFunctions.getPositionInBox(count)
       
     # v0speedVector = MDFunctions.getRandomNormalizedDirection()
    for i in range(0, count):
        particle = MDParticle.Particle(i, atomType, xPositions[i], v0[i], N, mass)
        particles.append(particle)
    
    return particles


def PopulateCreateCrystalStructure(n, periodicCode, distance, N, mass, temperature = None):
    print(f"Computes {n}x{n}x{n} tetrahedral crystal lattice with cell size {distance}")
    faceCenteredCubicLattice = np.array([
        np.array([0, 0, 0]),                # corner
         np.array([0.5,0.5, 0]),            # floor
         np.array([0.5, 0, 0.5]),           # hither wall,
         np.array([0, 0.5, 0.5])            # tither wall
        ])
    
    
    particles  = []
    if temperature is not None:
        v0 = MDFunctions.CreateRandomVelocityFromTemperature(4*n**3, temperature)
    else:
        v0 = np.array([[0., 0., 0.] for i in range(4*n**3)])  # Velocity
    index = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                cellPos = np.array([i, j, k])* distance
                for lstr in faceCenteredCubicLattice:
                    particles.append(MDParticle.Particle(f"{i}-{j}-{k}", periodicCode, distance * lstr + cellPos, v0[index], N, mass))
                    index +=1
    
    print(f"Done popluating, created {index} particles")
    density = mass * 4 / distance**3
    print(f"Density of this structure is {density}")
    
    return particles


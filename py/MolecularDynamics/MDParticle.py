# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:01:09 2021

@author: jenst
"""

import numpy as np
from MolecularDynamics import constants
class Particle:
    
    def __init__(self, ids, typ,  x0, v0, N, mass):
      
        
        self.ids = ids
        self.mass = mass
        self.type = typ
    
     
        self.v = np.array([[0., 0., 0.] for i in range(N)])# Velocity
        self.x = np.array([[0., 0., 0.] for i in range(N)]) # Position
        self.a = np.array([[0., 0., 0.] for i in range(N)])# Acceleration
        self.d = np.array([[0., 0., 0.] for i in range(N)])# Displacement from initial position
        
        self.p = np.zeros(N) # Potential energy
        self.v[0] = v0  # Initial velocity
        self.x[0] = x0  # Initial position
        self.d[0] = [0., 0., 0.]  # Displacement vector from initialPosition
        self.a[0] = np.array([0., 0., 0.]) #
        print(f"Creating {self.type} particle initialized at {self.x[0]},, speed {self.v[0]}")
       
    ### Will reflect the particle of the wall
    def swap(self, n):
        #print (f"Current position f{self.x[n]}")
        self.x[n] = self.x[n] % constants.box
            ### Will reflect the particle of the wall
    
    def reflect(self, n):
        #print (f"Current position f{self.x[n]}")
        if self.x[n][1] > constants.box:
           print(f"Y reflection {self.x[n]} {self.v[n]}")
           self.x[n][1] = constants.box
           self.v[n][1] *= -1
           
        if  self.x[n][1] < -constants.box:
           print(f"Y reflection {self.x[n]} {self.v[n]}")
           self.x[n][1] = -constants.box
           self.v[n][1] *= -1
           
        if self.x[n][0] > constants.box:
           print(f"X reflection {self.x[n]} {self.v[n]}")
           self.x[n][0] = constants.box
           self.v[n][0] *= -1
           
        if  self.x[n][0] < -constants.box:
           print(f"X reflection {self.x[n]} {self.v[n]}")
           self.x[n][0] = -constants.box
           self.v[n][0] *= -1
           
        if self.x[n][2] > constants.box:
           print(f"Z reflection {self.x[n]} {self.v[n]}")
           self.x[n][2] = constants.box
           self.v[n][2] *= -1
           
        if  self.x[n][2] < -constants.box:
           print(f"Z reflection {self.x[n]} {self.v[n]}")
           self.x[n][2] = -constants.box
           self.v[n][2] *= -1
    
    ## Return mass times velocity squared for the time
    def getKinetic(i):
        return  self.mass * self.v[i]**2
    
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 12:31:17 2021

@author: jenstandstad
""" 
import numpy as np


from MolecularDynamics import MDParticle
from MolecularDynamics import MDPlot
from MolecularDynamics import MDFunctions
from MolecularDynamics import MDGenerators
from MolecularDynamics import MDSolver
from MolecularDynamics import MDFileWriter
from MolecularDynamics import constants


M=1                                         
T = 5
increments = 0.005
def SetTime(T, increments, beginTime = 0):
    
    N =  int((T) / increments)+1
    t = np.linspace(beginTime, T+beginTime, N)
    return t, N, T

t, N, T= SetTime(T, increments)

particles = []

basedt = T/N  # defining delta t  

def Task2b_i():
    print("TASK : Task2b_i")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])

    particles = [
        MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), 
        MDParticle.Particle('right', 'Ar', np.array([1.5, 0. , 0.]), v0, N, constants.mass)
        ]

    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)

def Task2b_ii():
    print("TASK : Task2b_ii")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    MDPlot.plotDistanceBetween(t,"Distance, sig = 1.5", particles[0], particles[1], "two-particle-distance-15") 

def Task2b_iv():
    print("TASK : Task2b_iv")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments) 
    v0 = np.array([0, 0, 0])
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([0.95, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    MDPlot.plotDistanceBetween(t, "Distance, sig = 0.95" , particles[0], particles[1], "two-particle-distance-095") 

def Task2c_i():
    print("TASK : Task2c_i")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Euler Cromer) 1.5 sigma","2ci_allEnergiesCromer095")
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([0.95, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Euler Cromer) 0.95 sigma", "2ci_allEnergiesCromer15")




    
def Task2c_iv():
    print("TASK : Task2c_iv")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)  
    MDSolver.solve(MDSolver.euler, particles, T, N)
    MDPlot.plotAllEnergies(particles, t,"Total Energy (Euler) 1.5 sigma", "2civ_allEnergiesEuler15")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesEuler15", N)
    
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]   
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Euler Cromer) 1.5 sigma" , "2civ_allEnergiesEulerCromer15")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesEulerCromer15", N)
    

    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.velocityVerlet, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 1.5 sigma",  "2civ_allEnergiesVelocityVerlet15")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVer15", N)
    
       
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([0.95, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Euler) 0.95 sigma", "2civ_allEnergiesEuler095")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesEuler095", N)
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([0.95, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDFileWriter.WriteXYZFile(particles,"twoparticlesEulerCromer095", N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Euler Cromer) 0.95 sigma", "2civ_allEnergiesEulerCromer095")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesEulerCromer095", N)
   
    
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([0.95, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.velocityVerlet, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 0.95 sigma",  "2civ_allEnergiesVelocityVerlet095")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVer095", N)
    
    
    
def Task3a_iiii():
    print("TASK : Task3a_iiii")
    print("Verify the shifted potential")
    MDPlot.plotPotential(0.9, 5, "Plot potential (shifted and unmodified)",  "plotPotential")
    
    MDPlot.plotShiftedPotentialDiff(0.90, 5, "Plot potential (shifted and unmodified)",  "plotPotentialDiff")

def Task3b_i():
    print("TASK : Task3b_i")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    print("Verification by reproduction of the two atom simulation")
    ## Create particless
    #createParticles(3, "AR", 1, particles, T, N)
   
    # solveFaster(velocityVerletFast2, particles, T, N, constants.cutOffRange)
    #WriteXYZFile(particles,"threeparticlesVelFast")
    
        
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solve( MDSolver.velocityVerlet, particles, T, N)
    #MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 1.5 sigma",  "allEnergies15velocityVerlet15")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVer15", N)
       
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solveFaster( MDSolver.velocityVerletFast2, particles, T, N)
    #MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet)(fast) 1.5 sigma",  "allEnergies15velocityVerlet15fast")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVer15Fast", N)
    

def Task3b_ii():
    print("TASK : Task3b_ii")
    maxTime = 5
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    print("Simulating four atom")
    ## Create particless

    particles = [
        MDParticle.Particle('one',     'Ar', np.array([1.,0.0,0.] ), v0, N, constants.mass),
        MDParticle.Particle('two',     'Ar', np.array([0., 1. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('three',   'Ar', np.array([-1., 0. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('four',    'Ar', np.array([0.,-1. , 0.]), v0, N, constants.mass)
        ]

    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet)(fast) 1.0 sigma",  "allEnergies15velocityVerlet15fast")
    MDFileWriter.WriteXYZFile(particles,"fourparticlesVelVer1Fast", N)


def Task3b_iv():
   print("TASK : Task3b_iv")
   maxTime = 5
   t, N, T = SetTime(maxTime, increments)
   v0 = np.array([0, 0, 0])
   particles = [
        MDParticle.Particle('one',     'Ar', np.array([1.,0.0,0.] ), v0, N, constants.mass),
        MDParticle.Particle('two',     'Ar', np.array([0., 1. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('three',   'Ar', np.array([-1., 0. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('four',    'Ar', np.array([0.,-1. , 0.]), v0, N, constants.mass)
    ]
   MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, T, N)
   MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 4 particles",  "3b_iv")
   
   MDFileWriter.WriteXYZFile3(particles, "3b_iv", N)


    


def Task3b_v():
    print("TASK : Task3b_v")
    print("Simulating four atom")
    nthFrame = 10
    v0 = np.array([0, 0, 0])
        
    particles = [
        MDParticle.Particle('one',     'Ar', np.array([1.,0.1,0.] ), v0, N, constants.mass),
        MDParticle.Particle('two',     'Ar', np.array([0., 1. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('three',   'Ar', np.array([-1., 0. , 0.]), v0, N, constants.mass),
        MDParticle.Particle('four',    'Ar', np.array([0.,-1. , 0.]), v0, N, constants.mass)
        ]

    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 4 particles, perturbed",  "3b_v")
    MDFileWriter.WriteXYZFile(particles,"fourparticlesVelVer1FastPerturb", N, nthFrame)



def Task3c_i():
    print("TASK : Task3c_i")
    maxTime = 3 #should be 3
    
    t, N, T = SetTime(maxTime, increments)
    nthFrame = 10
    particles = MDGenerators.PopulateCreateCrystalStructure(3,"AR", 1, N, constants.mass)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDFileWriter.WriteXYZFile(particles,"crystal", N, nthFrame)



def Task3c_ii():
    print("TASK : Task3c_ii")
    maxTime = 3 #should be 3
    t, N, T = SetTime(maxTime, increments)
    nthFrame = 10
    particles = MDGenerators.PopulateCreateCrystalStructure(3,"AR", 1, N, constants.mass, 300)
    
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDFileWriter.WriteXYZFile(particles,"crystal", N, nthFrame)
    print("Finished")


def Task3d_i():
    print("TASK : Task3d_i")
    maxTime = 3 #should be 3
    t, N, T = SetTime(maxTime, increments)
    nthFrame = 10
    particles = MDGenerators.PopulateCreateCrystalStructure(4, "AR", 1.7, N, constants.mass)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 256 particles, random",  "3d_ii")
    MDFileWriter.WriteXYZFile(particles, "256_random", N, nthFrame)
   

def Task3d_ii():
    print("TASK : Task3d_i")
    maxTime = 5 #should be 3
    t, N, T = SetTime(maxTime, increments)
    nthFrame = 10
    dist = 1.7
    v0 = np.array([0, 0, 0])
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.,0.,0.] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([dist, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDFunctions.MoveToCentreOfBox(particles)
    
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 1.5 sigma",  "3dii_energy_two")
    MDFileWriter.WriteXYZFile(particles,"3dii_two", N)

    
    
    particles = MDGenerators.PopulateCreateCrystalStructure(1, "AR", dist, N, constants.mass)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 4 particles, random",  "3dii_energy_four")
    MDFileWriter.WriteXYZFile(particles, "3dii_four", N, nthFrame)
   


def Task3e_ii():
    print("TASK : Task3e_ii")
    maxTime = 1 #should be 3
    nthFrame = 1
    t, N, T = SetTime(maxTime, increments)
    particles = MDGenerators.PopulateCreateCrystalStructure(2, "AR", 1.7,  N, constants.mass )
    MDFunctions.MoveToCentreOfBox(particles)
    #MDFileWriter.WriteXYZFile(particles, "108_initial", N, nthFrame)

    #MDSolver.solve(MDSolver.velocityVerlet, particles, T, N)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    MDFileWriter.WriteXYZFile(particles, "108_initial", N, nthFrame)
     

def Task4a_ii():
    print("TASK : Task4a_ii")
    maxTime = 10
    t, N, T = SetTime(maxTime, increments)
    t_initial = 300
    particles = MDGenerators.PopulateCreateCrystalStructure(3, "Ar", 1.7,  N, constants.mass, t_initial)
    MDFileWriter.WriteFirstFrame(particles,"firstFrame", N)
    #MDFunctions.MoveToCentreOfBox(particles)
    #MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDSolver.solve(MDSolver.velocityVerlet, particles, maxTime, N, constants.eps, constants.sig)
    MDFileWriter.WriteXYZFile(particles,"checkthis_slow", N)
    chartName = f"4a) ii Temp ({t_initial} K) (slow)"
    MDPlot.plotTemperature(particles, t, N,  chartName, "4a_ii")

def Task4a_ii_fast():
    print("TASK : Task4a_ii")
    maxTime = 1.5
    t, N, T = SetTime(maxTime, increments)
    t_initial = 300
    particles = MDGenerators.PopulateCreateCrystalStructure(3, "Ar", 1.7,  N, constants.mass,  t_initial )
    MDFileWriter.WriteFirstFrame(particles,"firstFrame", N)
    #MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    #MDSolver.solve(MDSolver.velocityVerlet, particles, maxTime, N, constants.eps, constants.sig)
    MDFileWriter.WriteXYZFile(particles,"checkthis_fast", N)
    chartName = f"4a) ii Temp ({t_initial} K) (fast)"
    MDPlot.plotTemperature(particles, t, N, chartName, "4a_ii_fast")



def Task4a_iii():
    print("TASK : Task4a_iii")
    maxTime = 10
    t_initial = 250
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(3, "Ar", 1.7,  N, constants.mass, t_initial)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.velocityVerlet, particles, maxTime, N, constants.eps, constants.sig)
    
    #temp = MDFunctions.ComputeTemperatureReduced(N-1, particles, constants.eps)
    #print(f"AT time {maxTime} temperature is {temp} kelvin")
    chartName = f"4a) iii Temp ({t_initial} K) (slow)"
    MDPlot.plotTemperature(particles, t, N,  chartName, "4a_iii_slow")


def Task4a_iii_fast():
    print("TASK : Task4a_iii")
    maxTime = 10
    t, N, T = SetTime(maxTime, increments)
    t_initial = 190
    particles = MDGenerators.PopulateCreateCrystalStructure(3, "Ar", 1.7,  N, constants.mass, t_initial)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    #temp = MDFunctions.ComputeTemperatureReduced(N-1, particles, constants.eps)
    #print(f"AT time {maxTime} temperature is {temp} kelvin")
    chartName = f"4a) iii Temp ({t_initial} K) (fast)"
    MDPlot.plotTemperature(particles, t, N, chartName, "4a_iii_fast")


def Task4b_simAndWrite():
    print("TASK : Task4b_ii_write")
    maxTime = 1.5
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(4, "Ar", 1.7,  N, constants.mass, 190)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDFileWriter.WriteEverything(particles, "4b_ii_everything", N)

def Task4b_ii():
    print("TASK : Task4b_ii")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    
    #particles = MDGenerators.PopulateCreateCrystalStructure(4, "Ar", 1.7,  N, constants.mass, 190)
    #MDFunctions.MoveToCentreOfBox(particles)
    #MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    particles = MDFileWriter.ReadEverything("4b_ii_everything", N)
    MDPlot.plotVelocityAutoCorrelation(particles, t, N,  "4b) ii", True, "4b_ii")
    MDFileWriter.WriteXYZFile(particles, "4b_ii", N)
    MDFileWriter.WriteLastFrame(particles, "4b_ii_lastFrame", N)
    


def Task4b_iii():
    print("TASK : Task4b_iii")
    beginTime = 0.7
    maxTime = 1.5
    t, N, T = SetTime(maxTime, increments, beginTime)
    
    particles = MDFileWriter.ReadInitialConditions("4b_ii_lastFrame", N, constants.mass)
    #MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDPlot.plotVelocityAutoCorrelation(particles, t , N,  "4b) iii - Velocity autocorrelation plot", "4b_iii")
    MDFileWriter.WriteXYZFile(particles, "4b_iii", N)
    
    MDFileWriter.WriteLastFrame(particles, "4b_iii_lastFrame", N)

def Task4b_iv_avg():
    print("TASK : Task4b_iv_avg")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    vac = []
    sampling = 5
    for i in range(0, sampling):
        particles = MDGenerators.PopulateCreateCrystalStructure(4, "Ar", 1.7,  N, constants.mass, 190)
        MDFunctions.MoveToCentreOfBox(particles)
        MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
        vag, diffu_coeff = MDFunctions.ComputeVelocityAutoCorrelation(t, particles)
        vac.append( vag)
        
    avgVag = vac[0]
    for i in range(0, sampling):
        avgVag += vac[i]
    
    avgVag[0] /= sampling
    chartName = f"Avg Vel AUtocorrelation {sampling}"
    MDPlot.plot([avgVag], t, N, chartName, ["vac"], "time (s)", "4biv_avg")
    MDFileWriter.WriteXYZFile(particles, "4b_ii", N)
    MDFileWriter.WriteLastFrame(particles, "4b_iiii_lastFrame", N)


def Task4b_v():
    print("TASK : Task4b_v")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(4, "Ar", 1.7,  N, constants.mass, 190)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    MDPlot.plotVelocityAutoCorrelation(particles, t, N, "4b) ii - Vel.Corr / Diff.coeff ", "4b_v")
    MDFileWriter.WriteXYZFile(particles, "4b_v", N)
    MDFileWriter.WriteLastFrame(particles, "4b_v_lastFrame", N)
    
def Task4c_i():
    print("TASK : Task4c_i")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(1, "Ar", 1.7,  N, constants.mass, 190)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    MDPlot.plotMeanSquaredDisplacement(particles, t, N, constants.eps, "4c) i - Mean squared displacement", False, "4c_i")
    MDFileWriter.WriteXYZFile(particles, "4c_i", N)
    MDFileWriter.WriteLastFrame(particles, "4c_i_lastFrame", N)
    
  
def Task4c_ii():
    print("TASK : Task4c_ii")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(6, "Ar", 1.7,  N, constants.mass, 190)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    
    MDPlot.plotMeanSquaredDisplacement(particles, t, N, constants.eps,  "4c) i - Mean squared displacement", True, "4c_ii")
    MDFileWriter.WriteXYZFile(particles, "4c_ii", N)
    MDFileWriter.WriteLastFrame(particles, "4c_ii_lastFrame", N)
    


def Task2FasterCompare():
    print("TASK : Task2FasterCompare")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solve(MDSolver.velocityVerlet, particles, T, N)
    #MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 0.95 sigma",  "allEnergies15velocityVerlet95")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVerSlow", N, 50)
    
    
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0 , 0.0]), v0, N, constants.mass)]
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, T, N)
   # MDPlot.plotAllEnergies(particles, t, "Total Energy (Velocity Verlet) 0.95 sigma",  "allEnergies15velocityVerlet95")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVerFast", N, 50)


def Task2FasterCompare2():
    print("TASK : Task2FasterCompare2")
    maxTime = 0.7
    t, N, T = SetTime(maxTime, increments)
    v0 = np.array([0, 0, 0])
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar', np.array([1.5, 0.0, 0.0]), v0, N, constants.mass)]
    
    #particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass)]
    
    MDFunctions.MoveToCentreOfBox(particles)
    
    MDSolver.solve(MDSolver.velocityVerlet, particles, T, N)
    #MDSolver.solve(MDSolver.eulerCromer, particles, T, N)
   
    MDPlot.plotAllEnergies(particles, t, "Test slow",  "testSlow")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVerSlow", N, 50)
    
    
    
    particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass), MDParticle.Particle('right', 'Ar',  np.array([1.5, 0.0, 0.0]), v0, N, constants.mass)]
    #particles = [MDParticle.Particle('left', 'Ar', np.array([0.0,0.0,0.0] ), v0, N, constants.mass)]
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, T, N)
    MDPlot.plotAllEnergies(particles, t, "Test fast",  "testFast")
    MDFileWriter.WriteXYZFile(particles,"twoparticlesVelVerFast", N, 50)
  
    
def Task4d_i_write():
    
    print("TASK : Task4b_ii_write")
    maxTime= 3
    t, N, T = SetTime(maxTime, increments)
    
    particles = MDGenerators.PopulateCreateCrystalStructure(6, "Ar", 1.7,  N, constants.mass, 190)
    MDFunctions.MoveToCentreOfBox(particles)
    MDSolver.solveFaster(MDSolver.velocityVerletFast2, particles, maxTime, N)
    MDFileWriter.WriteEverything(particles, "Task4_di_everything", N)
    

def Task4_di_read(t  = None):
    num_bins = 15
  
    bin_edges =  np.linspace(0, constants.box, num_bins+1)
    bin_centres = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    V = constants.box**3
    
    #Read the file prepared by the previous 
    particles, maxTindex = MDFileWriter.ReadEverything("Task4_di_everything", N)
    
    #Convert x-axis units from sigma to Ångstrøm
    bin_centres = MDFunctions.SigmaToAngstrom(bin_centres)
    
    if t is None:
        t = maxTindex
    elif t > maxTindex:
        raise Exception("Trying to read out rdf from frame {t} which only simulated {maxTindex} frames")
    
    tInSeconds = MDFunctions.ArgonTimeToSeconds(t)
   
    count = len(particles)
    x = np.array([[0., 0., 0.] for i in range(count)])# Velocity
   
   
   
    j = 0
    counter = 0

    for j in range(0, count):
         x[counter] = particles[j].x[t]
         counter +=1
    
    #x = np.where(x > 0)
        
    #Bin the distance counts
    data = MDFunctions.rdf(bin_edges, x, V)
    
    
    chartName = f"Ar N={count} RDF after {tInSeconds[0]:11.6} {tInSeconds[1]}"
    # def barchart(data,  xlabels, chartName, filename = "barChart"):
    MDPlot.barchart(data, bin_centres, "RDF", chartName, "4di_rdf_bar")
    
    strLabels = []
    for lb in xLabel[0]:
        strLabels.append(f"{lb:11.1f}")
   
    MDPlot.plot(data, bin_centres, count, chartName,strLabels, "4di_rdf_plot")         
    

def TemperatureStuff():
    bltz = constants.bltz
    print(bltz)
    print(f"Temperature {1/bltz}")
    
#Task2FasterCompare2()

# Task2b_i()    #OK
# Task2b_ii()     #OK
# Task2b_iv()   #OK
# Task2c_i()      #OK
# Task2c_iv()   #OK

# Task3a_iiii()
# Task3b_i()
# Task3b_ii()
# Task3b_iv()
# Task3b_v()


# Task3c_i()
# Task3c_ii()
# Task3d_i()
# Task3d_ii()

# Task3e_ii()
#Task4a_ii_fast()

# Task4a_iii_fast()
# Task4b_simAndWrite()
# Task4b_ii()
# Task4b_v()
# Task4b_iii()
# Task4b_iv_avg()
#Task4c_i()
#Task4c_ii()


Task4d_i_write()
Task4_di_read()
#TemperatureStuff()

#TwoParticleTest()
#for p in particles:
#    plotParticle(t, p)
    



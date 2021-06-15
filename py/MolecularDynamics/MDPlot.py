# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:09:52 2021

@author: jenst
"""

import numpy as np
import matplotlib.pyplot as plt
import math as math


from MolecularDynamics import MDFunctions
from MolecularDynamics import constants




def plot(series, t, N,  chartName, xlabels, yLabel, filename = "genericPLot"):
    if len(series) > 0:
        if len(series) != len(xlabels):
            raise Exception("Generic plot requires label array of same length as series array")
    elif  len(series)==0:
        raise Exception("Generic plot requires at least one series")
        
    plt.title(chartName)
    for i in range(0, len(series)):
        plt.plot(t, series[i], label=xlabels[i])
        
 
    plt.ylabel(yLabel)
    plt.legend()
    
    plt.savefig(f'./img/{filename}.jpg')
    plt.show()



def barchart(data,  xLabel, yLabel, chartName, filename = "barChart"):
    print(f"Making barchart {chartName}")
    #objects = ('Python', 'C++', 'Java', 'Perl', 'Scala', 'Lisp')
    #y_pos = np.arange(len(objects))
    #performance = [10,8,6,4,2,1]
    strLabels = []
    for lb in xLabel[0]:
        strLabels.append(f"{lb:11.1f}")
    plt.bar(strLabels, data, align='center')
    plt.title(chartName)
   
    #plt.xticks(xLabel[0], strLabels)
   
    plt.ylabel(yLabel)
    plt.xlabel(xLabel[1])
    
    
  
    
    print(f"Saving file {filename}.jpg")
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()
   

def plotMeanSquaredDisplacement(particles, t, N, eps, chartName, plotDiffusionConstant = False,  filename = "meanSquaredDisplamcent"):
    

    msd, diffusion_constant = MDFunctions.ComputeMeanSquaredDisplacement(t, particles)
  
    msd, yLabel = MDFunctions.SigmaToAngstrom(msd)
    t, xLabel = MDFunctions.ArgonTimeToSeconds(t)
    
    plt.title(chartName)
    plt.plot(t, msd, label=f'Mean squared displacement')
    plt.xlabel(xLabel)
    plt.ylabel(f"Mean sq. disp {yLabel}")
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()
    if plotDiffusionConstant :
        plt.title("Diffusion constant")
        plt.plot(t, diffusion_constant, label=f'Diffusion constant')
        plt.xlabel(xLabel)
        plt.ylabel('D')
        plt.legend()
        
        plt.savefig(f'{filename}_diffusionConstant.jpg')
        plt.show()
        
def plotTemperature(particles, t, N, chartName, filename = "plotTemperature"):
    
    
    tPrime = MDFunctions.ComputeTemperatureReduced(N, particles)    
    
    tKelvin = MDFunctions.ArgonTPrimeToKelvin(tPrime)
    
    tToSeconds, xLabel = MDFunctions.ArgonTimeToSeconds(t)
        
    plt.title(chartName)
    plt.plot(tToSeconds, tKelvin, label=f'Temperature (K)')
    plt.xlabel(xLabel)
    plt.ylabel('Kelvin (K)')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()
    

def plotVelocityAutoCorrelation(particles, t, N, chartName, includeDiffusionCoefficientSubPlot = False, filename = "velocityAutoCorrelation"):
    

    vac, diffusion_coefficcient = MDFunctions.ComputeVelocityAutoCorrelation(t, particles)
    
    t, timeUnits= MDFunctions.ArgonTimeToSeconds(t)
        
    plt.subplot(1, 1, 1)
    plt.title(f"{chartName} Velocity correlation")
    plt.plot(t, vac, label=f'Velocity auto correlation')
    plt.xlabel(f"time ({timeUnits})")
    plt.ylabel('Vel. Autocorrelation')
    plt.legend()
  
    plt.savefig(f'{filename}_velcorr.jpg')
    plt.show()
    if includeDiffusionCoefficientSubPlot:
        plt.subplot(2, 1, 1)
        plt.title(f"{chartName} Diffusion Coefficient")
        plt.plot(t, diffusion_coefficcient, label=f'Diffusion coefficient')
        plt.xlabel('time(s)')
        plt.ylabel('Diffusion coefficient')
   
    plt.legend()
    
    plt.savefig(f'{filename}_diffcoeff.jpg')
    plt.show()
    
def plotTemperature(particles, t, N, chartName, filename = "plotTemperature"):
    
    
    tPrime = MDFunctions.ComputeTemperatureReduced(N, particles)    
    
    tKelvin = MDFunctions.ArgonTPrimeToKelvin(tPrime)
    
    tToSeconds, xLabel= MDFunctions.ArgonTimeToSeconds(t)
        
    plt.title(chartName)
    plt.plot(tToSeconds, tKelvin, label=f'Temperature (K)')
    plt.xlabel(xLabel)
    plt.ylabel('Kelvin (K)')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()
    print(f"Finished plotting temperature filename {filename}")
    
    
## Plot all energies
def plotAllEnergies(particles, t, chartName, filename = "allEnergies"):
  
    kinetic = np.zeros(len(t))
    potential = np.zeros(len(t))
    total = np.zeros(len(t))
    for i in range(0, len(t)-1):
        for p in particles:
            kinetic[i] += ( MDFunctions.getLength(p.v[i])**2 )/2
            potential[i] += p.p[i] / 2
            
    
  
    total = np.add(kinetic , potential)
    
    plt.title(chartName)
    plt.plot(t, total, label=f'Total ')
    plt.plot(t, potential, label=f'Potential')
    plt.plot(t, kinetic, label=f'Kinetic')
    plt.xlabel('time(s)')
    plt.ylabel('energy')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()
    
    
       
## Plot all energies
def plotShiftedPotentialDiff(minDist, maxDist,  chartName, filename = "plotPotentialDiff"):
  
    shift = MDFunctions.LennardJonesPotential(3, constants.eps, constants.sig)
    
    series = np.linspace(2.8, 3.1, 100)
    pot =  MDFunctions.LennardJonesPotential(series, constants.eps, constants.sig, 0)
    pot2 = MDFunctions.LennardJonesPotential(series, constants.eps, constants.sig, shift)
    diff = pot - pot2 
 
    plt.subplot(1,1,1)
    plt.title(chartName)
    plt.plot(series, pot, label=f'Non shifted')
    plt.plot(series, pot2, label=f'Shifted')
    plt.plot(series, diff, label=f'Difference between shifted and non shifted')

    plt.xlabel('distance (sigmas)')
    plt.ylabel('force')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()


## Plot the unmodified and shifted potentail
def plotPotential(minDist, maxDist,  chartName, filename = "plotPotential"):
  

    series = np.linspace(minDist, maxDist, 100)
    pot = MDFunctions.LennardJonesPotential(series, constants.eps, constants.sig)
    shiftedPot = MDFunctions.LennardJonesPotential(series, constants.eps, constants.sig, constants.cutOffRange)
    
    
    #plt.subplot(1,1,1)
    plt.title(chartName)
    plt.plot(series, pot, label=f'Unmodified')
    plt.plot(series, shiftedPot, label=f'Shifted')
    plt.xlabel('distance (sigmas)')
    plt.ylabel('force')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    plt.show()

## Compute distance between two particles and plot
def plotDistanceBetween(t, chartName, p1, p2, filename = "two-particle-distance"):
    N = len(t)
    dist = np.zeros(N)
    distanceVectors = np.subtract(p1.x, p2.x)
    for i in range(N):
        dist[i] = MDFunctions.getDistance(p1.x[i], p2.x[i])
    
    #plt.subplot(3,1,1)
    plt.title(chartName)
    plt.plot(t, dist, label=f'x(t) ')
    plt.xlabel('time(s)')
    plt.ylabel('distance (sigmas)')
    plt.legend()
    
    plt.savefig(f'./fig/{filename}.jpg')
    
## Plot the x component of a specific particle 
def plotParticle(t, p):
    plt.subplot(3,1,1)
    plt.title(f"Particle p {p.ids}")
    plt.plot(t, p.x, label=f'x(t) - {p.ids}')
    plt.xlabel('time(s)')
    plt.ylabel('x(t)')
    plt.legend()
    
    
    plt.subplot(3,1,2)
    plt.plot(t, p.v, label=f'v(t) - {p.ids}')
    plt.xlabel('time(s)')
    plt.ylabel('v(t)')
    plt.legend()
    
    
    plt.subplot(3,1,3)
    
    plt.plot(t, p.a, label=f'a(t) - {p.ids}')
    plt.xlabel('time(s)')
    plt.ylabel('a(t)')
    plt.legend()
    
 
    plt.savefig(f'Paricle-{p.ids}.jpg')
    
    plt.show()

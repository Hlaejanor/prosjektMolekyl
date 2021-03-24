# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:48:02 2021

@author: jenst
"""

bltz =  1.380649e-23 
eps = 1
eps2 = 119.7  * bltz 
sig = 1
mass =1


#Physical constants metric
sigmaMetric= 3.405E-10 
sigmaAnsgtrom = 3.405 
sigmaScaler = sig / sigmaMetric
massMetric = 1.66E-27 * 39.95
massScaler = mass / massMetric
epsMetric = 1.0318E-2 *  1.602E-19 #Specificed in joules
epsScaler = eps / epsMetric
epsEV = 1.0318 * 10**(-2) #Argon Electron Volt
kBoltz = 8.6173*10**(-5)


box = 40

#Simulation specifics
cutOffRange = 3* sig #
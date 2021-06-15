# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:34:52 2021

@author: jenst
"""
from MolecularDynamics import MDParticle
## Write  YZ files


def UtoXYZFileBulk(u, gobblers, filename, N, nthFrame  = 1):
        print(f"Writing universe bulk file: {filename} {N} frames")
        if nthFrame > 1:
            print(f"Writing : {nthFrame}-th frames")
        frame = 0
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            size = len(u)
            while frame < N:
                
                xyz_file.write(f"{len(gobblers)}\n")
                xyz_file.write(f"Frame{frame}\n")
                for i in range(0, size):
                    for j in range(0, size):
                        for k in range(0, size):
                            que =  u[i][j][k]
                            if que is not None:
                                xyz_file.write(f"{que.getType()} {i:11.6f} {j:11.6f} {k:11.6f}\n")
                frame += nthFrame                    
                    


def UtoXYZFile(gobblers, filename, N, nthFrame  = 1):
        print(f"Writing universe by gobblers file: {filename} {N} frames")
        if nthFrame > 1:
            print(f"Writing : {nthFrame}-th frames")
        frame = 0
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            size = len(u)
            while frame < N:
                
                xyz_file.write(f"{len(gobblers)}\n")
                xyz_file.write(f"Frame{frame}\n")
                for g in gobblers:
                    q = g.currentQue
                    while q is not None:
                        xyz_file.write(f"{que.getType()} {i:11.6f} {j:11.6f} {k:11.6f}\n")
                frame += nthFrame                    
                  
            




def WriteXYZFile2(particles, filename, N, nthFrame  = 1):
        print(f"Writing file: {filename} {N} frames")
        if nthFrame > 1:
            print(f"Writing : {nthFrame}-th frames")
        frame = 0
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            while frame < N:
               
                xyz_file.write(f"{len(particles)}\n")
                xyz_file.write(f"Frame{frame}\n")
                for part in particles:
                    x, y, z = part.x[frame]
                    xyz_file.write(f"{part.type} {x:11.6f} {y:11.6f} {z:11.6f}\n")
                frame += nthFrame                    
                    
                
            


## Write  YZ files
def WriteXYZFile(particles, filename, N, nthFrame = 1):
        print(f"Writing file: {filename}")
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            frame  = 0
            while frame < N:
                xyz_file.write(f"{len(particles)}\n")
                xyz_file.write(f"F{frame}\n")
                for part in particles:
                    x, y, z = part.x[frame]
                    xyz_file.write(f"{part.type}  {x:11.6f}  {y:11.6f}  {z:11.6f}\n")
                frame += nthFrame
            
 
         
## Write  YZ files
def WriteXYZFile3(particles, filename, N):
        print(f"Writing last frame: {filename}")
        count = len(particles)
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            for i in range(0, 2):
                xyz_file.write(f"{count}\n")
                xyz_file.write("Frame 0\n")
                for part in particles:
                    x, y, z = part.x[0]
                    xyz_file.write(f"{part.type}  {x:11.6f}  {y:11.6f}  {z:11.6f}\n")
                
   
                     
 
## Write  YZ files
def WriteLastFrame(particles, filename, N):
        print(f"Writing last frame: {filename}")
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
     
            xyz_file.write(f"{len(particles)}\n")
            xyz_file.write(f"Comments for frame {N}\n")
            for part in particles:
                x, y, z = part.x[N-1]
                vx, vy, vz = part.v[N-1]
                xyz_file.write(f"{part.type}  {x:11.6f}  {y:11.6f}  {z:11.6f}  {vx:11.6f}  {vy:11.6f}  {vz:11.6f}\n")



 
## Write  YZ files
def WriteEverything(particles, filename, N):
        print(f"Writing everythin : {filename}")
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
            for i in range(0, N-1):
                xyz_file.write(f"{len(particles)}\n")
                xyz_file.write(f"Frame {i}\n")
                for part in particles:
                    x, y, z = part.x[i]
                    vx, vy, vz = part.v[i]
                    ax, ay, az = part.a[i]
                    pot = part.p[i]
                    mass = part.mass
                  
                    xyz_file.write(f"{part.type}   {part.mass}   {x:11.6f}  {y:11.6f}  {z:11.6f}  {vx:11.6f}  {vy:11.6f}  {vz:11.6f}  {ax:11.6f}  {ay:11.6f}  {az:11.6f}    {pot:11.6f}\n")

 
## Write  YZ files
def WriteFirstFrame(particles, filename, N):
        print(f"Writing first frame: {filename}")
        with open(f"./xyz/{filename}.xyz", 'w') as xyz_file:
          
                xyz_file.write(f"{len(particles)}\n")
                xyz_file.write(f"Comments for frame {N}\n")
                for part in particles:
                    x, y, z = part.x[0]
                    vx, vy, vz = part.v[0]
                  
                    xyz_file.write(f"{part.type}  {x:11.6f}  {y:11.6f}  {z:11.6f}  {vx:11.6f}  {vy:11.6f}  {vz:11.6f} {pot:11.6f}\n")




## Write  YZ files
def ReadEverything(filename, N, mass = 1):
    particles = []
    index = 0
    
    particles = ReadInitialConditions(filename, N, mass)
    
    print(f"Reading remaining positions, velocities and acceleration from : {filename}")
    i = 0
    t =0
    with open(f"./xyz/{filename}.xyz", 'r') as xyz_file:
        for line in xyz_file:
            line_data = line.split()
            values = len(line_data)
            if values == 1:
                print(f"Everything : Reading frame {t}")
                #first line we have
                partCount = int(line_data[0])
                i = 0
                t += 1
            if values == 12:
                # print(f"   F {t} particle {i}" )
                 atom, mass,  x, y, z, vx, vy, vz , ax, ay, az, pot = line_data
                 
                 particles[i].a[t] = [ ax, ay, az]
                 particles[i].mass = mass
                 particles[i].v[t] = [vx, vy, vz]
                 particles[i].x[t] = [x, y, z]
                 particles[i].p[t] = pot
                 i +=1
                 
             
                
                 
         
    return particles, t


## Write  YZ files
def ReadInitialConditions(filename, N, mass):
    particles = []
    index = 0
    startedReadingFirstFrame = False
    print(f"Reading positions fomr (will stop after first frame): {filename}")
    with open(f"./xyz/{filename}.xyz", 'r') as xyz_file:
        for line in xyz_file:
            line_data = line.split()
            values = len(line_data)
            if len(line_data) == 1:
                if startedReadingFirstFrame:
                     print(f"Finished reading first frame {filename}")
                     return particles
            elif len(line_data) == 7:

                 atom, x, y, z, vx, vy, vz = line_data
                 particles.append(MDParticle.Particle(index, atom, [x, y, z], [vx,vy,vz], N, 1))
            if values == 12:
                startedReadingFirstFrame = True
                atom, mass,  x, y, z, vx, vy, vz , ax, ay, az, pot = line_data
                part = MDParticle.Particle(index, atom, [x, y, z], [vx,vy,vz], N, 1)
                part.mass = mass
                
                particles.append(part)
    
         
    return particles
    
       
          
         
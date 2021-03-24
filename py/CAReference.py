# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 15:48:07 2021

@author: jenst
"""
import matplotlib.pyplot as plt
import numpy as np

from MolecularDynamics import MDFileWriter
from sys import exit

n = 100
class Que:
    
    def __init__(self, prevQue, a, b, u, T, initialPos = None, queName = ""):
        self.H = None
        self.queName = queName
        self.Ataken = True
        self.Btaken = True
        self.Htaken = True
        self.x = np.array([[0., 0., 0.] for i in range(T)]) # Position
        self.nextQue = None
        if prevQue is None and initialPos is not None:
            self.pos = initialPos
        else:
            self.prevQue = prevQue
            self.pos = prevQue.pos + prevQue.A
            prevQue.nextQue = self
        
        if a is None:
            raise Exception("attempting to set A as None in constructor")
        
        self.setA(a)
        self.setB(b)
     
        #Fix this to the worldsheet for collission detecting
        x,y, z = self.pos
        print(f"Creating que objet at {self.pos}")
        if u[x][y][x] != None :
            raise Exception(f"Cannot create queue at {[x, y, z]}, already occupied slot in universe")
        u[x][y][z] = self
    def printRecursive(self):
        if self.nextQue is None:
         
            print(f"HEAD {self.queName} A: {self.readA()} B:{self.readB()} H : {self.readH() }p: {self.pos}")
        else:
            self.nextQue.printRecursive()
            print(f"BODY {self.queName} A: {self.readA()} B:{self.readB()} p: {self.pos}")
            
    def readA(self):
        if self.Ataken:
            return None
        else:
            return self.A
       
    def readB(self):
        if self.Btaken:
            return None
        else:
            return self.B
    
    def readH(self):
        if self.Htaken:
            return None
        else:
            return self.H
    
    
    def setHold(self, v):
        if v is None:
            print(f"Slot Hold remains open")
            self.Htaken = True
        elif self.Htaken is True:
            self.H = v
            self.Htaken = False
        else:
            raise Exception("Information loss - trying to overwrite Hold")
    def getHold(self):
        if self.Htaken:
            return None
        else:
            self.Htaken = True
            return self.H
    def getA(self):
        if self.Ataken:
            return None
        else:
            self.Ataken = True
            return self.A.copy()
        
    def getB(self):
       if self.Btaken:
            return None
       else:
            self.Btaken =True
            return self.B.copy()
    def setB(self, v):
        if v is None:
            print(f"Slot B remains open")
            self.Btaken = True
        elif self.Btaken is True:
            if v is None:
                print(f"Problem : Trying to set b value as none")
                raise Exception("Should never try to set value to be None")
            self.B = v
            self.Btaken = False
        else:
            raise Exception(f"Information loss - trying to overwrite B {self.queName} at position {self.pos}")
    def setA(self, v):
        if v is None:
            print(f"Slot A remains open")
            self.Ataken = True
        elif self.Ataken is True:
            if v is None:
                print(f"Problem : Trying to set b value as none")
                raise Exception("Should never try to set value to be None")
            self.A = v
            self.Ataken = False
        else:
            raise Exception("Information loss - trying to overwrite A")
  
    def passItOn2(self, u,  b = None):
        
        x, y, z = self.pos + self.A
        uNextQue = u[x][y][z]
        #print(f"Calling passItOn on {self.queName} pos {self.pos}")
        if uNextQue is None:
            
            print(f"{self.queName} - Next is none, need que. Writing {b} to hold")
            self.setHold(b)
        else:
            if uNextQue is self :
                raise Exception("Linking error - unext is same as self")
            if uNextQue is not self.nextQue:
                print(f"Collission or linking error")
                if self.nextQue is None:
                    print(f"The problem is that the nextQue is not set. UnextQue is {uNextQue.queName}")
                elif uNextQue is None:
                    print(f"The problem is that the next thing in the universe is None")
                else:
                    print(f"Universe says next que is {uNextQue.queName}, but linking says next Que is {self.nextQue.queName}")
                exit()
            print(f"Passing from {self.queName} to {uNextQue.queName}")
            uNextQue.passItOn2(u, self.getB())
            print(f"Now setting B on {self.queName} {b}")
            self.setB(b)
    
          
    def secondPass2(self, u, t, recycleThisQue  = None,  b = None):
        #u is universe
        # t is current time
        # the gobbled que to be recycled
        # b the queue being pushed forward
        x, y, z = self.pos + self.A
        uNextQue = u[x][y][z]
        if recycleThisQue is None:
            #This is the initial call
            recycleThisQue = self
            uNextQue.secondPass2(u, t, self, self.getA())
       
        elif uNextQue is None:
            
            print("EVENT")
            print(f"Determines where the next que is laid {self.queName}")
            print(f"Laying the recylced {recycleThisQue.queName}")
            print(f"Current position of {self.queName} is {self.pos}, a is {self.A}")
            print(f"Next position of {recycleThisQue.queName} is { [x, y, z] }")
            xl, yl, zl = recycleThisQue.pos #The old coordinate for this que
            u[xl][yl][zl]  = None # Remove the reference to this que from universe at old position
            u[x][y][z] = recycleThisQue #Create new reference in u
            recycleThisQue.pos = [x, y, z]
            recycleThisQue.setA(self.getHold()) # Apply the values from memory into recycled que
            print(f"Setting b at {recycleThisQue.queName} at pos {recycleThisQue.pos}")
            recycleThisQue.setB(b) # Allow the queue to push into the new lay
            recycleThisQue.nextQue = None
            self.nextQue = recycleThisQue # Linked list maintenance
            recycleThisQue.prevQue = self # Linked list maintenance
            recycleThisQue.x[t] = [x, y, z]
        elif uNextQue is not None:
           print(f"We want to push from {self.queName} B:{self.readB()} into {self.nextQue.queName} setting {self.readB()}")
           if uNextQue is self:
               raise Exception("Linking error")
           self.nextQue.secondPass2(u, t, recycleThisQue, self.getB()) # TEmpties the B and passes it on
           self.x[t] = self.pos.copy()
           print(f"We want tor write B {b} to {self.queName} at pos {self.pos}")
           self.setB(b)  # notice that getB() is used in both branches of if, the B is always emptied
        print(f"{self.queName} - {self.readA()} , {self.readB()}")
        

    
    def passItOn(self, u, b = None):
        x, y, z = self.pos + self.A
        uNextQue = u[x][y][z]
        
        if uNextQue is None:
            print(f"{self.queName} finds no que at {[x, y, z]}, lays one:")
            que = Que(self, self.getB(), None, u, None, "newQue")
            self.setB(b)
           
            self.nextQue = que
        elif uNextQue is not self.nextQue:
            print(f"We have {uNextQue.queName} vs {self.nextQue.queName}")
            self.nextQue = uNextQue #                   Information is lost
            print(f"We have exchange. Will use the next value")
            self.nextQue.passItOn(u, self.getB())
            self.setB(b)
        else:
            self.nextQue.passItOn(u, self.getB())
            self.setB(b)
            
         
    def secondPass(self, u, recycleThisQue,  b = None):
       
        
        x, y, z = self.pos + self.A
        uNextQue = u[x][y][z]
       
        if uNextQue is not None:
           self.nextQue.secondPass(u, recycleThisQue, self.getB())
        print("is this right before?")
        self.setB(b)
        print(f"{self.queName} - {self.readA()} , {self.readB()}")
        

u =  [[[None for k in range(n)] for j in range(n)] for i in range(n)]

def getOrientationStrByVector(ns):
    if ns is None:
        return "none"
    if ns[0] == 1  and  ns[1] == 0  and  ns[2] == 0:
        return "right"
    if ns[0] == -1  and  ns[1] == 0  and  ns[2] == 0:
        return "left"
    if ns[0] == 0  and  ns[1] == 1  and  ns[2] == 0:
        return "up"
    if ns[0] == 0  and  ns[1] == 0  and  ns[2] == 1 :
        return "in"
    if ns[0] == 0  and  ns[1] == 0  and  ns[2] == -1 :
        return "out"      
      

def getOrientationByInt(ns):
    if ns == 0:
        return np.array([1, 0, 0])
    elif ns == 1:
        return  np.array([-1, 0, 0])
    elif ns == 2:
        return  np.array([0, 1, 0])
    elif ns == 3:
        return  np.array([0, -1, 0])
    elif ns == 4:
        return np.array([0, 0, 1])
    elif ns == 5:
        return np.array([0, 0, -1])
        
def getOrientationByName(ns):
    if ns == "right":
        return np.array([1, 0, 0])
    elif ns == "left":
        return  np.array([-1, 0, 0])
    elif ns == "up":
        return  np.array([0, 1, 0])
    elif ns == "down":
        return  np.array([0, -1, 0])
    elif ns == "in":
        return np.array([0, 0, 1])
    elif ns == "out":
        return np.array([0, 0, -1])
    


def setPath(avgSegmentLength, segments):
	
	startPos = [3, 3, 3]
	currPos = startPos.copy()
	segSymb = np.zeros(segments)
	segLen = int(np.zeros(np.random()*avgSegmentLength)+0.5)
	segSym += np.randInt(0, 6) 
	
	for i in range(0, segments):
		orient = getOrientation(segSymb[i])
		for j in range(0, segLen[i]):
			currPos += orient
			x, y, z = currpos
			if u[x][y][z] is not None:
				print(f"Self collission at {x}, {y}, {z}. ")
				currPos -= orient
				break
			u[x][y][z] = orient.copy()


class Gobbler:
    def __init__(self, que, T):
        self.currentPos = que.pos
        self.currentQue = que
        self.x = np.array([[0., 0., 0.] for i in range(T)]) # Position
        self.T = T
        self.x[0] = self.currentPos.copy()
        self.memory = np.array([0, 0, 0])
        print(f"Gobbler kommer! {self.currentPos}")
    	
    def advance(self, u, t):
        print(f"CAlling advance on Gobbler with t:{t}")
       
        x, y, z = self.currentPos
        self.currentQue = u[x][y][z]
        xn, yn, zn = self.currentPos + self.currentQue.readA()
        
        if self.currentQue is None:
            raise Exception("Gobbler ate into something that wasnt a que")
         
        print(f"Reading universe {self.currentPos}  is {u[x][y][z].queName}")
        
        #self.memory = copy.copy(u[x][y][z][0])
        #self.brickLayer(u)
        self.currentQue.passItOn2(u)
        
        print(f"GOBBLER SQUATS ON : {self.currentQue.queName}")
        self.currentQue.printRecursive()
        #|a = self.currentQue.getA()
        self.currentQue.secondPass2(u, t)
        if self.currentQue is None:
            raise Exception("Thadslkdj")
        
        
        #u[x][y][z] = None 
        print(f"--Gobbler moving from {self.currentQue.queName} at {self.currentPos} to {self.currentPos + self.currentQue.A}")
        
        if u[xn][yn][zn] is None:
            raise Exception(f"Expected que at {[x, y, z]} was {u[x][y][z]}. Gobbler moving into something that isnt a que")
        self.currentPos = [xn, yn, zn]
        self.currentQue  = u[xn][yn][zn] 
        self.x[t] = self.currentPos.copy()
        self.currentQue.printRecursive()

        
         
        
       
        
        

    # def brickLayer(self, u ):
    #     x, y, z = self.currentPos
    #                                                                               # Bricklayer, where do you begin now?
    #     xi, yi, zi = self.currentPos    
        
    #     cellCoord = np.array([xi, yi, zi])      
    #     #Start with two empty hands
    #     print(f"Begin bricklayer - CellCoord is {cellCoord} ")
    #     rightHand = None
    #     #1 Pickup the tile we are leaving in right hand
    #     rightHand = u[xi][yi][zi][0]
    #                                            #The tile before the  
    #     if rightHand is None :
    #         raise Exception("Impossibility : rightand cannot be none")
    #     while u[xi][yi][zi][0] is not None: 
            
          
    #         print(f"Walking to {cellCoord}")
    #         print(f"Contains {u[xi][yi][zi]}")
                                                                        
    #         leftHand = copy.copy(u[xi][yi][zi][1])                              # Take the next boot in left hand
    #         u[xi][yi][zi][1] = copy.copy(rightHand)                                       # Drop the old boot
    #         rightHand = copy.copy(leftHand)  
    #                                               # Move the boot from the left to right hand
    #         leftHand = None
            
    #         cellCoord += u[xi][yi][zi][0]                                      # Follow the arrow
    #         xi, yi, zi = cellCoord
        
    #     print(f"Laying brick on position {cellCoord} will lay {rightHand}")
    #     u[xi][yi][zi][0] = rightHand # Empty your hands, builder                   # The space ahead of you is empty. Empty your right hand into it
    #     print(f"Resetting the gobbled at {[x, y, z]}")
    #     u[x][y][z] = [None, None]  

    
    
    
    def plotPosition(self):
        t = np.linspace(0, self.T, self.T)
        
        xPosition = self.x[:,0]
        yPosition = self.x[:,1]
        zPosition = self.x[:,2]
        print(yPosition)
        
        fig, (ax1, ax2, ax3)  = plt.subplots(3)
        fig.suptitle('XYZ position over time')
        ax1.plot(t, xPosition)
#        ax1.setTitle("X")
        ax2.plot(t, yPosition)
 #       ax2.setTitle("X")
        ax3.plot(t, zPosition)
  #      ax3.setTitle("Z")
          
        plt.title("Gobbler position in time")
           
        plt.legend()
        plt.show()
    			


def test1():
    T= 20
    x0, y0, z0 = [15, 15, 15]
    
    que1 = Que(None, getOrientationByName("up"), getOrientationByName("up"), u, T, [x0, y0, z0], "InitialQue")
    que2 = Que(que1, getOrientationByName("up"), getOrientationByName("up"), u, T,  None, "seconQue")
    que3 = Que(que2, getOrientationByName("up"), getOrientationByName("up"), u, T,  None, "thirdque")
    que4 = Que(que3, getOrientationByName("up"), getOrientationByName("up"), u, T,  None, "fourthque")
    que5 = Que(que4, getOrientationByName("up"), getOrientationByName("up"), u, T,  None, "fifthque")
    
    queList =[que1, que2]

    assert getOrientationStrByVector(que1.A) == "up", f"Que 1 A should be [0, 1, 0] was {que1.A}"
     
    assert getOrientationStrByVector(que1.B) == "up", f"Que 1 B should be [0, 1, 0] was {que1.B}"
   
    assert getOrientationStrByVector(que2.A) == "up", f"Que 1 A should be [0, 1, 0] was {que1.A}"
    
    assert getOrientationStrByVector(que2.B) == "up", f"Que 1 B should be [0, 1, 0] was {que1.B}"

    que1.passItOn2(u, que1.getA())

   
    assert getOrientationStrByVector(que1.readA()) == "none", f"Que 1 A should be none was {que1.A}"
    assert getOrientationStrByVector(que1.readB()) == "up", f"Que 1 B should be up was {que1.B}"
   
    
    que1.nextQue.secondPass2(u, 1, que1, que1.getB())
    queRec = que5.nextQue
    
    assert que1 is queRec, f"Que 3 and 1 should be the same object"
    assert getOrientationStrByVector(que1.readA()) == "up", f"Que 1 A should be [0, 1, 0] was {que1.readA()}"
    assert getOrientationStrByVector(que1.readB()) == "up", f"Que 1 B should be [0, 1, 0] was {que1.readB()}"
   
    assert getOrientationStrByVector(que3.readB()) == "up", f"Que 3 B should up none was {que3.readB()}"




def test2():
    T= 20
    x0, y0, z0 = [15, 15, 15]
    
    que1 = Que(None, getOrientationByName("right"), getOrientationByName("up"), u, T, [x0, y0, z0], "InitialQue")
    que2 = Que(que1, getOrientationByName("right"), getOrientationByName("up"), u, T,  None, "seconQue")
    que3 = Que(que2, getOrientationByName("right"), getOrientationByName("in"), u, T,  None, "thirdque")
    que4 = Que(que3, getOrientationByName("up"), getOrientationByName("up"), u, T,  None, "fourthque")
    que5 = Que(que4, getOrientationByName("left"), getOrientationByName("up"), u, T,  None, "fifthque")
    
   
    queList =[que1, que2, que3, que4, que5]
    
    gobbler = Gobbler(que1, T)
    print(f" Starting simulation ")
    for i in range(1, T):
        print(f"Frame {i}")
        gobbler.advance(u, i)

    gobbler.plotPosition()
    


def test3():
  
    N = 3
    T= 1
    x0, y0, z0 = [15, 15, 15]
    queList = []
    path = [[5, "up"]]
    
    initialQue = Que(None, getOrientationByName("up"), getOrientationByName("up"), u, [x0, y0, z0])
    
    queList.append(initialQue)
    prevQue = initialQue
    for i in range(0, N):
        nQue = Que(prevQue, getOrientationByName("up"), getOrientationByName("up"), u)
        queList.append(nQue)
        prevQue = nQue
        
    gobbler = Gobbler(queList[0], T)
   
    for i in range(0, T-1):
        gobbler.advance(u)
     
        


    gobbler.plotPosition()
    
    
test2()
	
	


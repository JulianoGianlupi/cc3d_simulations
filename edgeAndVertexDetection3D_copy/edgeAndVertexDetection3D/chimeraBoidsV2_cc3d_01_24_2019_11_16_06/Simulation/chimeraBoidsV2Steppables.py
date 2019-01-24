
from PySteppables import *
import CompuCell
import sys

import numpy as np
import os


class chimeraBoidsV2Steppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = 64.
        self.lambdaVolume = 8.
        self.forceModulus = -20
        self.forceCounter = 0
        self.deltaTime = 10
        
        #boids parameters
        self.alphaBoids = 1.5
        self.betaBoids = 5.
        self.gammaBoids = .5
        self.deltaBoids = .1
        
        #assigning parameters
        for cell in self.cellList:
            if cell:
                cell.targetVolume = self.targetVolume
                cell.lambdaVolume = self.lambdaVolume
                
                
                cell.dict['forceAngle'] = np.random.uniform(-np.pi,np.pi)
                cell.dict['angle']=cell.dict['forceAngle']
                
                cell.dict['centerMassX'] = [cell.xCOM]
                cell.dict['centerMassY'] = [cell.yCOM]
                cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
                
                cell.dict['positionX'] = [cell.xCOM]
                cell.dict['positionY'] = [cell.yCOM]
                
                cell.dict['velocityX'] = 0.
                cell.dict['velocityY'] = 0.
                
                cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
    

    def positionTracking(self,mcs,cur_Cell):
        self.positionX = cur_Cell.dict['positionX']
        self.positionY = cur_Cell.dict['positionY']
        
        dx = cur_Cell.xCOM - cur_Cell.dict['centerMassX']
        dy = cur_Cell.yCOM - cur_Cell.dict['centerMassY'] 
        
        #if the cell "moved" more than half the space it whent trough the border
        if abs(dx)>0.5*self.dim.x:
            if dx < 0:
                dx +=self.dim.x
            else:
                dx -= self.dim.x
        if abs(dy)>0.5*self.dim.y:
            if dy<0:
                dy +=self.dim.y
            else:
                dy -=self.dim.y  
        
        #new absolute position is previous position + fixed dx
        cur_Cell.dict['positionX'].append(self.positionX[-1]+dx)
        cur_Cell.dict['positionY'].append(self.positionY[-1]+dy)
        
        #with a dt = 1mcs v = dx
        cur_Cell.dict['velocityX'] = dx
        cur_Cell.dict['velocityY'] = dy
        
        #angle of movement:
        cur_Cell.dict['angle'] = np.arctan2(dy,dx)
        
        cur_Cell.dict['centerMassX'] = cell.xCOM
        cur_Cell.dict['centerMassY'] = cell.yCOM
    
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        for cell in self.cellList:
            if cell:
                positionTracking(self,mcs,cell)
                print cell.dict['positionX']
#                 #tracking position trough border
#                 self.positionX = cell.dict['positionX']
#                 self.positionY = cell.dict['positionY']
                
#                 dx = cell.xCOM - cell.dict['centerMassX']
#                 dy = cell.yCOM - cell.dict['centerMassY'] 
                
#                 #if the cell "moved" more than half the space it whent trough the border
#                 if abs(dx)>0.5*self.dim.x:
#                     if dx < 0:
#                         dx +=self.dim.x
#                     else:
#                         dx -= self.dim.x
#                 if abs(dy)>0.5*self.dim.y:
#                     if dy<0:
#                         dy +=self.dim.y
#                     else:
#                         dy -=self.dim.y  
                
#                 #new absolute position is previous position + fixed dx
#                 cell.dict['positionX'].append(self.positionX[-1]+dx)
#                 cell.dict['positionY'].append(self.positionY[-1]+dy)
                
#                 #with a dt = 1mcs v = dx
#                 cell.dict['velocityX'] = dx
#                 cell.dict['velocityY'] = dy
                
#                 #angle of movement:
#                 cell.dict['angle'] = np.arctan2(dy,dx)
                
#                 cell.dict['centerMassX'] = cell.xCOM
#                 cell.dict['centerMassY'] = cell.yCOM
                
                
                
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        